#bring in LANDIS layers to generate new prediction surfaces

#proxies from LANDIS to initial inputs:
# biomass: biomass
# height: 98th percentile of cohort heights, estimated from regressions from FIA relating age ~ height
# fhd: shannon diversity of cohort heights
# understory ratio: proportion of cohort biomass that is below 5m, compared to total biomass. Calculated
#                   from heights and assumptions about 
# area_short: classify sites as "short" or not, and then take a moving average window approach
# prop_decid (and similar): classify sites as deciduous, then take a moving average window approach

library("terra")
library("sf")
library("tidyverse")
options(warn = 0)
input_dir <- "D:/SApps LANDIS/"

comm_output <- read.csv("./landis_analysis/landis_predictor_layers/initial landscape layers/community-input-file-0.csv")
comm_map <- terra::rast("./landis_analysis/landis_predictor_layers/initial landscape layers/output-community-0.img")
plot(comm_map)

climate_future <- read.csv(paste0("./landis_analysis/landis_predictor_layers/initial landscape layers/Climate-future-input-log.csv"))
ecoregions <- terra::rast(paste0(input_dir, "/Inputs/Basic_inputs/Ecos11_NCLD.tif"))
plot(ecoregions)

predictor_stack <- terra::rast("./landis_analysis/landis_predictor_layers/predictor_stack_bcr28.tif")
predictor_stack2 <- predictor_stack %>%
  project(ecoregions, mask = TRUE, method = "near", align = FALSE) %>%
  mask(ecoregions, maskvalues = 1)
predictor_stack_coarse <- terra::aggregate(predictor_stack2, 5, na.rm = TRUE)
samp <- sample(length(predictor_stack2[]), 10000)
samp_coarse <- sample(length(predictor_stack_coarse[]), 10000)

comm_rast <- terra::ifel(comm_map %in% comm_output$MapCode, comm_map, NA)


#----------------------------------------------------------
# estimate cohort attributes from FIA
library("tidyverse")
library("rFIA")


#what states should we use for the analysis?
states <- c("NC", "TN", "VA", "KY", "GA", "SC")

#this script uses the rFIA package to access the tables needed
tables <- c("TREE","SEEDLING","PLOT", "COND", "SITETREE")
directory <- "D:/Data/fia/rFIA_downloads"


#import fia data
#using rFIA package automatically knits the tables together; you could also
# use readr::read_csv() or import several csvs then rbind() 
fia <- readFIA(dir = directory,
               tables = tables,
               states = states)

trees <- fia$TREE
plot <- fia$PLOT %>%
  mutate(PLOT.YEAR = paste(CN, INVYR, sep="."))
cond <- fia$COND
seedlings <- fia$SEEDLING
sitetrees <- fia$SITETREE

rm(fia)
gc()

sp_ref <- read.csv("D:/Data/fia/FIADB_REFERENCE/REF_SPECIES.csv")


cond_to_use <- cond %>%
  filter(!(DSTRBCD1 %in% c(30,31,32,46,53,54,80)),
         !(DSTRBCD2 %in% c(30,31,32,46,53,54,80)),
         !(DSTRBCD3 %in% c(30,31,32,46,53,54,80)),
         TRTCD1 == 0 | is.na(TRTCD1),
         TRTCD2 == 0 | is.na(TRTCD2),
         TRTCD3 == 0 | is.na(TRTCD3)) %>%
  mutate(IS_FOREST = ifelse(FORTYPCD %in%(c(1:998)), 1, 0)) %>%
  group_by(PLT_CN) %>%
  summarise(total_cond = sum(CONDPROP_UNADJ),
            natural = sum(STDORGCD, na.rm = TRUE),
            treatment = sum(TRTCD1, na.rm = TRUE),
            proportion_forest = sum(CONDPROP_UNADJ * IS_FOREST)) %>%
  filter(total_cond > 0.95,
         proportion_forest > 0.95)

plots_to_use <- plot %>%
  filter(PLOT_STATUS_CD == 1) %>%
  left_join(cond_to_use, by = c("CN" = "PLT_CN")) %>%
  dplyr::select(CN, proportion_forest)


#species reference data
sp_data <- read.csv(paste0(input_dir, "/Inputs/NECN/NECN_Spp_Table-NECN_v6.csv"))
spp_to_use <- sp_data$SpeciesCode

# func_table <- read.csv("./Models/Input_files/necn_functional_params.csv")
# sp_table <- read.csv("./Models/Input_files/necn_species_params.csv") %>%
#   left_join(func_table, by = "FunctionalGroupIndex") 
spp_to_use_all <- sp_data$SpeciesCode

sp_ref$SpeciesCode <- paste0(substr(sp_ref$GENUS, 1, 4), substr(sp_ref$SPECIES, 1, 4) %>%
                               stringr::str_to_title())

spp_to_use_all[!(spp_to_use_all %in% sp_ref$SpeciesCode)]

spp_crosswalk <- sp_ref[sp_ref$SpeciesCode %in% spp_to_use_all, ] %>%
  dplyr::arrange(SpeciesCode) %>%
  dplyr::select(SpeciesCode, SPCD) %>%
  dplyr::filter(!(SPCD %in% c(318, 114, 8349, 845, 953, 952)))
spp_crosswalk <- rbind(spp_crosswalk,
                       c("AescBuck", 330),
                       c("CaryCodi", 402),
                       c("FrasFirr", 16),
                       c("PinuEnch", 110),
                       c("PlanOcid", 731),
                       c("PrunPenn", 761),
                       c("SassAlid", 931),
                       c("TiliAmhe", 951)) %>%
  mutate(SPCD = as.integer(SPCD))

sitetrees_bak <- sitetrees

sitetrees <- sitetrees_bak %>%
  # filter(PLT_CN %in% plots$CN) %>%
  mutate(SFTWD_HRDWD = sp_ref[match(SPCD, sp_ref$SPCD), "SFTWD_HRDWD"]) %>% #add a softwood/hardwood column for later 
  dplyr::rename(CohortAge = AGEDIA) %>%
  mutate(CohortAge = CohortAge + 5)

#fit a linear regression 
tree_regressions <- sitetrees %>% 
  dplyr::filter(!is.na(HT) & !is.na(CohortAge) & !is.na(SPCD)) %>%
  dplyr::filter(SPCD %in% spp_crosswalk$SPCD) %>%
  dplyr::group_by(SPCD) %>%
  dplyr::do(model = lm(log(HT) ~ log(CohortAge), data = .))

newdat <- data.frame(CohortAge = seq(1, 200, length.out = 1000))

map(tree_regressions$model, .f = function(x) summary(x))
map(tree_regressions$model, .f = function(x) plot(exp(predict(x, newdata = newdat)) ~ newdat$CohortAge))


#do models for hardwood and softwood -- 
hardwood_regression <- lm(HT ~ log(CohortAge) + 0, data = sitetrees[sitetrees$SFTWD_HRDWD == "H", ])
softwood_regression <- lm(HT ~ log(CohortAge) + 0, data = sitetrees[sitetrees$SFTWD_HRDWD == "S", ])

tree_regressions <- tibble(SPCD = c(1,2),
                           model = list(hardwood_regression, softwood_regression)) %>%
  bind_rows(tree_regressions, .)


comm_height <- comm_output %>%
  left_join(select(spp_crosswalk, SpeciesCode, SPCD), by = c("SpeciesName" = "SpeciesCode")) %>%
  left_join(select(sp_ref, SPCD, SFTWD_HRDWD)) %>%
  mutate(HT = NA)
  
  
spcd_in_comm <- unique(comm_height$SPCD)

for(i in 1:length(spcd_in_comm)){
  
  spcd = spcd_in_comm[i]
  newdat <- comm_height[which(comm_height$SPCD == spcd), ] %>% select(CohortAge)
  
  if(spcd %in% tree_regressions$SPCD){
    
    sp <- spp_crosswalk$SpeciesCode[spp_crosswalk$SPCD == spcd]
    
    comm_height[which(comm_height$SPCD == spcd), "HT"] <- exp(predict(tree_regressions$model[tree_regressions$SPCD == spcd][[1]],
                                                      newdata = newdat))
  
  }else if(spcd %in% sp_ref$SPCD[sp_ref$SFTWD_HRDWD == "H"]){
      comm_height[which(comm_height$SPCD == spcd), "HT"] <- predict(hardwood_regression,
                                                                    newdata = newdat)
  }else{
    comm_height[which(comm_height$SPCD == spcd), "HT"] <- predict(softwood_regression,
                                                                  newdata = newdat)
  }
}

comm_height$HT <- comm_height$HT / 3.28

height_rast_max <- comm_height %>%
  group_by(MapCode) %>%
  summarise(height = quantile(HT, 1))
height_rast_max <-  terra::subst(x = comm_map, 
                                from = height_rast_max$MapCode, 
                                to = height_rast_max$height) #convert to m
plot(height_rast_max)
plot(predictor_stack2$height[][samp] ~ height_rast_max[][samp])
abline(0,1)
summary(lm(predictor_stack2$height[][samp] ~ height_rast_max[][samp]))
abline(coef(lm(predictor_stack2$height[][samp] ~ height_rast_max[][samp])))
height_rast_max_coarse <- terra::aggregate(height_rast_max, 5, na.rm = TRUE)
plot(predictor_stack_coarse$height[][samp_coarse] ~ height_rast_max_coarse[][samp_coarse])
abline(0,1)
summary(lm(predictor_stack_coarse$height[][samp_coarse] ~ height_rast_max_coarse[][samp_coarse]))
abline(coef(lm(predictor_stack_coarse$height[][samp_coarse] ~ height_rast_max_coarse[][samp_coarse])))


height_rast_98 <- comm_height %>%
  group_by(MapCode) %>%
  summarise(height = quantile(HT, 0.98))
height_rast_98 <-  terra::subst(x = comm_map, 
                                from = height_rast_98$MapCode, 
                                to = height_rast_98$height) #convert to m
plot(height_rast_98)
plot(predictor_stack2$height[][samp] ~ height_rast_98[][samp])
abline(0,1)
summary(lm(predictor_stack2$height[][samp] ~ height_rast_98[][samp]))
abline(coef(lm(predictor_stack2$height[][samp] ~ height_rast_98[][samp])))
height_rast_98_coarse <- terra::aggregate(height_rast_98, 5, na.rm = TRUE)
plot(predictor_stack_coarse$height[][samp_coarse] ~ height_rast_98_coarse[][samp_coarse])
abline(0,1)
summary(lm(predictor_stack_coarse$height[][samp_coarse] ~ height_rast_98_coarse[][samp_coarse]))
abline(coef(lm(predictor_stack_coarse$height[][samp_coarse] ~ height_rast_98_coarse[][samp_coarse])))


height_rast_95 <- comm_height %>%
  group_by(MapCode) %>%
  summarise(height = quantile(HT, 0.95))
height_rast_95 <-  terra::subst(x = comm_map, 
                                from = height_rast_95$MapCode, 
                                to = height_rast_95$height) #convert to m
plot(height_rast_95)
plot(predictor_stack2$height[][samp] ~ height_rast_95[][samp])
abline(0,1)
summary(lm(predictor_stack2$height[][samp] ~ height_rast_95[][samp]))
abline(coef(lm(predictor_stack2$height[][samp] ~ height_rast_95[][samp])))
height_rast_95_coarse <- terra::aggregate(height_rast_95, 5, na.rm = TRUE)
plot(predictor_stack_coarse$height[][samp_coarse] ~ height_rast_95_coarse[][samp_coarse])
abline(0,1)
summary(lm(predictor_stack_coarse$height[][samp_coarse] ~ height_rast_95_coarse[][samp_coarse]))
abline(coef(lm(predictor_stack_coarse$height[][samp_coarse] ~ height_rast_95_coarse[][samp_coarse])))



height_rast_90 <- comm_height %>%
  group_by(MapCode) %>%
  summarise(height = quantile(HT, 0.90))
height_rast_90 <-  terra::subst(x = comm_map, 
                                from = height_rast_90$MapCode, 
                                to = height_rast_90$height)  #convert to m

plot(height_rast_90)
plot(predictor_stack2$height[][samp] ~ height_rast_90[][samp])
abline(0,1)
summary(lm(predictor_stack2$height[][samp] ~ height_rast_90[][samp]))
abline(coef(lm(predictor_stack2$height[][samp] ~ height_rast_90[][samp])))

height_rast_50 <- comm_height %>%
  group_by(MapCode) %>%
  summarise(height = quantile(HT, 0.50))

height_rast_50 <-  terra::subst(x = comm_map, 
                                from = height_rast_50$MapCode, 
                                to = height_rast_50$height)  #convert to m

plot(height_rast_50)
plot(predictor_stack2$height[][samp] ~ height_rast_50[][samp])
abline(0,1)
summary(lm(predictor_stack2$height[][samp] ~ height_rast_50[][samp]))
abline(coef(lm(predictor_stack2$height[][samp] ~ height_rast_50[][samp])))


#---------------------------------------------------------
# understory_ratio 
# uses height data from above

under_ratio <- comm_height %>%
  mutate(bin = cut(HT, breaks = seq(0, 100, 5), labels = FALSE)) %>%
  group_by(MapCode, bin) %>%
  mutate(bio_height = sum(CohortBiomass)) %>%
  group_by(MapCode, bin) %>%
  summarise(biomass_sum_bin = sum(bio_height)) %>%
  ungroup() %>%
  tidyr::complete(MapCode, bin, fill = list(biomass_sum_bin = 0)) %>%
  group_by(MapCode) %>%
  summarise(understory_2 = biomass_sum_bin[bin <= 2]/sum(biomass_sum_bin),
            understory_3 = sum(biomass_sum_bin[bin <= 3])/sum(biomass_sum_bin),
            understory_4 = sum(biomass_sum_bin[bin <= 4])/sum(biomass_sum_bin))
  
# under_rast <- terra::ifel(comm_map %in% under_ratio$MapCode, comm_map, NA)
under_ratio_10m <- comm_height %>%
  mutate(prop_under_10m = 10 / HT) %>% 
  mutate(understory = prop_under_10m) %>%
  mutate(understory = ifelse(SpeciesName %in% c("CornFlor", "AmelArbo", "AcerPens",
                                                "HaleDipt", "OxydArbo", "SassAlid",
                                                "PrunPenn"),
                             1, understory)) %>%
  mutate(understory_biomass = understory * CohortBiomass) %>%
  group_by(MapCode) %>%
  summarise(understory_ratio = sum(understory_biomass) / sum(CohortBiomass))

# under_rast <- terra::ifel(comm_map %in% under_ratio$MapCode, comm_map, NA)
under_rast_4 <-  terra::subst(x = comm_map,
                            from = under_ratio$MapCode,
                            to = under_ratio$understory_4)
plot(under_rast_4)
plot(I(under_rast_4[][samp]) ~ predictor_stack2$understory[][samp])
abline(0,1)
summary(lm(under_rast_4[][samp] ~ predictor_stack2$understory[][samp]))

under_rast_coarse <- terra::aggregate(under_rast_4, 5, na.rm = TRUE)
plot(predictor_stack_coarse$understory[][samp_coarse] ~ under_rast_coarse[][samp_coarse])
abline(0,1)
summary(lm(predictor_stack_coarse$understory[][samp_coarse] ~ under_rast_coarse[][samp_coarse]))
abline(coef(lm(predictor_stack_coarse$understory[][samp_coarse] ~ under_rast_coarse[][samp_coarse])))


# under_rast <- terra::ifel(comm_map %in% under_ratio$MapCode, comm_map, NA)
under_rast_10 <-  terra::subst(x = comm_map,
                              from = under_ratio_10m$MapCode,
                              to = under_ratio_10m$understory_ratio)
plot(under_rast_10)
samp <- sample(1:length(under_rast_10[]), 10000)
plot(I(under_rast_10[][samp]) ~ predictor_stack2$understory[][samp])
abline(0,1)
summary(lm(under_rast_10[][samp] ~ predictor_stack2$understory[][samp]))


#mean_age
mean_age <- comm_output %>%
  group_by(MapCode) %>%
  summarise(mean_age = mean(CohortAge))
age_rast <-  terra::subst(x = comm_map,
                          from = mean_age$MapCode,
                          to = mean_age$mean_age)
plot(age_rast)
plot(predictor_stack2$understory[][samp] ~ age_rast[][samp])
summary(lm(predictor_stack2$understory[][samp] ~ age_rast[][samp] + under_rast_10[][samp]))

# age_coarse <- 
# 
# # saveRDS(under_rast, "understory_raster_60.RDS")

#---------------------------------------------------------


#---------------------------------------------------------
#foliage height diversity
#---------------------------------------------------------
fhd <- comm_height %>%
  mutate(bin = cut(HT, breaks = seq(0, 100, 5), labels = FALSE)) %>%
  group_by(MapCode, bin) %>%
  summarize(bio_height = sum(CohortBiomass)) %>%
  group_by(MapCode) %>%
  mutate(shan = -1*(bio_height/sum(bio_height)) * log(bio_height/sum(bio_height))) %>%
  summarize(ht_diversity = sum(shan))
fhd_rast <- comm_map
fhd_rast <-  terra::subst(x = fhd_rast, 
                          from = fhd$MapCode, 
                          to = fhd$ht_diversity)


plot(predictor_stack2$fhd[][samp] ~ fhd_rast[][samp])
summary(lm(predictor_stack2$fhd[][samp] ~ fhd_rast[][samp]))
abline(coef(lm(predictor_stack2$fhd[][samp] ~ fhd_rast[][samp])))

fhd_coarse <- terra::aggregate(fhd_rast, 5, na.rm = TRUE)
pred_fhd_coarse <- terra::aggregate(predictor_stack2$fhd, 5, na.rm = TRUE)
plot(pred_fhd_coarse[][samp2] ~ fhd_coarse[][samp2])
summary(lm(pred_fhd_coarse[][samp2] ~ fhd_coarse[][samp2]))
abline(lm(pred_fhd_coarse[][samp2] ~ fhd_coarse[][samp2]))
abline(0,1)

#----------------------------------------------------------
#biomass
comm_output[comm_output[] == 0] <- NA

total_biomass <- comm_output %>%
  group_by(MapCode) %>%
  summarise(biomass = sum(CohortBiomass))
biomass_rast <- comm_rast
biomass_rast <- terra::subst(x = biomass_rast, 
                             from = total_biomass$MapCode, 
                             to = total_biomass$biomass) / 100 #convert to tonnes ha-1
plot(biomass_rast)
plot(predictor_stack2$biomass)
plot(predictor_stack2$biomass[][samp]~ I(biomass_rast[][samp]))
abline(0,1)
summary(lm(predictor_stack2$biomass[][samp]~ I(biomass_rast[][samp])))
abline(coef(lm(predictor_stack2$biomass[][samp]~ I(biomass_rast[][samp]))))

# saveRDS(biomass_rast, "biomass_raster_60.RDS")
lai <- terra::rast("C:/Users/Sam/Documents/Research/LANDIS Bird Habitat Change/landis_analysis/landis_predictor_layers/initial landscape layers/LAI-5.img")
plot(predictor_stack2$biomass[][samp] ~ I(lai[][samp]))
mod <- lm(predictor_stack2$biomass[][samp] ~ lai[][samp] * biomass_rast[][samp])
summary(mod)
abline(coef(lm(predictor_stack2$biomass[][samp] ~ I(lai[][samp]))))

biomass_coarse <- terra::aggregate(biomass_rast, 5, na.rm = TRUE)
pred_biomass_coarse <- terra::aggregate(predictor_stack2$biomass, 5, na.rm = TRUE)
plot(biomass_coarse)
plot(pred_biomass_coarse)
samp2 <- sample(1:length(biomass_coarse[]), 10000)
plot(biomass_coarse[][samp2] ~ pred_biomass_coarse[][samp2])
summary(lm(biomass_coarse[][samp2] ~ pred_biomass_coarse[][samp2]))
abline(0,1)

#deciduousness
#areas dominated by trees generally greater than 5 meters tall, 
#and greater than 20% of total vegetation cover. 
#More than 75% of the tree species shed foliage simultaneously in response to seasonal change.

decid_all <- comm_output
decid_all$decid = !grepl("Pinu|Tsug", comm_output$SpeciesName)
decid_sum <- decid_all %>%
  group_by(MapCode, decid) %>%
  summarise(biomass = sum(CohortBiomass))
deciduous <- decid_sum %>% 
  group_by(MapCode) %>%
  tidyr::pivot_wider(names_from = decid, values_from = biomass, values_fill = 0) %>%
  mutate(decid_prop = `TRUE`/(`TRUE`+ `FALSE`)) %>%
  mutate(decid = ifelse(decid_prop > 0.75, 1, 0)) %>%
  dplyr::select(MapCode, decid)
decid_rast <- comm_rast
decid_rast <-  terra::subst(x = decid_rast, 
                            from = deciduous$MapCode, 
                            to = deciduous$decid)
# decid_rast[is.na(decid_rast[])] <- 0
# decid_smooth <- terra::focal(decid_rast, 5, fun = mean, na.rm = TRUE)
decid_coarse <- terra::aggregate(decid_rast, 5, na.rm = TRUE)
pred_decid_coarse <- terra::aggregate(predictor_stack2$prop_decid, 5, na.rm = TRUE)
pred_hardwood_coarse <- terra::aggregate(predictor_stack2$prop_hardwood, 5, na.rm = TRUE)
plot(decid_coarse)
plot(pred_decid_coarse)
plot(decid_coarse[][samp2] ~ pred_hardwood_coarse[][samp2])


#--------------------------------------
# compare climate variables
#--------------------------------------




#climate variables
if(year > 0){
  summer_clim <- climate_future %>%
    filter(Year > 2005 + year - 10 & Year <= 2005 + year, 
           Timestep > 182 & Timestep < 273) %>%
    group_by(EcoregionIndex) %>%
    summarise(tmax = median(max_airtemp)*10,#rescale to match GEE output
              tmin = median(min_airtemp)*10,
              precip = median(ppt)*1000)
  tmax_rast <- ecoregions %>%
    terra::subst(from = summer_clim$EcoregionIndex + 2, #add 2 to line up with ecoregions; different for each model!
                 to = summer_clim$tmax,
                 others = NA)
  tmin_rast <- ecoregions %>%
    terra::subst(from = summer_clim$EcoregionIndex + 2, #add 2 to line up with ecoregions; different for each model!
                 to = summer_clim$tmin,
                 others = NA)
  precip_rast <- ecoregions %>%
    terra::subst(from = summer_clim$EcoregionIndex + 2, #add 2 to line up with ecoregions; different for each model!
                 to = summer_clim$precip,
                 others = NA)
}



#stack all the predictors
if(year == 0){
  pred_stack <- c(biomass2, 
                  under2, 
                  decid2,
                  # predictor_stack2$aet,
                  # predictor_stack2$pet,
                  # predictor_stack2$def,
                  # predictor_stack2$soil,
                  # predictor_stack2$vpd,
                  # predictor_stack2$pdsi,
                  predictor_stack2$tmmn,
                  predictor_stack2$tmmx,
                  predictor_stack2$pr,
                  predictor_stack2$slope,
                  predictor_stack2$chili,
                  predictor_stack2$tpi)
} else{
  pred_stack <- c(biomass2, 
                  under2, 
                  decid2,
                  # predictor_stack2$aet,
                  # predictor_stack2$pet,
                  # predictor_stack2$def,
                  # predictor_stack2$soil,
                  # predictor_stack2$vpd,
                  # predictor_stack2$pdsi,
                  tmin_rast,
                  tmax_rast,
                  precip_rast,
                  predictor_stack2$slope,
                  predictor_stack2$chili,
                  predictor_stack2$tpi)
}


names(pred_stack) <- c("biomass", "understory_ratio", "prop_decid",
                       "tmmn", "tmmx", "pr", "slope", "chili", "tpi")

# writeRaster(pred_stack, 
#             filename = paste0("./landis_predictor_layers/pred_stack_", 
#                               model_name, "_", year, ".tif"),
#             overwrite = TRUE)


















library("terra")
library("tidyverse")
library("sf")

template <- terra::rast("C:/Users/Sam/Documents/Research/Project-Southern-Appalachians-2018/Models/LANDIS_Sapps_Ecosystem_Papers/Ecos11_NCLD.tif")
template <- terra::classify(template, rcl = matrix(data = c(2, 12, 1),
                                             ncol = 3),
                            others = NA,
                            include.lowest = TRUE,
                            right = FALSE)
plot(template)

boundary <- sf::st_read("./maps/study_area.gpkg") %>%
  st_transform(crs(template))
# plot(landis_vars[[1]])
plot(vect(boundary), add = TRUE) 

#these are the original layers from the SDM analysis, we'll crop them to the study area
pred_stack <- rast("./landis_analysis/landis_predictor_layers/predictor_stack_bcr28.tif") %>%
  terra::project(template) %>%
  terra::crop(vect(boundary)) %>%
  terra::mask(template)

landis_vars <- rast("./landis_analysis/landis_predictor_layers/pred_stack_LowTLowV BAU_0.tif")


samp <- sample(length(pred_stack[[1]][]), 10000)
plot(pred_stack[[1]])
plot(landis_vars[[1]])
plot(pred_stack[[1]][][samp] ~ landis_vars[[1]][][samp])
abline(0,1)

plot(pred_stack[[1]] - landis_vars[[1]])
hist(pred_stack[[1]] - landis_vars[[1]])
global(pred_stack[[1]] - landis_vars[[1]], mean, na.rm = TRUE)

plot(pred_stack[[2]])
plot(landis_vars[[2]])
plot(pred_stack[[2]][][samp] ~ landis_vars[[2]][][samp])
summary(lm(pred_stack[[2]][][samp] ~ landis_vars[[2]][][samp]))
abline(0,1)

plot(pred_stack[[1]] - landis_vars[[1]])
hist(pred_stack[[1]] - landis_vars[[1]])
global(pred_stack[[1]] - landis_vars[[1]], mean, na.rm = TRUE)


