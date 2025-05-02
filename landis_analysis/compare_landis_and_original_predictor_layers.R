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

map_attr_to_comm <- function(attrib_table, comm_map, var){
  attr_map <- comm_map
  attr_map[] <- attrib_table[match(comm_map[], attrib_table$MapCode), var][[1]]
  return(attr_map)
}

map_attr_to_ecoregions <- function(attrib_table, comm_map, var){
  attr_map <- comm_map
  attr_map[] <- attrib_table[match(comm_map[], attrib_table$EcoregionIndex), var][[1]]
  return(attr_map)
}


#make a tibble to catch the models we'll be making
landis_crosswalk_mods <- tibble(var = character(),
                                mod = list())


comm_output <- read.csv("./landis_analysis/reference/initial landscape layers/community-input-file-0.csv")
comm_map <- terra::rast("./landis_analysis/reference/initial landscape layers/output-community-0.img")
plot(comm_map)

climate_future <- read.csv(paste0("./landis_analysis/reference/initial landscape layers/Climate-future-input-log.csv"))
ecoregions <- terra::rast(paste0("./landis_analysis/reference/Ecos11_NCLD.tif"))
# plot(ecoregions)

necn_monthly <- read.csv("./landis_analysis/reference/initial landscape layers/NECN-succession-monthly-log.csv")

predictor_stack <- terra::rast("./landis_analysis/landis_predictor_layers/predictor_stack_bcr28_train.tif")
predictor_stack2 <- predictor_stack %>%
  project(ecoregions, mask = TRUE, method = "near", align = FALSE) %>%
  mask(ecoregions, maskvalues = 1)
predictor_stack_coarse <- terra::aggregate(predictor_stack2, 5, na.rm = TRUE)
samp <- sample(length(predictor_stack2[[1]][]), 10000)
samp_coarse <- sample(length(predictor_stack_coarse[[1]][]), 10000)

comm_map[] <- ifelse(comm_map[] %in% comm_output$MapCode, comm_map[], NA)
comm_rast <- comm_map

#species reference data
sp_ref <- read.csv("./landis_analysis/reference/REF_SPECIES.csv")

sp_data <- read.csv("./landis_analysis/reference/NECN_Spp_Table-NECN_v6.csv")

spp_to_use <- sp_data$SpeciesCode
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


#height regressions fit in create_models_height_from_cohort_attr.R
tree_regressions <- readRDS("./landis_analysis/reference/tree_height_regressions.RDS")

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
      comm_height[which(comm_height$SPCD == spcd), "HT"] <- exp(predict(tree_regressions$model[tree_regressions$SPCD == 1][[1]], 
                                                                    newdata = newdat)) #SPCD set in create_models script; not a real FIA code
  }else{
    comm_height[which(comm_height$SPCD == spcd), "HT"] <- exp(predict(tree_regressions$model[tree_regressions$SPCD == 2][[1]],
                                                                     newdata = newdat)) #SPCD set in create_models script; not a real FIA code
  }
}

comm_height$HT <- comm_height$HT / 3.28


#----------------------------
#Compare rasters and fit crosswalk models
#---------------------------
height_rast_98 <- comm_height %>%
  group_by(MapCode) %>%
  summarise(height = quantile(HT, 0.98))
height_rast_98 <- map_attr_to_comm(height_rast_98, comm_map, "height")
plot(height_rast_98)
plot(predictor_stack2$height[][samp] ~ height_rast_98[][samp])
abline(0,1)
summary(lm(predictor_stack2$height[][samp] ~ height_rast_98[][samp]))
abline(coef(lm(predictor_stack2$height[][samp] ~ height_rast_98[][samp])))
height_rast_98_coarse <- terra::aggregate(height_rast_98, 5, na.rm = TRUE)
plot(log(predictor_stack_coarse$height[][samp_coarse]) ~ log(height_rast_98_coarse[][samp_coarse]))
abline(0,1)
summary(lm(log(predictor_stack_coarse$height[][samp_coarse]) ~ log(height_rast_98_coarse[][samp_coarse])))
abline(coef(lm(log(predictor_stack_coarse$height[][samp_coarse]) ~ log(height_rast_98_coarse[][samp_coarse]))))

pred_vals <- log(predictor_stack_coarse$height[][samp_coarse])
height_vals <- log(height_rast_98_coarse[][samp_coarse])

height_mod <- lm(pred_vals ~ height_vals)

landis_crosswalk_mods <- bind_rows(landis_crosswalk_mods,
                                   tibble(var = "height", 
                                          mod = list(height_mod)))

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
            understory_4 = sum(biomass_sum_bin[bin <= 4])/sum(biomass_sum_bin),
            understory_5 = sum(biomass_sum_bin[bin <= 5])/sum(biomass_sum_bin))
  

# under_ratio_10m <- comm_height %>%
#   mutate(prop_under_10m = 10 / HT) %>%
#   mutate(understory = prop_under_10m) %>%
#   mutate(understory = ifelse(SpeciesName %in% c("CornFlor", "AmelArbo", "AcerPens",
#                                                 "HaleDipt", "OxydArbo", "SassAlid",
#                                                 "PrunPenn"),
#                              1, understory)) %>%
#   mutate(understory_biomass = understory * CohortBiomass) %>%
#   group_by(MapCode) %>%
#   summarise(understory_ratio = sum(understory_biomass) / sum(CohortBiomass))
# under_rast_10 <- map_attr_to_comm(under_ratio_10m, comm_map, "understory_ratio")
# plot(under_rast_10)
# plot(I(under_rast_10[][samp]) ~ predictor_stack2$understory[][samp])
# abline(0,1)
# summary(lm(under_rast_10[][samp] ~ predictor_stack2$understory[][samp]))


under_rast_4 <- map_attr_to_comm(under_ratio, comm_map, "understory_4")
plot(under_rast_4)
plot(I(under_rast_4[][samp]) ~ predictor_stack2$understory[][samp])
abline(0,1)
summary(lm(under_rast_4[][samp] ~ predictor_stack2$understory[][samp]))
under_rast_coarse <- terra::aggregate(under_rast_4, 5, na.rm = TRUE)
plot(predictor_stack_coarse$understory[][samp_coarse] ~ under_rast_coarse[][samp_coarse])
abline(0,1)
summary(lm(predictor_stack_coarse$understory[][samp_coarse] ~ under_rast_coarse[][samp_coarse]))
abline(coef(lm(predictor_stack_coarse$understory[][samp_coarse] ~ under_rast_coarse[][samp_coarse])))



# #mean_age
mean_age <- comm_output %>%
  group_by(MapCode) %>%
  summarise(mean_age = mean(CohortAge))
age_rast <- map_attr_to_comm(mean_age, comm_map, "mean_age")
plot(age_rast)
plot(predictor_stack2$understory[][samp] ~ age_rast[][samp])
age_coarse <- terra::aggregate(age_rast, 5, na.rm = TRUE)

lai <- terra::rast("./landis_analysis/reference/initial landscape layers/LAI-5.img")
lai_coarse <- terra::aggregate(lai, 5, na.rm = TRUE)
lai_vals <- lai_coarse[][samp_coarse]
under_vals <- under_rast_coarse[][samp_coarse]
pred_vals <- predictor_stack_coarse$understory[][samp_coarse]
# age_vals <- age_coarse[][samp_coarse]
# biomass_vals <- biomass_coarse[][samp_coarse]

under_mod <- lm(pred_vals ~  lai_vals * under_vals )
summary(under_mod)
plot(effects::allEffects(under_mod))

landis_crosswalk_mods <- bind_rows(landis_crosswalk_mods,
                                   tibble(var = "understory", 
                                          mod = list(under_mod)))

#---------------------------------------------------------
#foliage height diversity
#---------------------------------------------------------
fhd <- comm_height %>%
  mutate(bin = cut(HT, breaks = seq(0, 100, 2), labels = FALSE)) %>%
  group_by(MapCode, bin) %>%
  summarize(bio_height = sum(log(CohortBiomass))) %>%
  group_by(MapCode) %>%
  mutate(shan = -1*(bio_height/sum(bio_height)) * log(bio_height/sum(bio_height))) %>%
  summarize(ht_diversity = sum(shan))

fhd_rast <-  map_attr_to_comm(fhd, comm_map, "ht_diversity")
plot(predictor_stack2$fhd[][samp] ~ fhd_rast[][samp])
summary(lm(predictor_stack2$fhd[][samp] ~ fhd_rast[][samp]))
abline(coef(lm(predictor_stack2$fhd[][samp] ~ fhd_rast[][samp])))

fhd_coarse <- terra::aggregate(fhd_rast, 5, na.rm = TRUE)
pred_fhd_coarse <- terra::aggregate(predictor_stack2$fhd, 5, na.rm = TRUE)
plot(pred_fhd_coarse[][samp_coarse] ~ fhd_coarse[][samp_coarse])
summary(lm(pred_fhd_coarse[][samp_coarse] ~ fhd_coarse[][samp_coarse]))
abline(lm(pred_fhd_coarse[][samp_coarse] ~ fhd_coarse[][samp_coarse]))
abline(0,1)

fhd_vals <- fhd_coarse[][samp_coarse]
pred_vals <- predictor_stack_coarse$fhd[][samp_coarse]
fhd_mod <- lm(pred_vals ~  fhd_vals)
summary(fhd_mod)
plot(effects::allEffects(fhd_mod))

landis_crosswalk_mods <- bind_rows(landis_crosswalk_mods,
                                   tibble(var = "fhd", 
                                          mod = list(fhd_mod)))


#---------------------------------------------
# area_short
#---------------------------------------------
area_short <- comm_height %>%
  group_by(MapCode) %>%
  summarize(short = as.numeric(quantile(HT, 1) < 20))
  
area_short_rast <- map_attr_to_comm(area_short, comm_map, "short")
area_short_rast2 <- terra::focal(area_short_rast, 7, fun = "mean", 
                                 na.policy = "omit", na.rm = TRUE, 
                                 expand = TRUE, fillvalue = 0)
plot(area_short_rast2)
plot(predictor_stack2$area_short)

#this is pretty good; probably should just leave it alone
mean(area_short$short, na.rm = TRUE)
mean(predictor_stack2$area_short[], na.rm = TRUE)

lai_radius <- terra::focal(lai, 5, fun = "mean")
plot(predictor_stack2$area_short[][samp] ~ lai_radius[][samp])
summary(lm(predictor_stack2$area_short[][samp] ~ lai_radius[][samp]))
abline(lm(predictor_stack2$area_short[][samp] ~ lai_radius[][samp]))

plot(predictor_stack2$area_short[][samp] ~ area_short_rast2[][samp])
summary(lm(predictor_stack2$area_short[][samp] ~ area_short_rast2[][samp] + lai_radius[][samp])) # prtty good!
abline(coef(lm(predictor_stack2$area_short[][samp] ~ area_short_rast2[][samp])))


short_coarse <- terra::aggregate(area_short_rast2, 5, na.rm = TRUE)
lai_coarse <- terra::aggregate(lai_radius, 5, na.rm = TRUE)
plot(predictor_stack_coarse$area_short[][samp_coarse] ~ short_coarse[][samp_coarse])
abline(lm(predictor_stack_coarse$area_short[][samp_coarse] ~ short_coarse[][samp_coarse]))
abline(0,1)
plot(predictor_stack_coarse$area_short[][samp_coarse] ~ lai_coarse[][samp_coarse])
summary(lm(predictor_stack_coarse$area_short[][samp_coarse] ~ short_coarse[][samp_coarse] + lai_coarse[][samp_coarse]))

mod <- lm(predictor_stack_coarse$area_short[][samp_coarse] ~ short_coarse[][samp_coarse] + lai_coarse[][samp_coarse])
plot(residuals(mod) ~ fitted(mod))
plot(mod)

area_short_vals <- short_coarse[][samp_coarse]
lai_vals <- lai_coarse[][samp_coarse]
pred_vals <- predictor_stack_coarse$area_short[][samp_coarse]
area_mod <- lm(pred_vals ~  area_short_vals + lai_vals)
summary(area_mod)
plot(effects::allEffects(area_mod))

landis_crosswalk_mods <- bind_rows(landis_crosswalk_mods,
                                   tibble(var = "area_short", 
                                          mod = list(area_mod)))
#----------------------------------------------------------
#biomass
comm_output[comm_output[] == 0] <- NA

total_biomass <- comm_output %>%
  group_by(MapCode) %>%
  summarise(biomass = sum(CohortBiomass))
biomass_rast <- comm_rast
biomass_rast <- map_attr_to_comm(total_biomass, comm_map, "biomass") / 100 #convert to tonnes ha-1
plot(biomass_rast)
plot(predictor_stack2$biomass)
plot(predictor_stack2$biomass[][samp]~ I(biomass_rast[][samp]))
abline(0,1)
summary(lm(predictor_stack2$biomass[][samp]~ I(biomass_rast[][samp])))
abline(coef(lm(predictor_stack2$biomass[][samp]~ I(biomass_rast[][samp]))))

biomass_coarse <- terra::aggregate(biomass_rast, 5, na.rm = TRUE)
pred_biomass_coarse <- terra::aggregate(predictor_stack2$biomass, 5, na.rm = TRUE)
plot(biomass_coarse)
plot(pred_biomass_coarse)
plot(biomass_coarse[][samp_coarse] ~ pred_biomass_coarse[][samp_coarse])
summary(lm(biomass_coarse[][samp_coarse] ~ log(pred_biomass_coarse[][samp_coarse])))
abline(0,1)


biomass_vals <- biomass_coarse[][samp_coarse]
pred_vals <- predictor_stack_coarse$biomass[][samp_coarse]
biomass_mod <- lm(pred_vals ~  log(biomass_vals))
summary(biomass_mod)
plot(effects::allEffects(biomass_mod))

landis_crosswalk_mods <- bind_rows(landis_crosswalk_mods,
                                   tibble(var = "biomass", 
                                          mod = list(biomass_mod)))

#---------------------------------------
#proportion hardwood -- trianing data from Treemap
decid_all <- comm_output
decid_all$decid = !grepl("Pinu|Tsug|Pice|FrasFirr", comm_output$SpeciesName)
decid_sum <- decid_all %>%
  group_by(MapCode, decid) %>%
  summarise(biomass = sum(CohortBiomass))
deciduous <- decid_sum %>% 
  group_by(MapCode) %>%
  tidyr::pivot_wider(names_from = decid, values_from = biomass, values_fill = 0) %>%
  mutate(decid_prop = `TRUE`/(`TRUE`+ `FALSE`)) %>%
  mutate(decid = ifelse(decid_prop > 0.75, 1, 0)) %>%
  dplyr::select(MapCode, decid)
hardwood <- decid_sum %>%
  group_by(MapCode) %>%
  tidyr::pivot_wider(names_from = decid, values_from = biomass, values_fill = 0) %>%
  mutate(decid = `TRUE`/(`TRUE`+ `FALSE`)) %>%
  dplyr::select(MapCode, decid)
decid_rast <-  map_attr_to_comm(hardwood, comm_map, "decid")
# decid_rast[is.na(decid_rast[])] <- 0
# decid_smooth <- terra::focal(decid_rast, 5, fun = mean, na.rm = TRUE)
decid_coarse <- terra::aggregate(decid_rast, 5, na.rm = TRUE)
pred_decid_coarse <- terra::aggregate(predictor_stack2$prop_decid, 5, na.rm = TRUE)
pred_hardwood_coarse <- terra::aggregate(predictor_stack2$prop_hardwood, 5, na.rm = TRUE)
plot(decid_coarse)
plot(pred_hardwood_coarse)
plot(pred_hardwood_coarse[][samp_coarse] ~ decid_coarse[][samp_coarse])
abline(0,1)

hardwood_vals <- decid_coarse[][samp_coarse]
pred_vals <- predictor_stack_coarse$prop_hardwood[][samp_coarse]
hardwood_mod <- lm(pred_vals ~  hardwood_vals)
summary(hardwood_mod)
plot(effects::allEffects(hardwood_mod))

landis_crosswalk_mods <- bind_rows(landis_crosswalk_mods,
                                   tibble(var = "prop_hardwood", 
                                          mod = list(hardwood_mod)))

#---------------------------------------
# spruce
spruce <- comm_output
spruce$spruce = grepl("Pice|FrasFirr", comm_output$SpeciesName)
spruce_sum <- spruce %>%
  group_by(MapCode, spruce) %>%
  summarise(biomass = sum(CohortBiomass)) %>% 
  group_by(MapCode) %>%
  tidyr::pivot_wider(names_from = spruce, values_from = biomass, values_fill = 0) %>%
  mutate(spruce_prop = `TRUE`/(`TRUE`+ `FALSE`)) %>%
  dplyr::select(MapCode, spruce_prop)

spruce_rast <-  map_attr_to_comm(spruce_sum, comm_map, "spruce_prop")
spruce_coarse <- terra::aggregate(spruce_rast, 5, na.rm = TRUE)

spruce_vals <- spruce_coarse[][samp_coarse]
pred_vals <- predictor_stack_coarse$prop_spruce[][samp_coarse]
spruce_mod <- lm(pred_vals ~  spruce_vals)
summary(spruce_mod)
plot(effects::allEffects(spruce_mod))

landis_crosswalk_mods <- bind_rows(landis_crosswalk_mods,
                                   tibble(var = "prop_spruce", 
                                          mod = list(spruce_mod)))
#------------
#- oak
#-----------
oak <- comm_output
oak$oak = grepl("Quer|Cary", comm_output$SpeciesName)
oak_sum <- oak %>%
  group_by(MapCode, oak) %>%
  summarise(biomass = sum(CohortBiomass)) %>% 
  group_by(MapCode) %>%
  tidyr::pivot_wider(names_from = oak, values_from = biomass, values_fill = 0) %>%
  mutate(oak_prop = `TRUE`/(`TRUE`+ `FALSE`)) %>%
  dplyr::select(MapCode, oak_prop)

oak_rast <-  map_attr_to_comm(oak_sum, comm_map, "oak_prop")
oak_coarse <- terra::aggregate(oak_rast, 5, na.rm = TRUE)

oak_vals <- oak_coarse[][samp_coarse]
pred_vals <- predictor_stack_coarse$prop_oak[][samp_coarse]
plot(pred_vals ~ oak_vals)
oak_mod <- lm(pred_vals ~  oak_vals)
summary(oak_mod)
plot(effects::allEffects(oak_mod))

landis_crosswalk_mods <- bind_rows(landis_crosswalk_mods,
                                   tibble(var = "prop_oak", 
                                          mod = list(oak_mod)))

#------------------------------------------------
# Prop_forest
#-----------------------------------------------
forest <- lai > 0.5
prop_forest <- terra::focal(forest, 5, mean) %>% terra::aggregate(5, na.rm = TRUE)

plot(predictor_stack_coarse$prop_forest)
plot(prop_forest)

prop_forest_vals <- prop_forest[][samp_coarse]
pred_vals <- predictor_stack_coarse$prop_forest[][samp_coarse]
plot(pred_forest_coarse_vals ~ prop_forest_vals)
forest_mod <- lm(pred_vals ~ prop_forest_vals + 0)
summary(forest_mod)
plot(effects::allEffects(forest_mod))

landis_crosswalk_mods <- bind_rows(landis_crosswalk_mods,
                                   tibble(var = "prop_forest", 
                                          mod = list(forest_mod)))
        


#--------------------------------------
# compare climate variables
#--------------------------------------
#climate variables
year <- 18
if(year > 0){
  summer_clim <- climate_future %>%
    filter(Year > 2005 + year - 10 & Year <= 2005 + year, 
           Timestep > 182 & Timestep < 273) %>%
    group_by(EcoregionIndex) %>%
    summarise(tmax = median(max_airtemp)*10,#rescale to match GEE output
              tmin = median(min_airtemp)*10,
              precip = median(ppt)*1000) %>%
    mutate(EcoregionIndex = EcoregionIndex + 2)
  tmax_rast <- map_attr_to_ecoregions(summer_clim, ecoregions, "tmax")/10
  tmin_rast <- map_attr_to_ecoregions(summer_clim, ecoregions, "tmin")/10
  precip_rast <- map_attr_to_ecoregions(summer_clim, ecoregions, "precip")
}


plot(tmax_rast)
plot(predictor_stack2$tmmx)
plot(predictor_stack2$tmmx[][samp] ~ tmax_rast[][samp])
abline(0,1)
cor.test(predictor_stack2$tmmx[][samp], tmax_rast[][samp])

plot(tmin_rast)
plot(predictor_stack2$tmmn)
plot(predictor_stack2$tmmn[][samp] ~ tmin_rast[][samp])
abline(0,1)
cor.test(predictor_stack2$tmmn[][samp], tmin_rast[][samp])


plot(precip_rast)
plot(predictor_stack2$pr)
plot(predictor_stack2$pr[][samp] ~ precip_rast[][samp])
abline(0,1)
cor.test(predictor_stack2$pr[][samp], precip_rast[][samp])

ppt_vals <- precip_rast[][samp]
pred_vals <- predictor_stack$pr[][samp]
plot(pred_vals ~ ppt_vals)
ppt_mod <- lm(pred_vals ~ prop_forest_vals + 0)
summary(forest_mod)
plot(effects::allEffects(forest_mod))

ppt_mod <- lm(pred_vals ~ ppt_vals)
abline(lm(predictor_stack2$pr[][samp] ~ precip_rast[][samp]))
landis_crosswalk_mods <- bind_rows(landis_crosswalk_mods,
                                   tibble(var = "pr", 
                                          mod = list(ppt_mod)))


#-------------------------------
# pet, aet

necn_month <- necn_monthly %>%
                # filter(#Time > 2005 + year - 10 & Year <= 2005 + year, 
                #        Month %in% c(6,7,8)) %>%
                mutate(EcoregionIndex = ClimateRegionIndex + 2) %>%
                group_by(EcoregionIndex) %>%
                summarise(aet = median(avgAET)*10, #rescale to match GEE output
                          pet = median(PET)*10,
                          soil = median(SoilWaterContent)*10,
                          npp = median(AvgTotalNPP_C)) 

pet_rast <- map_attr_to_ecoregions(necn_month, ecoregions, "pet") %>% terra::aggregate(5, na.rm = TRUE)
aet_rast <- map_attr_to_ecoregions(necn_month, ecoregions, "aet") %>% terra::aggregate(5, na.rm = TRUE)
npp_rast <- map_attr_to_ecoregions(necn_month, ecoregions, "npp")
plot(npp_rast)

plot(predictor_stack_coarse$pet)
plot(pet_rast)
pet_vals <- pet_rast[][samp_coarse]
pred_vals <- predictor_stack_coarse$pet[][samp_coarse]
plot(pred_vals ~ pet_vals)
pet_mod <- lm(pred_vals ~ log(pet_vals))
summary(pet_mod)
plot(effects::allEffects(pet_mod))

landis_crosswalk_mods <- bind_rows(landis_crosswalk_mods,
                                   tibble(var = "pet", 
                                          mod = list(pet_mod)))

plot(predictor_stack_coarse$aet)
plot(aet_rast)
aet_vals <- aet_rast[][samp_coarse]
pred_vals <- predictor_stack_coarse$aet[][samp_coarse]
plot(pred_vals ~ aet_vals)
aet_mod <- lm(pred_vals ~ log(aet_vals))
summary(aet_mod)
plot(effects::allEffects(aet_mod))

landis_crosswalk_mods <- bind_rows(landis_crosswalk_mods,
                                   tibble(var = "aet", 
                                          mod = list(aet_mod)))


plot(npp_rast)
npp_rast_coarse <- terra::aggregate(npp_rast, 5, na.rm = TRUE)
plot(predictor_stack_coarse$npp)
plot(npp_rast_coarse)

npp_vals <- npp_rast_coarse[][samp_coarse]
pred_vals <- predictor_stack_coarse$npp[][samp_coarse]

plot(pred_vals ~ npp_vals)
npp_mod <- lm(pred_vals ~ npp_vals)
summary(npp_mod)
plot(effects::allEffects(npp_mod))

landis_crosswalk_mods <- bind_rows(landis_crosswalk_mods,
                                   tibble(var = "npp",
                                          mod = list(npp_mod)))

#----------------------
# soil moisture
#---------------------
soil_rast <- list.files("landis_analysis/reference/initial landscape layers", full.names = TRUE) %>% 
  `[`(grep("AvailableWater", .)) %>%
  terra::rast() %>%
  mean() * 10
plot(soil_rast)
soil_rast_coarse <- terra::aggregate(soil_rast, 5)
plot(predictor_stack_coarse$soil)
plot(soil_rast_coarse)

soil_vals <- soil_rast_coarse[][samp_coarse]
pred_vals <- predictor_stack_coarse$soil[][samp_coarse]

plot(pred_vals ~ soil_vals)
soil_mod <- lm(pred_vals ~ log(soil_vals))
summary(soil_mod)
plot(effects::allEffects(soil_mod))

landis_crosswalk_mods <- bind_rows(landis_crosswalk_mods,
                                   tibble(var = "soil", 
                                          mod = list(soil_mod)))



#----------------------
# npp (from rasters)
#---------------------
# npp_rast <- list.files("landis_analysis/reference/initial landscape layers", full.names = TRUE) %>%
#   `[`(grep("NPP", .)) %>%
#   terra::rast() %>%
#   mean() * 10
# plot(npp_rast)
# npp_rast_coarse <- terra::aggregate(npp_rast, 5, na.rm = TRUE)
# plot(predictor_stack_coarse$npp)
# plot(npp_rast_coarse)
# 
# npp_vals <- npp_rast_coarse[][samp_coarse]
# pred_vals <- predictor_stack_coarse$npp[][samp_coarse]
# 
# plot(pred_vals ~ npp_vals)
# npp_mod <- lm(pred_vals ~ npp_vals)
# summary(npp_mod)
# plot(effects::allEffects(npp_mod))
# 
# landis_crosswalk_mods <- bind_rows(landis_crosswalk_mods,
#                                    tibble(var = "prop_forest",
#                                           mod = list(forest_mod)))


saveRDS(landis_crosswalk_mods, "./landis_analysis/landis_predictor_layers/landis_crosswalk_mods.RDS")
