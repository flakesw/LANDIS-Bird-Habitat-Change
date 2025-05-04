# Make new predictor layers from LANDIS-II outputs
# Used to make predictions in make_maps_from_landis_outputs.R

#TODO for v2, revise initial climate

library("terra")
library("sf")
library("tidyverse")


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


#"HighTHighV BAU extra Rx"
model_parent_dir <- "C:/Users/swflake/Documents/SApps-LANDIS/Model runs/"
# model_names <- c("LowTLowV BAU", "HighTHighV BAU")
model_names <- list.dirs(model_parent_dir, recursive = FALSE, full.names = FALSE)
model_names <- model_names[-c(6,7)]

input_dir <- "C:/Users/swflake/Documents/SApps LANDIS/Inputs/" #"D:/SApps LANDIS/Inputs/"
years <- c(0, 10, 20, 30, 40, 50, 60, 70, 80)

#----------------------------------------------
#species reference data
#-----------------------------------------------
sp_ref <- read.csv("./landis_analysis/reference/REF_SPECIES.csv")

sp_data <- read.csv("./landis_analysis/reference/NECN_Spp_Table-NECN_v6.csv")
spp_to_use <- sp_data$SpeciesCode
spp_to_use_all <- sp_data$SpeciesCode
sp_ref$SpeciesCode <- paste0(substr(sp_ref$GENUS, 1, 4), substr(sp_ref$SPECIES, 1, 4) %>%
                               stringr::str_to_title())

spp_crosswalk <- sp_ref[sp_ref$SpeciesCode %in% spp_to_use_all, ] %>%
  dplyr::arrange(SpeciesCode) %>%
  dplyr::select(SpeciesCode, SPCD) %>%
  dplyr::filter(!(SPCD %in% c(318, 114, 8349, 845, 953, 952)))
spp_crosswalk <- rbind(spp_crosswalk, #fix some names which don't match SPCDs
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

crosswalk_mods <- readRDS("./landis_analysis/landis_predictor_layers/landis_crosswalk_mods.RDS")




#-----------------------------------
# the big function that does everything
#---------------------------------

make_predictor_stack_from_landis <- function(model_name, input_dir, year)
{

  model_dir <- paste0(model_parent_dir, model_name, "/")
  
  comm_output <- read.csv(paste0(model_dir, "/community-input-file-", year,".csv"))
  comm_output[comm_output[] == 0] <- NA
  comm_map <- ecoregions
  comm_map[] <- terra::rast(paste0(model_dir, "/output-community-", year, ".img"))[]
  plot(comm_map)
  
  climate_future <- read.csv(paste0(model_dir, "Climate-future-input-log.csv"))
  ecoregions <- terra::rast("./landis_analysis/reference/Ecos11_NCLD.tif")
  
  necn_monthly <- read.csv(paste0(model_dir, "/NECN-succession-monthly-log.csv"))
  
  predictor_stack <- terra::rast("./landis_analysis/landis_predictor_layers/predictor_stack_bcr28_train.tif")
  predictor_stack2 <- predictor_stack %>%
    project(ecoregions, mask = TRUE, method = "near", align = FALSE) %>%
    mask(ecoregions, maskvalues = 1)

  #------------------
  # height
  # used for several of the variables
  #------------------
  
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
  
  #---------------------------------------
  # rh98
  #---------------------------------------
  height_rast_98 <- comm_height %>%
    group_by(MapCode) %>%
    summarise(height = quantile(HT, 0.98))
  height_rast_98 <- map_attr_to_comm(height_rast_98, comm_map, "height")
  plot(height_rast_98)
  plot(predictor_stack2$height)
  mod <-  crosswalk_mods[crosswalk_mods$var == "height", "mod"][[1]][[1]]
  new_vals <- predict(mod, newdata = list(height_vals = as.numeric(height_rast_98[])))
  new_vals <- ifelse(new_vals < 0, 0, new_vals)
  height2 <- height_rast_98
  height2[] <- new_vals[]
  
  #--------------------------------------------
  #biomass
  #--------------------------------------------
  total_biomass <- comm_output %>%
    group_by(MapCode) %>%
    summarise(biomass = sum(CohortBiomass))
  biomass_rast <- map_attr_to_comm(total_biomass, comm_map, "biomass") / 100 #convert to tonnes ha-1
  plot(biomass_rast)
  plot(predictor_stack2$biomass)
  mod <-  crosswalk_mods[crosswalk_mods$var == "biomass", "mod"][[1]][[1]]
  new_vals <- predict(mod, newdata = list(biomass_vals = as.numeric(biomass_rast[])))
  new_vals <- ifelse(new_vals < 0, 0, new_vals)
  biomass2 <- biomass_rast
  biomass2[] <- new_vals[]
  
  #----------------------------------------
  # understory_ratio
  #-----------------------------------------
  under_ratio <- comm_height %>%
    mutate(bin = cut(HT, breaks = seq(0, 100, 5), labels = FALSE)) %>%
    group_by(MapCode, bin) %>%
    mutate(bio_height = sum(CohortBiomass)) %>%
    group_by(MapCode, bin) %>%
    summarise(biomass_sum_bin = sum(bio_height)) %>%
    ungroup() %>%
    tidyr::complete(MapCode, bin, fill = list(biomass_sum_bin = 0)) %>%
    group_by(MapCode) %>%
    summarise(understory_4 = sum(biomass_sum_bin[bin <= 4])/sum(biomass_sum_bin))
  under_rast_4 <- map_attr_to_comm(under_ratio, comm_map, "understory_4")
  
  year_lai <- ifelse(year < 5, 5, year)
  lai <- terra::rast(paste0(model_dir, "NECN/LAI-", year_lai, ".img"))
  
  mod <-  crosswalk_mods[crosswalk_mods$var == "understory", "mod"][[1]][[1]]
  new_vals <- predict(mod, newdata = list(under_vals = as.numeric(under_rast_4[]),
                                          lai_vals = as.numeric(lai[])))
  under2 <- under_rast_4
  under2[] <- new_vals[]
  
  #-------------------------------------
  #foliage height diversity
  #------------------------------------
  fhd <- comm_height %>%
    mutate(bin = cut(HT, breaks = seq(0, 100, 2), labels = FALSE)) %>%
    group_by(MapCode, bin) %>%
    summarize(bio_height = sum(log(CohortBiomass))) %>%
    group_by(MapCode) %>%
    mutate(shan = -1*(bio_height/sum(bio_height)) * log(bio_height/sum(bio_height))) %>%
    summarize(ht_diversity = sum(shan))
  fhd_rast <-  map_attr_to_comm(fhd, comm_map, "ht_diversity")
  mod <-  crosswalk_mods[crosswalk_mods$var == "fhd", "mod"][[1]][[1]]
  new_vals <- predict(mod, newdata = list(fhd_vals = as.numeric(fhd_rast[])))
  fhd2 <- fhd_rast
  fhd2[] <- new_vals[]
  
  
  #---------------------------------------------
  # area_short
  #---------------------------------------------
  area_short <- comm_height %>%
    group_by(MapCode) %>%
    summarize(short = as.numeric(quantile(HT, 1) < 20))
  
  area_short_rast <- map_attr_to_comm(area_short, comm_map, "short")
  area_short_rast <- terra::focal(area_short_rast, 7, fun = "mean", 
                                   na.policy = "omit", na.rm = TRUE, 
                                   expand = TRUE, fillvalue = 0)
  # lai_radius <- terra::focal(lai, 5, fun = "mean")
  # 
  # mod <-  crosswalk_mods[crosswalk_mods$var == "area_short", "mod"][[1]][[1]]
  # new_vals <- predict(mod, newdata = list(area_short_vals = as.numeric(area_short_rast[]),
  #                                         lai_vals = as.numeric(lai_radius[])))
  area_short2 <- area_short_rast
  # area_short2[] <- new_vals
  
  
  #----------------------------------------
  # Prop hardwood
  #----------------------------------------
  #proportion hardwood -- trianing data from Treemap
  decid_all <- comm_output
  decid_all$decid = !grepl("Pinu|Tsug|Pice|FrasFirr", comm_output$SpeciesName)
  decid_sum <- decid_all %>%
    group_by(MapCode, decid) %>%
    summarise(biomass = sum(CohortBiomass))
  hardwood <- decid_sum %>%
    group_by(MapCode) %>%
    tidyr::pivot_wider(names_from = decid, values_from = biomass, values_fill = 0) %>%
    mutate(decid = `TRUE`/(`TRUE`+ `FALSE`)) %>%
    dplyr::select(MapCode, decid)
  hardwood_rast <-  map_attr_to_comm(hardwood, comm_map, "decid")
  mod <-  crosswalk_mods[crosswalk_mods$var == "prop_hardwood", "mod"][[1]][[1]]
  new_vals <- predict(mod, newdata = list(hardwood_vals = as.numeric(hardwood_rast[])))
  hardwood2 <- hardwood_rast
  hardwood2[] <- new_vals[]
  
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
  mod <-  crosswalk_mods[crosswalk_mods$var == "prop_spruce", "mod"][[1]][[1]]
  new_vals <- predict(mod, newdata = list(spruce_vals = as.numeric(spruce_rast[])))
  spruce2 <- spruce_rast
  spruce2[] <- new_vals[]
  
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
  mod <-  crosswalk_mods[crosswalk_mods$var == "prop_oak", "mod"][[1]][[1]]
  new_vals <- predict(mod, newdata = list(oak_vals = as.numeric(oak_rast[])))
  oak2 <- oak_rast
  oak2[] <- new_vals[]
  
  #------------------------------------------------
  # Prop_forest
  #-----------------------------------------------
  #pretty good -- no crosswalk needed
  forest <- lai > 0.5
  prop_forest <- comm_map
  prop_forest[] <-  terra::focal(forest, 5, mean, na.rm = TRUE)[]
  prop_forest <- mask(prop_forest, comm_map, maskvalue = 0)
  
  # mod <-  crosswalk_mods[crosswalk_mods$var == "prop_forest", "mod"][[1]][[1]]
  # new_vals <- predict(mod, newdata = list(prop_forest_vals = as.numeric(prop_forest[])))
  forest2 <- ecoregions
  forest2[] <- prop_forest[]
  # forest2[] <- new_vals[]
  
  #--------------------------------------
  #climate variables
  
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
  
  
  #--------------------------------------
  #monthly NECN data
  # pet, aet, npp
  #---------------------------------------
  
  if(year > 0){
    necn_month <- necn_monthly %>%
      filter(Time > year - 10 & Time <= year,
             Month %in% c(6,7,8)) %>%
      mutate(EcoregionIndex = ClimateRegionIndex + 2) %>%
      group_by(EcoregionIndex) %>%
      summarise(aet = median(avgAET)*10, #rescale to match GEE output
                pet = median(PET)*10,
                soil = median(SoilWaterContent)*10,
                npp = median(AvgTotalNPP_C)) 
    
    pet_rast <- map_attr_to_ecoregions(necn_month, ecoregions, "pet")
    mod <-  crosswalk_mods[crosswalk_mods$var == "pet", "mod"][[1]][[1]]
    new_vals <- predict(mod, newdata = list(pet_vals = as.numeric(pet_rast[])))
    pet2 <- pet_rast
    pet2[] <- new_vals[]
    
    aet_rast <- map_attr_to_ecoregions(necn_month, ecoregions, "aet")
    mod <-  crosswalk_mods[crosswalk_mods$var == "aet", "mod"][[1]][[1]]
    new_vals <- predict(mod, newdata = list(aet_vals = as.numeric(aet_rast[])))
    aet2 <- aet_rast
    aet2[] <- new_vals[]
    
    npp_rast <- map_attr_to_ecoregions(necn_month, ecoregions, "npp")
    mod <-  crosswalk_mods[crosswalk_mods$var == "npp", "mod"][[1]][[1]]
    new_vals <- predict(mod, newdata = list(npp_vals = as.numeric(npp_rast[])))
    npp2 <- npp_rast
    npp2[] <- new_vals[]
  }
  
  # if(year == 0){
  #   necn_month <- necn_monthly %>%
  #     filter(Time > year - 10 & Time <= year,
  #            Month %in% c(6,7,8)) %>%
  #     mutate(EcoregionIndex = ClimateRegionIndex + 2) %>%
  #     group_by(EcoregionIndex) %>%
  #     summarise(aet = median(avgAET)*10, #rescale to match GEE output
  #               pet = median(PET)*10,
  #               soil = median(SoilWaterContent)*10,
  #               npp = median(AvgTotalNPP_C)) 
  #   
  #   pet_rast <- map_attr_to_ecoregions(necn_month, ecoregions, "pet")
  #   mod <-  crosswalk_mods[crosswalk_mods$var == "pet", "mod"][[1]][[1]]
  #   new_vals <- predict(mod, newdata = list(pet_vals = as.numeric(pet_rast[])))
  #   pet2 <- pet_rast
  #   pet2[] <- new_vals[]
  #   
  #   aet_rast <- map_attr_to_ecoregions(necn_month, ecoregions, "aet")
  #   mod <-  crosswalk_mods[crosswalk_mods$var == "aet", "mod"][[1]][[1]]
  #   new_vals <- predict(mod, newdata = list(aet_vals = as.numeric(aet_rast[])))
  #   aet2 <- aet_rast
  #   aet2[] <- new_vals[]
  #   
  #   npp_rast <- map_attr_to_ecoregions(necn_month, ecoregions, "npp")
  #   mod <-  crosswalk_mods[crosswalk_mods$var == "npp", "mod"][[1]][[1]]
  #   new_vals <- predict(mod, newdata = list(npp_vals = as.numeric(npp_rast[])))
  #   npp2 <- npp_rast
  #   npp2[] <- new_vals[]
  # }
  
  #----------------------
  # soil moisture
  #---------------------
  soil_rast <- terra::rast(paste0(model_dir, "NECN/AvailableWater-", year_lai, ".img")) * 10
  plot(soil_rast)
  mod <-  crosswalk_mods[crosswalk_mods$var == "soil", "mod"][[1]][[1]]
  new_vals <- predict(mod, newdata = list(soil_vals = as.numeric(soil_rast[])))
  soil2 <- ecoregions
  soil2[] <- new_vals[]
  
  
  

  #stack all the predictors
  if(year == 0){ #year 0 -- use initial communities from LANDIS but training data for climate
    pred_stack <- c(biomass2, 
                    under2, 
                    area_short2,
                    fhd2,
                    height2,
                    forest2,
                    spruce2,
                    oak2,
                    hardwood2,
                    predictor_stack2$aet,
                    predictor_stack2$pet,
                    predictor_stack2$soil,
                    predictor_stack2$tmmn,
                    predictor_stack2$tmmx,
                    predictor_stack2$pr,
                    predictor_stack2$slope,
                    predictor_stack2$chili,
                    predictor_stack2$tpi,
                    predictor_stack2$npp)
  } else{
    pred_stack <- c(biomass2, 
                    under2, 
                    area_short2,
                    fhd2,
                    height2,
                    forest2,
                    spruce2,
                    oak2,
                    hardwood2,
                    aet2,
                    pet2,
                    soil2,
                    tmin_rast,
                    tmax_rast,
                    precip_rast,
                    predictor_stack2$slope,
                    predictor_stack2$chili,
                    predictor_stack2$tpi,
                    predictor_stack2$npp)
  }
  
  names(pred_stack) <- c("biomass", "understory", "area_short", "fhd", "height",
                         "prop_forest", "prop_spruce", "prop_oak", "prop_hardwood", 
                         "aet", "pet", "soil", 
                         "tmmn", "tmmx", "pr", "slope", "chili", "tpi", "Npp")
  
  
  writeRaster(pred_stack, 
              filename = paste0("./landis_analysis/landis_predictor_layers/pred_stack_", 
                                            model_name, "_", year, ".tif"),
              overwrite = TRUE)
}


#--------------------------------------------------
# run the script
#-------------------------------------------------
for(model in model_names){
  for(year in years){
    make_predictor_stack_from_landis(model_name = model, input_dir = input_dir, year)
  }
}
