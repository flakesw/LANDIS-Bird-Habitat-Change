# compare prediction maps between timesteps
library("terra")
library("tidyverse")
library("tidyterra")
library("colorspace")
library("cowplot")
library("lemon")
library("gridExtra")
library("sf")

theme_set(theme_minimal())
theme_set(theme(legend.position = "bottom"))

states <- sf::st_read("./maps/cb_2018_us_state_5m/cb_2018_us_state_5m.shp")
study_area <- sf::st_read("./maps/study_area.gpkg")


pred_dir <- "./landis_analysis/landis_predictions/"
all_pred_layers <- list.files(pred_dir) %>% `[`(grep("tif", .))
all_pred_locs <- list.files(pred_dir, full.names = TRUE) %>% `[`(grep("tif", .))
pred_info <- tibble(pred_loc = all_pred_locs,
                    pred_run = str_split(all_pred_layers, pattern = "[_.]+") %>% map(pluck(4)) %>% unlist(),
                    species = str_split(all_pred_layers, pattern = "_") %>% map(pluck(3)) %>% unlist()) %>%
  mutate(climate = str_split(all_pred_layers, pattern = "[_. ]+") %>% map(pluck(4)) %>% unlist(),
         replicate = str_split(all_pred_layers, pattern = "[_. ]+") %>% map(pluck(7)) %>% unlist(),
         year = str_split(all_pred_layers, pattern = "[_. ]+") %>% map(pluck(8)) %>% unlist(),
         type = str_split(all_pred_layers, pattern = "[_. ]+") %>% map(pluck(1)) %>% unlist()) %>%
  filter(year %in% c(10,80)) %>%
  filter(replicate == "Run1")

pred_results <- tibble(pred_loc = all_pred_locs,
                    pred_run = str_split(all_pred_layers, pattern = "[_.]+") %>% map(pluck(4)) %>% unlist(),
                    species = str_split(all_pred_layers, pattern = "_") %>% map(pluck(3)) %>% unlist()) %>%
  mutate(climate = str_split(all_pred_layers, pattern = "[_. ]+") %>% map(pluck(4)) %>% unlist(),
         replicate = str_split(all_pred_layers, pattern = "[_. ]+") %>% map(pluck(7)) %>% unlist(),
         year = str_split(all_pred_layers, pattern = "[_. ]+") %>% map(pluck(8)) %>% unlist(),
         type = str_split(all_pred_layers, pattern = "[_. ]+") %>% map(pluck(1)) %>% unlist()) %>%
  filter(year %in% c(80)) %>%
  filter(replicate == "Run1") %>%
  mutate(change_in_hi = vector(mode = "list", length = n()))

species_list <- unique(pred_info$species)

for(sp in species_list){
  
  pred_info_sp <- pred_info %>% filter(species == sp)
  
  for(clim in unique(pred_info_sp$climate)){
    pred_info_init <- pred_info_sp[pred_info_sp$type == "full" & pred_info_sp$climate == clim & pred_info_sp$year == 10, ]
    pred_info_end <- pred_info_sp[pred_info_sp$type == "full" & pred_info_sp$climate == clim & pred_info_sp$year == 80, ]
                                
    for(run in unique(pred_info_init$replicate)){
      preds_10 <- rast(pred_info_init[pred_info_init$replicate == run, "pred_loc"][[1]]) 
      # plot(preds_10)
      preds_80 <- rast(pred_info_end[pred_info_init$replicate == run, "pred_loc"][[1]]) 
      # plot(preds_80)
      change_in_hi <- preds_80 - preds_10
      # plot(change_in_hi, main = sp)
      
      pred_results[pred_results$replicate == run & pred_results$climate == clim & pred_results$species== sp, ] $change_in_hi <- list(change_in_hi)
      
    }
  }
}


  
  
mfrow(c(5,4))
pred_select <- pred_results %>% filter(type == "full",
                                       climate == "HighTHighV") 
change_in_hi_stack <- rast(pred_select$change_in_hi)
names(change_in_hi_stack) <- pred_select$species %>% toupper()

states <- sf::st_transform(states, crs = crs(change_in_hi_stack))
states_study_area <- st_bbox(buffer(change_in_hi_stack, 100)) %>% 
  st_as_sfc()%>%
  st_intersection(states, .)

change_in_hi_high_full <- ggplot() + 
  # geom_spatraster(data = change_in_hi_stack, show.legend = TRUE) +
  geom_sf(data = study_area, fill = "white", linewidth = 0) +
  geom_spatraster(data = change_in_hi_stack, show.legend = TRUE) +
  geom_sf(data = study_area, fill = NA, linewidth = 1) +
  # geom_sf(data = states_study_area, fill = NA, alpha = 0.5)+
  facet_wrap(~lyr)+
  scale_fill_continuous_divergingx(palette = 'RdBu',
                                   rev = FALSE, mid = 0,
                                   na.value = NA,
                                   guide = "colourbar",
                                   limits = c(-1,1)) +
  labs(fill = "Change in Habitat Index") +
  scale_x_continuous(breaks=scales::pretty_breaks(n=3))

plot(change_in_hi_high_full)


#-----------------------
# contribution of vegetation
# TODO
for(sp in species_list){
  
  pred_info_sp <- pred_info %>% filter(species == sp)
  
  for(clim in unique(pred_info_sp$climate)){
    pred_info_init <- pred_info_sp[pred_info_sp$type == "full" & pred_info_sp$climate == clim & pred_info_sp$year == 10, ]
    pred_info_end <- pred_info_sp[pred_info_sp$type == "full" & pred_info_sp$climate == clim & pred_info_sp$year == 80, ]
    
    for(run in unique(pred_info_init$replicate)){
      preds_10 <- rast(pred_info_init[pred_info_init$replicate == run, "pred_loc"][[1]]) 
      # plot(preds_10)
      preds_80 <- rast(pred_info_end[pred_info_init$replicate == run, "pred_loc"][[1]]) 
      # plot(preds_80)
      change_in_hi <- preds_80 - preds_10
      # plot(change_in_hi, main = sp)
      
      pred_results[pred_results$replicate == run & pred_results$climate == clim & pred_results$species== sp, ] $change_in_hi <- list(change_in_hi)
      
    }
  }
}

change_in_hi_stack <- rast(pred_select$change_in_hi)
names(change_in_hi_stack) <- pred_select$species %>% toupper()

states <- sf::st_transform(states, crs = crs(change_in_hi_stack))
states_study_area <- st_bbox(buffer(change_in_hi_stack, 100)) %>% 
  st_as_sfc()%>%
  st_intersection(states, .)

change_in_hi_high_full <- ggplot() + 
  # geom_spatraster(data = change_in_hi_stack, show.legend = TRUE) +
  geom_sf(data = study_area, fill = "white", linewidth = 0) +
  geom_spatraster(data = change_in_hi_stack, show.legend = TRUE) +
  geom_sf(data = study_area, fill = NA, linewidth = 1) +
  # geom_sf(data = states_study_area, fill = NA, alpha = 0.5)+
  facet_wrap(~lyr)+
  scale_fill_continuous_divergingx(palette = 'RdBu',
                                   rev = FALSE, mid = 0,
                                   na.value = NA,
                                   guide = "colourbar",
                                   limits = c(-1,1)) +
  labs(fill = "Change in Habitat Index")

plot(change_in_hi_high_full)


#----------------------------------
# find refugia locations

pred_results_refugia <- tibble(pred_loc = all_pred_locs,
                               pred_run = str_split(all_pred_layers, pattern = "[_.]+") %>% map(pluck(4)) %>% unlist(),
                               species = str_split(all_pred_layers, pattern = "_") %>% map(pluck(3)) %>% unlist()) %>%
                      mutate(climate = str_split(all_pred_layers, pattern = "[_. ]+") %>% map(pluck(4)) %>% unlist(),
                             replicate = str_split(all_pred_layers, pattern = "[_. ]+") %>% map(pluck(7)) %>% unlist(),
                             year = str_split(all_pred_layers, pattern = "[_. ]+") %>% map(pluck(8)) %>% unlist(),
                             type = str_split(all_pred_layers, pattern = "[_. ]+") %>% map(pluck(1)) %>% unlist()) %>%
                      filter(year %in% c(80)) %>%
                      filter(replicate == "Run1") %>%
                      mutate(in_situ_ref = vector(mode = "list", length = n()),
                             ex_situ_ref = vector(mode = "list", length = n()),
                             lost_habitat = vector(mode = "list", length = n()),
                             combined_ref = vector(mode = "list", length = n()))

for(sp in species_list){
    
    pred_info_sp <- pred_info %>% filter(species == sp)
    
    for(clim in unique(pred_info_sp$climate)){
      pred_info_init <- pred_info_sp[pred_info_sp$type == "full" & pred_info_sp$climate == clim & pred_info_sp$year == 10, ]
      pred_info_end <- pred_info_sp[pred_info_sp$type == "full" & pred_info_sp$climate == clim & pred_info_sp$year == 80, ]
      pred_info_init_clim <- pred_info_sp[pred_info_sp$type == "clim" & pred_info_sp$climate == clim & pred_info_sp$year == 10, ]
      pred_info_end_clim <- pred_info_sp[pred_info_sp$type == "clim" & pred_info_sp$climate == clim & pred_info_sp$year == 80, ]
      pred_info_init_veg <- pred_info_sp[pred_info_sp$type == "veg" & pred_info_sp$climate == clim & pred_info_sp$year == 10, ]
      pred_info_end_veg <- pred_info_sp[pred_info_sp$type == "veg" & pred_info_sp$climate == clim & pred_info_sp$year == 80, ]
      
      for(run in unique(pred_info_init$replicate)){
        preds_10 <- rast(pred_info_init[pred_info_init$replicate == run, "pred_loc"][[1]]) 
        preds_10[] <- boot::logit(preds_10[])
        preds_10_clim <- rast(pred_info_init_clim[pred_info_init_clim$replicate == run, "pred_loc"][[1]]) 
        preds_10_clim[] <- boot::logit(preds_10_clim[])
        preds_10_veg <- rast(pred_info_init_veg[pred_info_init_veg$replicate == run, "pred_loc"][[1]]) 
        preds_10_veg[] <- boot::logit(preds_10_veg[])
        # plot(preds_10)
        preds_80 <- rast(pred_info_end[pred_info_end$replicate == run, "pred_loc"][[1]]) 
        preds_80[] <- boot::logit(preds_80[])
        preds_80_clim <- rast(pred_info_end_clim[pred_info_end_clim$replicate == run, "pred_loc"][[1]]) 
        preds_80_clim[] <- boot::logit(preds_80_clim[])
        preds_80_veg <- rast(pred_info_end_veg[pred_info_end_veg$replicate == run, "pred_loc"][[1]]) 
        preds_80_veg[] <- boot::logit(preds_80_veg[])
        
        in_situ_ref <- preds_10 > 0 & preds_80 > 0
        in_situ_ref[which(in_situ_ref[] == TRUE)] <- 1
        ex_situ_ref <- preds_10 < 0 & preds_80 > 0
        ex_situ_ref[which(ex_situ_ref[] == TRUE)] <- 2
        lost_habitat <- preds_10 > 0 & preds_80 < 0
        lost_habitat[which(lost_habitat[] == TRUE)] <- 3
        
        combined_ref <- in_situ_ref + ex_situ_ref + lost_habitat
        combined_ref[combined_ref[] == 0] <- NA
        combined_ref <- as.factor(combined_ref)
        cls <- data.frame(id=1:3, cover=c("In-situ refugium", "Ex-situ refugium", "Habitat loss"))
        levels(combined_ref) <- cls
        
        # 
        # change_due_to_veg <- (preds_80_clim - preds_10_clim)
        # change_due_to_clim <- (preds_80_veg - preds_10_veg)
        # 
        # end_veg_suitability <- preds_10 + change_due_to_veg
        # end_veg_suitability[] <- boot::inv.logit(end_veg_suitability[])
        # end_clim_suitability <- preds_10 + change_due_to_clim
        # end_clim_suitability[] <- boot::inv.logit(end_clim_suitability[])
        # potential_ref <- preds_80_veg > 0.5 & preds_80_clim < 0.5
        
        pred_results_refugia[pred_results_refugia$replicate == run & pred_results_refugia$climate == clim & pred_results_refugia$species== sp, ]$in_situ_ref <- list(in_situ_ref)
        pred_results_refugia[pred_results_refugia$replicate == run & pred_results_refugia$climate == clim & pred_results_refugia$species== sp, ]$ex_situ_ref <- list(ex_situ_ref)
        pred_results_refugia[pred_results_refugia$replicate == run & pred_results_refugia$climate == clim & pred_results_refugia$species== sp, ]$lost_habitat <- list(lost_habitat)
        pred_results_refugia[pred_results_refugia$replicate == run & pred_results_refugia$climate == clim & pred_results_refugia$species== sp, ]$combined_ref <- list(combined_ref)
        
        }
    }
  }
  

pred_select_ref <- pred_results_refugia %>% filter(type == "full",
                                       climate == "HighTHighV") 
refugia_stack <- rast(pred_select_ref$combined_ref)
names(refugia_stack) <- pred_select_ref$species %>% toupper()

# states <- sf::st_transform(states, crs = crs(change_in_hi_stack))
# states_study_area <- st_bbox(buffer(change_in_hi_stack, 100)) %>% 
#   st_as_sfc()%>%
#   st_intersection(states, .)

refugia<- ggplot() + 
  # geom_spatraster(data = change_in_hi_stack, show.legend = TRUE) +
  geom_sf(data = study_area, fill = "lightgrey", linewidth = 0) +
  geom_spatraster(data = refugia_stack, show.legend = TRUE) +
  geom_sf(data = study_area, fill = NA, linewidth = 1) +
  # geom_sf(data = states_study_area, fill = NA, alpha = 0.5)+
  facet_wrap(~lyr) +
  scale_fill_viridis_d(na.value = "NA", na.translate = F)+
  labs(fill = "") 

  # scale_fill_continuous_divergingx(palette = 'RdBu',
  #                                  rev = FALSE, mid = 0,
  #                                  na.value = NA,
  #                                  guide = "colourbar",
  #                                  limits = c(-1,1)) +
  # labs(fill = "Change in Habitat Index")

plot(refugia)

  
in_situ_ref_cerw <- preds_0_cerw > 0.5 & preds_60_cerw > 0.5
in_situ_ref_cerw <- as.numeric(in_situ_ref_cerw)
plot(in_situ_ref_cerw)
ex_situ_ref_cerw <- as.numeric(preds_0_cerw < 0.5 & preds_60_cerw > 0.5)
plot(ex_situ_ref_cerw)
ex_situ_ref_cerw[ex_situ_ref_cerw[] == TRUE] <- 2
plot(in_situ_ref_cerw + ex_situ_ref_cerw)

# #potential refugia
# pot_ref_cerw <- as.numeric(preds_60_cerw < 0.5 & preds_60_veg_cerw > 0.5)
# pot_ref_cerw[pot_ref_cerw[] == TRUE] <- 3
# plot(pot_ref_cerw)
# 
# pot_ref_gwwa <- as.numeric(preds_60_gwwa < 0.5 & preds_60_veg_gwwa > 0.5)
# pot_ref_gwwa[pot_ref_gwwa[] == TRUE] <- 3
# plot(pot_ref_gwwa)

# lost habitat
lost_cerw <- as.numeric(preds_0_cerw > 0.5 & preds_60_cerw < 0.5)
lost_cerw[lost_cerw[] == TRUE] <- 3
plot(lost_cerw)
plot(lost_cerw + in_situ_ref_cerw + ex_situ_ref_cerw)


#-----------------------------------
# draft figure
coltb <- data.frame(value=0:3, col=c("lightgray", "#1E6D60", "#1E88E5", "#FF609A"))


ref_cerw <- in_situ_ref_cerw + ex_situ_ref_cerw + lost_cerw
coltab(ref_cerw) <- coltb
plot(ref_cerw)


#----------------------------------
#characteristics of refugia?

pred_stack_0 <- terra::rast("./landis_analysis/landis_predictor_layers/pred_stack_study_area.tif")

in_situ_refugia_predictors <- terra::mask(pred_stack_0, ref_cerw, maskvalue = 1, inverse = TRUE)
ex_situ_refugia_predictors <- terra::mask(pred_stack_0, ref_cerw, maskvalue = 2, inverse = TRUE)
lost_predictors <- terra::mask(pred_stack_0, ref_cerw, maskvalue = 3, inverse = TRUE)

ref_df <- data.frame(in_situ_refugia_predictors) %>%
  slice_sample(n = 10000) %>%
  mutate(refugia = "in_situ")
n_ref_df <- data.frame(ex_situ_refugia_predictors) %>%
  slice_sample(n = 10000)%>%
  mutate(refugia = "ex_situ")
lost_df <- data.frame(lost_predictors) %>%
  slice_sample(n = 10000)%>%
  mutate(refugia = "lost")
ref_preds_data <- rbind(ref_df, n_ref_df, lost_df)

test <- glm(refugia ~ ., data = ref_preds_data)
summary(test)

for(i in 1:22){
  formula <- paste0(names(ref_preds_data)[i], "~ refugia")

  vioplot::vioplot(as.formula(formula), 
                   data = ref_preds_data, 
                   main =names(ref_preds_data)[i], 
                 col = c("#1E88E5", "#1E6D60", "#FF609A"))
  
}
