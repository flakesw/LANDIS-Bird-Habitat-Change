#-----------------------------
#make predictions surfaces for bird species using predictors from LANDIS
#requires predictor layers created in create_predictors_from_landis_all_mods.R

#full model
library("terra")
library("sf")
library("tidyverse")
library("dismo")
library("gbm")
library("tidyterra")

diverging_color_ramp <- function(ras){
  the_palette_fc <- leaflet::colorNumeric(palette = "RdBu", 
                                          # domain = c(-max(abs(ras[]), na.rm = TRUE), max(abs(ras[]), na.rm = TRUE)),
                                          domain = c(-1, 1),
                                          reverse = FALSE)
  # the_colors <- the_palette_fc(seq((-max(abs(ras[]), na.rm = TRUE)), max(abs(ras[]), na.rm = TRUE), length.out = 50))
  the_colors <- the_palette_fc(seq(-1, 1, length.out = 50))

}

theme_set(theme_gray())
# theme_update(panel.border = element_blank(),
#              panel.grid.major = element_blank(),
#              panel.grid.minor = element_blank(),
#              legend.position = "none")
theme_update(legend.position = "none")

ecoregions <- terra::rast("./landis_analysis/reference/Ecos11_NCLD.tif")


pred_dir <- "./landis_analysis/landis_predictor_layers/landis_predictors/"
all_pred_stacks <- list.files(pred_dir)
pred_info <- tibble(model_loc = all_pred_stacks,
                    model = str_split(all_pred_stacks, pattern = "_") %>% map(pluck(3)) %>% unlist()) %>%
  mutate(climate = str_split(model, " ") %>% map(pluck(2)) %>% unlist(),
         replicate = str_split(model, " ") %>% map(pluck(4)) %>% unlist(),
         year = str_split(all_pred_stacks, pattern = "[_.]+") %>% map(pluck(4)) %>% unlist())

predictor_stack_orig <- rast("./landis_analysis/landis_predictor_layers/predictor_stack_bcr28_train.tif") %>%
  project(ecoregions, mask = TRUE, method = "near", align = FALSE) %>%
  mask(ecoregions, maskvalues = 1)
names(predictor_stack_orig)[30] <- "Npp"


species_list <- c("bcch", "blbw",  "cerw", "gwwa", "heth", "kewa", "lowa",
                  "praw", "prow", "recr", "swwa", "veer", "wewa","ytwa",
                  "ewpw", "rugr", "nswo",  "acfl", "alfl", "bhnu")


for(i in 1:length(species_list)){
  species <- species_list[i]
  mod <- list.files("./landis_analysis/models_for_landis/", full.names = TRUE) %>%
    `[`(grep(species, .)) %>%
    readRDS()
  mod <- mod$models[[1]]

  for(j in 1:nrow(pred_info)){#----------------------------------------------------------------
  #updated predictions for given year/model
    
    predictor_stack <- paste0(pred_dir, pred_info$model_loc[j]) %>%
      terra::rast() %>%
      terra::project(ecoregions) %>%
      mask(ecoregions, maskvalues = 1)
    
    
    preds <- terra::predict(object = predictor_stack, 
                            model = mod,
                            # ext = st_bbox(bcr_albers),
                            const = data.frame(time_observations_started = 7,
                                               duration_minutes = 60)
                            )
    
    values(preds) <- boot::inv.logit(values(preds))
    preds <- crop(preds, ecoregions) %>%
      terra::project(ecoregions) %>%
      mask(ecoregions, maskvalues = 1)
    plot(preds)
    
    writeRaster(preds, filename = paste0("./landis_analysis/landis_predictions/full_predictions_", 
                                         species, "_", pred_info$model[j], "_", pred_info$year[j], ".tif"))
    #--------------------
    #climate variables held constant
    clim_stack <- c(predictor_stack_orig[[c(6:14)]], predictor_stack[[-c(10:15)]])
    
    clim_preds <- terra::predict(object = clim_stack,  #TODO fix this like I did below
                                          model = mod,
                                          # ext = st_bbox(bcr_albers),
                                          const = data.frame(time_observations_started = 7,
                                                             duration_minutes = 60)
            )
    
    values(clim_preds) <- boot::inv.logit(values(clim_preds)) 
    
    clim_preds <- clim_preds %>%
      crop(ecoregions) %>%
      terra::project(ecoregions) %>%
      mask(ecoregions, maskvalues = 1)
    
    writeRaster(clim_preds, filename = paste0("./landis_analysis/landis_predictions/clim_predictions_", 
                                              species, "_", pred_info$model[j], "_", pred_info$year[j], ".tif"))
    
    #----------------
    #vegetation variables held constant
    veg_stack <- c(predictor_stack_orig[[c(1,2,3,4,5,18:18:30)]], predictor_stack[[-c(1:9, 19)]])
    veg_preds <- terra::predict(object = veg_stack, 
                                 model = mod,
                                 # ext = st_bbox(bcr_albers),
                                 const = data.frame(time_observations_started = 7,
                                                    duration_minutes = 60)
                                  )
    
    values(veg_preds) <- boot::inv.logit(values(veg_preds)) 
    
    veg_preds <- veg_preds %>%
      crop(ecoregions) %>%
      terra::project(ecoregions) %>%
      mask(ecoregions, maskvalues = 1)
    
    writeRaster(veg_preds, filename = paste0("./landis_analysis/landis_predictions/veg_predictions_", 
                                             species, "_", pred_info$model[j], "_", pred_info$year[j], ".tif"))
  }
}
