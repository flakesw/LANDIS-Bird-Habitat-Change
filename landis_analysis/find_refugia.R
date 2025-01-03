#make maps from LANDIS outputs
#full model
library("terra")
library("sf")
library("tidyverse")
library("dismo")
library("gbm")
library("tidyverse")
library("tidyterra")

diverging_color_ramp <- function(ras){
  the_palette_fc <- leaflet::colorNumeric(palette = "RdBu", 
                                          # domain = c(-max(abs(ras[]), na.rm = TRUE), max(abs(ras[]), na.rm = TRUE)),
                                          domain = c(-1, 1),
                                          reverse = FALSE)
  # the_colors <- the_palette_fc(seq((-max(abs(ras[]), na.rm = TRUE)), max(abs(ras[]), na.rm = TRUE), length.out = 50))
  the_colors <- data.frame(values = seq(-1, 1, length.out = 50), col = the_palette_fc(seq(-1, 1, length.out = 50)))
}

states <- sf::st_read("C:/Users/Sam/Documents/Maps/Basic maps/state boundaries/cb_2018_us_state_5m/cb_2018_us_state_5m.shp") %>%
  sf::st_transform(crs = crs(preds_0)) %>%
  st_crop(preds_0)

species <- "gwwa"
model_name <- "LowTLowV BAU nofire"
model_dir <- paste0("D:/SApps LANDIS/Model templates/", model_name, "/")
input_dir <- "D:/SApps LANDIS/Inputs/"
year <- 0

ecoregions <- terra::rast("C:/Users/Sam/Documents/Research/Project-Southern-Appalachians-2018/Models/LANDIS_Sapps_Ecosystem_Papers/Ecos11_NCLD.tif")

sdm <- readRDS(paste0("./landis_analysis/models_for_landis/", species, "_dist_model_full_landis.RDS"))
sdm_clim <- readRDS(paste0("./landis_analysis/models_for_landis/", species, "_dist_model_clim_landis.RDS"))

# predictor_stack <- terra::rast(paste0("./landis_analysis/landis_predictor_layers/pred_stack_", 
#                                       model_name, "_", year, ".tif"))%>%
#   terra::project(ecoregions) %>%
#   mask(ecoregions, maskvalues = 1)
# predictor_stack$understory_ratio <- 1 - predictor_stack$understory_ratio
# 
# preds_60 <- terra::predict(object = predictor_stack,
#                         model = sdm_clim,
#                         # ext = st_bbox(bcr_albers),
#                         const = data.frame(time_observations_started = 7,
#                                            duration_minutes = 60)
# )
# 
# values(preds_60) <- boot::inv.logit(values(preds_60))
# preds_60 <- crop(preds_60, ecoregions) %>%
#   terra::project(ecoregions) %>%
#   mask(ecoregions, maskvalues = 1)
# plot(preds_60)
# writeRaster(preds_60, "./landis_analysis/landis_predictions/prediction_gwwa_LowTLowV BAU_60.tif", overwrite = TRUE)
# 
# 
# 
predictor_stack_init <- terra::rast("./landis_analysis/landis_predictor_layers/pred_stack_study_area.tif")
preds_init <- terra::predict(object = predictor_stack_init,
                             model = sdm_clim,
                             # ext = st_bbox(bcr_albers),
                             const = data.frame(time_observations_started = 7,
                                                duration_minutes = 60)
)

values(preds_init) <- boot::inv.logit(values(preds_init))
preds_init <- crop(preds_init, ecoregions) %>%
  terra::project(ecoregions) %>%
  mask(ecoregions, maskvalues = 1)
plot(preds_init)
preds_0 <- preds_init
# writeRaster(preds_0, "./landis_analysis/landis_predictions/prediction_cerw_LowTLowV BAU_0.tif", overwrite = TRUE)


preds_60_cerw <- rast("./landis_analysis/landis_predictions/prediction_cerw_LowTLowV BAU_60.tif")
preds_60_gwwa <- rast("./landis_analysis/landis_predictions/prediction_gwwa_LowTLowV BAU_60.tif")
preds_60_clim_cerw <- rast("./landis_analysis/landis_predictions/prediction_cerw_LowTLowV BAU_60 clim.tif")
preds_60_gwwa_cerw <- rast("./landis_analysis/landis_predictions/prediction_gwwa_LowTLowV BAU_60 clim.tif")

preds_0_cerw <- rast("./landis_analysis/landis_predictions/prediction_cerw_LowTLowV BAU_0.tif")
preds_0_gwwa <- rast("./landis_analysis/landis_predictions/prediction_gwwa_LowTLowV BAU_0.tif")

diff <- preds_60 - preds_0

# plot(diff, col = diverging_color_ramp(diff)[, 2], breaks = diverging_color_ramp(diff)[, 1])

ggplot() +
  geom_spatraster(data = diff) +
  scale_fill_distiller(palette = "RdBu", limits = c(-1, 1), 
                       direction = 1, na.value = "transparent")
# ggplot() +
#   geom_spatraster(data = diff/preds_0) +
#   scale_fill_distiller(palette = "RdBu", limits = c(-1, 1), 
#                        direction = 1, na.value = "transparent")



#----------------------------------
# find refugia locations

in_situ_ref_cerw <- preds_0_cerw > 0.5 & preds_60_cerw > 0.5
plot(in_situ_ref_cerw)
ex_situ_ref_cerw <- preds_0_cerw < 0.5 & preds_60_cerw > 0.5
plot(ex_situ_ref_cerw)

in_situ_ref_gwwa <- preds_0_gwwa > 0.5 & preds_60_gwwa > 0.5
plot(in_situ_ref_gwwa)
ex_situ_ref_gwwa <- preds_0_gwwa < 0.5 & preds_60_gwwa > 0.5
plot(ex_situ_ref_gwwa)

#potential refugia
pot_ref_cerw <- preds_60_cerw < 0.5 & preds_60_clim_cerw > 0.5
plot(pot_ref_cerw)
pot_ref_gwwa <- preds_60_gwwa < 0.5 & preds_60_clim_gwwa > 0.5
plot(pot_ref_gwwa)

#------------------------------------
#quick draft figure
layout(mat = matrix(c(1,2,3,4,5,6), ncol = 2, byrow = FALSE))


plot(in_situ_ref_cerw)
plot(vect(states), add = TRUE)

plot(ex_situ_ref_cerw)
plot(vect(states), add = TRUE)

plot(pot_ref_cerw)
plot(vect(states), add = TRUE)

plot(in_situ_ref_gwwa)
plot(vect(states), add = TRUE)

plot(ex_situ_ref_gwwa)
plot(vect(states), add = TRUE)

plot(pot_ref_gwwa)
plot(vect(states), add = TRUE)
#amounts of habitat and refugia


#----------------------------------
#characteristics of refugia?

pred_stack_0 <- terra::rast("./landis_analysis/landis_predictor_layers/pred_stack_study_area.tif")

refugia_predictors <- terra::mask(pred_stack_0, ref, maskvalue = FALSE)
non_refugia_predictors <- terra::mask(pred_stack_0, ref, maskvalue = TRUE)
boxplot(refugia_predictors)
boxplot(non_refugia_predictors)

ref_df <- data.frame(refugia_predictors) %>%
  slice_sample(n = 10000) %>%
  mutate(refugia = TRUE)
n_ref_df <- data.frame(non_refugia_predictors) %>%
  slice_sample(n = 10000)%>%
  mutate(refugia = FALSE)
ref_preds_data <- rbind(ref_df, n_ref_df)

test <- glm(refugia ~ ., data = ref_preds_data)
summary(test)

boxplot(biomass ~ refugia, data = ref_preds_data)
vioplot::vioplot(tmmn ~ refugia, data = ref_preds_data, 
                 col = c("violet", "gold"))
