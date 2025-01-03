library("terra")
library("tidyverse")
library("sf")

landis_vars <- rast("./landis_analysis/landis_predictor_layers/pred_stack_study_area.tif")
landis_vars

template <- terra::rast("C:/Users/Sam/Documents/Research/Project-Southern-Appalachians-2018/Models/LANDIS_Sapps_Ecosystem_Papers/Ecos11_NCLD.tif")
template <- terra::classify(template, rcl = matrix(data = c(2, 12, 1),
                                             ncol = 3),
                            others = NA,
                            include.lowest = TRUE,
                            right = FALSE)
plot(template)

boundary <- sf::st_read("./maps/blue_ridge.gpkg") %>%
  st_transform(crs(landis_vars))
plot(landis_vars[[1]])
plot(vect(boundary), add = TRUE) 


pred_stack <- rast("./predictor_layers/predictor_stack_bcr28.tif") %>%
  terra::project(landis_vars) %>%
  terra::crop(vect(boundary)) %>%
  terra::mask(template) %>%
  terra::clamp(upper = 300, values = TRUE)

plot(pred_stack[[1]])
plot(landis_vars[[1]])

plot(pred_stack[[1]] - landis_vars[[1]])
hist(pred_stack[[1]] - landis_vars[[1]])
global(pred_stack[[1]] - landis_vars[[1]], mean, na.rm = TRUE)

plot(pred_stack[[2]])
plot(landis_vars[[2]])

plot(pred_stack[[1]] - landis_vars[[1]])
hist(pred_stack[[1]] - landis_vars[[1]])
global(pred_stack[[1]] - landis_vars[[1]], mean, na.rm = TRUE)


