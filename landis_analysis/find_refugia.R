# compare prediction maps between timesteps
library(terra)

preds_0_cerw <- rast("./landis_analysis/landis_predictions/prediction_cerw_LowTLowV BAU_0.tif") 
plot(preds_0_cerw)
preds_60_cerw <- rast("./landis_analysis/landis_predictions/prediction_cerw_LowTLowV BAU_60.tif") 
plot(preds_60_cerw)
preds_60_veg_cerw <- rast("./landis_analysis/landis_predictions/prediction_cerw_LowTLowV BAU_60 constant_veg.tif")
plot(preds_60_veg_cerw)
preds_60_clim_cerw <- rast("./landis_analysis/landis_predictions/prediction_cerw_LowTLowV BAU_60 constant_clim.tif")
plot(preds_60_clim_cerw)

pred_change_due_to_clim <- preds_60_veg_cerw - preds_0_cerw



preds_0_gwwa <- rast("./landis_analysis/landis_predictions/prediction_gwwa_LowTLowV BAU_0.tif") 
plot(preds_0_gwwa)
preds_60_gwwa <- rast("./landis_analysis/landis_predictions/prediction_gwwa_LowTLowV BAU_60.tif") 
plot(preds_60_gwwa)
preds_60_veg_gwwa <- rast("./landis_analysis/landis_predictions/prediction_gwwa_LowTLowV BAU_60 constant_veg.tif")
plot(preds_60_veg_gwwa)
preds_60_clim_gwwa <- rast("./landis_analysis/landis_predictions/prediction_gwwa_LowTLowV BAU_60 constant_clim.tif")
plot(preds_60_clim_gwwa)

#----------------------------------
# find refugia locations

in_situ_ref_cerw <- preds_0_cerw > 0.5 & preds_60_cerw > 0.5
in_situ_ref_cerw <- as.numeric(in_situ_ref_cerw)
plot(in_situ_ref_cerw)
ex_situ_ref_cerw <- as.numeric(preds_0_cerw < 0.5 & preds_60_cerw > 0.5)
plot(ex_situ_ref_cerw)
ex_situ_ref_cerw[ex_situ_ref_cerw[] == TRUE] <- 2
plot(in_situ_ref_cerw + ex_situ_ref_cerw)

in_situ_ref_gwwa <- as.numeric(preds_0_gwwa > 0.5 & preds_60_gwwa > 0.5)
plot(in_situ_ref_gwwa)
ex_situ_ref_gwwa <- as.numeric(preds_0_gwwa < 0.5 & preds_60_gwwa > 0.5)
plot(ex_situ_ref_gwwa)
ex_situ_ref_gwwa[ex_situ_ref_gwwa[] == TRUE] <- 2
plot(in_situ_ref_gwwa + ex_situ_ref_gwwa)
# 
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

lost_gwwa <- as.numeric(preds_0_gwwa > 0.5 & preds_60_gwwa < 0.5)
lost_gwwa[lost_gwwa[] == TRUE] <- 3
plot(lost_gwwa)

plot(lost_gwwa + in_situ_ref_gwwa + ex_situ_ref_gwwa)


plot(preds_60_gwwa > 0.5)
plot(preds_60_veg_gwwa > 0.5 & preds_60_clim_gwwa > 0.5)

both_good <- preds_60_veg_gwwa > 0.5 & preds_60_clim_gwwa > 0.5


#-----------------------------------
# draft figure
coltb <- data.frame(value=0:3, col=c("lightgray", "#1E6D60", "#1E88E5", "#FF609A"))

ref_gwwa <- in_situ_ref_gwwa + ex_situ_ref_gwwa + lost_gwwa
coltab(ref_gwwa) <- coltb
plot(ref_gwwa)

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
