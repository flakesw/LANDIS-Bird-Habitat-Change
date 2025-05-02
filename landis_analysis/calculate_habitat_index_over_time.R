#get habitat index for each raster

#just one rep for demonstration purposes

library(terra)
library(tidyverse)



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
  filter(replicate == "Run1")

pred_info2 <- pred_info %>%
  mutate(mean_HI = global(rast(pred_loc), mean, na.rm = TRUE)$mean) %>%
  mutate(year = as.numeric(year))
pred_info3 <- pred_info %>%
  mutate(high_quality_habitat = global(rast(pred_loc) > 0.75, mean, na.rm = TRUE)$mean) %>%
  mutate(year = as.numeric(year))

test <- rast(pred_info[192, 1][[1]])
plot(test)

theme_set(theme_light())
theme_set(theme(strip.background = element_blank()))
          
ggplot(pred_info2 %>% filter(type == "full"), aes(x = year, y = mean_HI, 
                                                  color = climate, shape = climate, 
                                                  linetype = climate, group = climate)) +
  geom_line() +
  geom_point() +
  facet_wrap(facets = "species", labeller = as_labeller(toupper)) +
  ylab("Mean Habitat Index") +
  xlab("Timestep")+
  scale_x_continuous(breaks=scales::pretty_breaks(n=5))

ggplot(pred_info2 %>% filter(climate == "HighTHighV"), aes(x = year, y = mean_HI, color = type, shape = type, linetype = type )) +
  geom_line() +
  geom_point() +
  facet_wrap(facets = "species", labeller = as_labeller(toupper))+
  ylab("Mean Habitat Index") +
  xlab("Timestep")+
  scale_x_continuous(breaks=scales::pretty_breaks(n=5))
# group = interaction(type, climate) ##save this trick for later

ggplot(pred_info3, aes(x = year, y = high_quality_habitat, color = type, shape = climate, group = interaction(type, climate))) +
  geom_line() +
  geom_point() +
  facet_wrap(facets = "species", labeller = as_labeller(toupper))



pred_info4 <- pred_info2 %>%
  filter(climate == "HighTHighV") %>%
  group_by(species, year) %>%
  summarize(diff_due_to_clim_change = mean_HI[type == "full"] - mean_HI[type == "clim"],
            diff_due_to_veg_change = mean_HI[type == "full"] - mean_HI[type == "veg"])

ggplot(pred_info4, aes(x = year, y = diff_due_to_veg_change)) +
  geom_line(aes(group = "species")) +
  geom_point() +
  facet_wrap(facets = "species") +
  geom_hline(yintercept = 0) + 
  geom_line(aes(group = "species", y = diff_due_to_clim_change)) + 
  geom_point(aes(group = "species", y = diff_due_to_clim_change))
  

