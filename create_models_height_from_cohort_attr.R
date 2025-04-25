#----------------------------------------------------------
# estimate cohort attributes from FIA

# this just needs to be run once, to create the regression inputs we use later to 
# translate LANDIS outputs to the SDM data


library("tidyverse")
library("rFIA")

input_dir <- "D:/SApps LANDIS/"


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
hardwood_regression <- lm(log(HT) ~ log(CohortAge), data = sitetrees[sitetrees$SFTWD_HRDWD == "H", ])
softwood_regression <- lm(log(HT) ~ log(CohortAge), data = sitetrees[sitetrees$SFTWD_HRDWD == "S", ])

tree_regressions <- tibble(SPCD = c(1,2),
                           model = list(hardwood_regression, softwood_regression)) %>%
  bind_rows(tree_regressions, .)

saveRDS(tree_regressions, file = "./landis_analysis/landis_predictor_layers/tree_height_regressions.RDS")
