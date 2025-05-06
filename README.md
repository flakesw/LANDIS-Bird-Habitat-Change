This repository contains code, inputs, and data products for the project, "Climate adaptive forest management for breeding bird habitat in the southern Appalachians", funded by the USGS Southeast Climate Adaptation Science Center. Analysis was performed by Sam Flake. Please contact Sam at [swflake@ncsu.edu](mailto:swflake@ncsu.edu) or [sflake@gmail.com](mailto:sflake@gmail.com) with questions.

This dataset provides projections for future occurrence of several bird species of conservation concern, for the southern Appalachians of SC, GA, NC, TN, VA, and KY. Predictions were created from a LANDIS-II simulation model linked to a species distribution model. SDMs were trained on eBird data and using a suite of topographic, climatic, and vegetation variables, including spaceborne LiDAR. The LANDIS-II model was parameterized using inventory data, soils, topography, and species trait data. 

The LANDIS-II model was run for 80 years using two climate scenarios. Each 10 years, vegetation and climate data was output from the LANDIS-II model and used generate predictions from the SDM. Three sets of predictions were created: a "full" prediction using all of the variables; a "clim" prediction that holds climate constant (i.e., uses year-0 climate) and vegetation data from the present timestep; and "veg" prediction that holds vegetation constant (i.e. year-0 vegetation state) with varying climate. There are five replicate runs for each climate scenario. 

The naming convention is as follows: type of prediction _ species _ climate scenario _ management scenario _ replicate _ timestep

For example, clim_predictions_acfl_HighTHighV BAU - Run1_30 refers to the "clim_predictions" that keep climate constant, for Alder Flycatcher, for the HadGEM2 ES365 RCP 8.5 climate scenario. It is from the first replicate model run (Run1), and data from timestep 30 (i.e., year 2050).

These data are marked with a Creative Common CC0 1.0 Universal License. These data are in the public domain and do not have any use constraints. Users are advised to read the dataset's metadata thoroughly to understand appropriate use and data limitations.
