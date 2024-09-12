# Catchment selection  --------------------------------------------------------------

rm(list=ls())
cat("\014")  

library("terra")

setwd("/home/rooda/OneDrive/DeepHydro/")
path_pmet <- "/home/rooda/OneDrive/PatagoniaMet/"

# threshold
threshold <- 0.75

# time period
period_calib    <- seq(as.POSIXct("1990-01-01", tz = "UTC"), 
                       as.POSIXct("2019-12-31", tz = "UTC"), by = "day")

# Basin data
q_metadata <- read.csv(paste0(path_pmet, "Zenodo/v11/Q_PMETobs_v11_metadata.csv"), row.names = "gauge_id")
q_metadata$gauge_code <- as.character(q_metadata$gauge_code)
q_obs      <- read.csv(paste0(path_pmet, "Zenodo/v11/Q_PMETobs_1950_2020_v11d.csv"), row.names = "Date")
q_obs      <- q_obs[rownames(q_obs) %in% as.character(period_calib), ] 
q_shape    <- vect(paste0(path_pmet, "GIS/Basins_PMETobs_points_dev.shp"))

# subset
basin_select <- colSums(!is.na(q_obs))/nrow(q_obs)
q_metadata   <- q_metadata[(basin_select > threshold),] 
q_shape   <- q_shape[(basin_select > threshold),] 

# TODO. Remove Futaleufu river

# save
writeVector(q_shape, "GIS/Basins_PMETobs_points_subset.shp", overwrite = TRUE)
