# Code for the baseline climate  -----------------------------------------------------------
# Developed by Rodrigo Aguayo (2023)

rm(list=ls())
cat("\014")  

library("foreach")
library("doParallel")
library("exactextractr")
library("terra")
library("sf")

setwd("/home/rooda/Dropbox/Patagonia/")
terra::gdalCache(30000) # config terra cache for better performance

# function to extract time series for each catchment
ts_extract <- function(stack, shape) {
  ts <- round(t(exact_extract(stack, st_as_sf(shape), "mean", progress = T)), 2)
  colnames(ts) <- shape$gauge_id
  rownames(ts) <- as.character(time(stack))
  return(ts)
}

rgi7_hydro      <- vect("GIS South/Glaciers/RGI7_Hydro.shp")
catchments_pmet <- vect("GIS South/Basins_PMETobs.shp")
catchments_pmet_ng <- erase(catchments_pmet, rgi7_hydro)
catchments_all  <- vect("GIS South/Basins_Patagonia_all.shp")
catchments_all_ng  <- erase(catchments_all, rgi7_hydro) # this erase if 100% is glacier cover 

# this does not an impact since the glacier runoff is weight by the glacier area
to_recover <- catchments_all$gauge_id[!(catchments_all$gauge_id %in% catchments_all_ng$gauge_id)]
to_recover <- catchments_all[catchments_all$gauge_id %in% to_recover]
catchments_all_ng <- rbind(catchments_all_ng, to_recover)

# Reference climate: PMET -------------------------------------------------------------------------

pp   <- rast("Data/Precipitation/PP_PMETsim_1980_2020_v10d.nc")   
pp   <- pp[[time(pp) <= as.Date("2019-12-31")]]

write.csv(ts_extract(pp, catchments_pmet_ng), "/home/rooda/Hydro_results/climate_catchments/PP_ref_PMET_basins.csv")
write.csv(ts_extract(pp, catchments_all_ng), "/home/rooda/Hydro_results/climate_catchments/PP_ref_all_basins.csv")

t2m  <- rast("Patagonia/Data/Temperature/Tavg_PMETsim_1980_2020_v10d.nc", subds = "t2m") 
t2m  <- t2m[[time(t2m) <= as.Date("2019-12-31")]]
write.csv(ts_extract(t2m, catchments_pmet_ng), "/home/rooda/Hydro_results/climate_catchments/T2M_ref_PMET_basins.csv")
write.csv(ts_extract(t2m, catchments_all_ng), "/home/rooda/Hydro_results/climate_catchments/T2M_ref_all_basins.csv")

pet  <- rast("Patagonia/Data/Evapotranspiration/Ep_PMET_1980_2020d.nc")
pet  <- pet[[time(pet) <= as.Date("2019-12-31")]]
write.csv(ts_extract(pet, catchments_pmet_ng), "/home/rooda/Hydro_results/climate_catchments/PET_ref_PMET_basins.csv")
write.csv(ts_extract(pet, catchments_all_ng), "/home/rooda/Hydro_results/climate_catchments/PMET_ref_all_basins.csv")

# Future climate ----------------------------------------------------------------------------------
setwd("/home/rooda/Hydro_results/")

gcms  <- c("CMCC-ESM2", "GFDL-ESM4", "INM-CM5-0", "KACE-1-0-G",  "MIROC6", "MPI-ESM1-2-HR", "MPI-ESM1-2-LR", "MRI-ESM2-0")
gcms  <- c(paste(gcms, "ssp126", sep = "_"), paste(gcms, "ssp585", sep = "_"))

for (gcm in gcms) {
    
  pattern <- paste("PP", gcm, sep = "_")
  pp   <- rast(list.files("future_corrected", pattern, full.names = T))
  names(pp) <- time(pp)
  write.csv(ts_extract(pp, catchments_all_ng), paste("climate_catchments/PP", gcm,  "all_basins.csv", sep = "_"))
  
  pattern <- paste("T2M", gcm, sep = "_")
  t2m  <- rast(list.files("future_corrected", pattern, full.names = T))
  names(t2m) <- time(t2m)
  write.csv(ts_extract(t2m, catchments_all_ng), paste("climate_catchments/T2M", gcm,  "all_basins.csv", sep = "_"))
  
  pattern <- paste("PET", gcm, sep = "_")
  pet  <- rast(list.files("future_corrected", pattern, full.names = T))
  names(pet) <- time(pet)
  write.csv(ts_extract(pet, catchments_all_ng), paste("climate_catchments/PET", gcm,  "all_basins.csv", sep = "_"))
  cat(gcm)
}
  


