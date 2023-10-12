# Code for the baseline climate  -----------------------------------------------------------
# Developed by Rodrigo Aguayo (2023)

rm(list=ls())
cat("\014")  

library("foreach")
library("doParallel")
library("exactextractr")
source("elevationZones.R")
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

# glaciers that are modeled
rgi7_hydro      <- vect("GIS South/Glaciers/RGI7_Hydro.shp")

# dem 
dem <- rast("GIS South/dem_patagonia3f.tif")
dem <- aggregate(dem, fact=5, fun="mean")
dem <- disagg(dem, fact=2, method = "bilinear") # avoid bug in elevation bands

# catchments: PMET dataset
catchments_pmet <- vect("GIS South/Basins_PMETobs.shp")
catchments_pmet_ng <- erase(catchments_pmet, rgi7_hydro)

catchments_pmet_eb <- vect()
for(basin in 1:length(catchments_pmet_ng)) {
  
  # get the zones
  basin_ds <- elevationZones(x=catchments_pmet_ng[basin,], dem=dem, max.zones = 5, min.elevZ = 300)
  basin_i  <- as.polygons(basin_ds$zonesRaster, na.rm = TRUE, dissolve = TRUE)
  
  # assign id and zone
  basin_i$gauge_id <- paste0(catchments_pmet_ng[basin,]$gauge_id, "_", basin_ds$table$No.zone)
  catchments_pmet_eb <- rbind(catchments_pmet_eb, basin_i)
}
  
# catchment: all 
catchments_all  <- vect("GIS South/Basins_Patagonia_all.shp")
catchments_all_ng  <- erase(catchments_all, rgi7_hydro) #

# this does not an impact since the glacier runoff is weight by the glacier area
to_recover <- catchments_all$gauge_id[!(catchments_all$gauge_id %in% catchments_all_ng$gauge_id)]
to_recover <- catchments_all[catchments_all$gauge_id %in% to_recover]
catchments_all_ng <- rbind(catchments_all_ng, to_recover)

# Reference climate: PMET -------------------------------------------------------------------------

path <- "/home/rooda/Hydro_results/climate_catchments/"
pp   <- rast("Data/Precipitation/PP_PMETsim_1980_2020_v10d.nc")   
pp   <- pp[[time(pp) <= as.Date("2019-12-31")]]

write.csv(ts_extract(pp, catchments_pmet_ng), paste0(path, "PP_ref_PMET_basins.csv"))
write.csv(ts_extract(pp, catchments_pmet_eb), paste0(path, "PP_ref_PMET_basins_eb.csv"))
write.csv(ts_extract(pp, catchments_all_ng),  paste0(path, "PP_ref_all_basins.csv"))

t2m  <- rast("Data/Temperature/Tavg_PMETsim_1980_2020_v10d.nc", subds = "t2m") 
t2m  <- t2m[[time(t2m) <= as.Date("2019-12-31")]]
write.csv(ts_extract(t2m, catchments_pmet_ng), paste0(path, "T2M_ref_PMET_basins.csv"))
write.csv(ts_extract(t2m, catchments_pmet_eb), paste0(path, "T2M_ref_PMET_basins_eb.csv"))
write.csv(ts_extract(t2m, catchments_all_ng),  paste0(path, "T2M_ref_all_basins.csv"))

pet  <- rast("Data/Evapotranspiration/Ep_PMET_1980_2020d.nc")
pet  <- pet[[time(pet) <= as.Date("2019-12-31")]]
write.csv(ts_extract(pet, catchments_pmet_ng), paste0(path, "PET_ref_PMET_basins.csv"))
write.csv(ts_extract(pet, catchments_pmet_eb), paste0(path, "PET_ref_PMET_basins_eb.csv"))
write.csv(ts_extract(pet, catchments_all_ng),  paste0(path, "PET_ref_all_basins.csv"))

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
  


