# Historical and future climate for all catchments  -----------------------------------------------------------

rm(list=ls())
cat("\014")  

library("foreach")
library("doParallel")
library("exactextractr")
source("elevationZones.R")
library("terra")
library("sf")

# config terra cache for better performance
terra::gdalCache(30000) 

setwd("/home/rooda/OneDrive/DeepHydro/")
path_pmet <- "/home/rooda/OneDrive/PatagoniaMet/"
path_disk <- "/home/rooda/Hydro_results/climate_catchments/"
path_future <- "/home/rooda/Hydro_results/future_corrected/"

# function to extract time series for each catchment
ts_extract <- function(stack, shape) {
  ts <- round(t(exact_extract(stack, st_as_sf(shape), "mean", progress = T, max_cells_in_memory = 3e+08)), 2)
  colnames(ts) <- shape$gauge_id
  rownames(ts) <- as.character(time(stack))
  return(ts)
}

# glaciers that are modeled
rgi7_hydro      <- vect("GIS/RGI7_Hydro.shp")

# dem 
dem <- rast(paste0(path_pmet, "GIS/dem_patagonia1.tif"))
dem <- subst(dem, NA, 0)
dem <- aggregate(dem, fact=3, fun="mean", cores = 20)
dem <- dem + init(dem, runif) # add random noise

# catchments: PMET dataset 
catchments_pmet <- vect(paste0(path_pmet, "GIS/Basins_PMETobs_dev.shp"))
catchments_pmet_ng <- erase(catchments_pmet, rgi7_hydro)

catchments_pmet_eb <- vect()
for(basin in 1:length(catchments_pmet_ng)) {
  
  # get the zones
  basin_ds <- elevationZones(x=catchments_pmet_ng[basin,], dem=dem, max.zones = 10, min.elevZ = 300)
  basin_i  <- as.polygons(basin_ds$zonesRaster, na.rm = TRUE, dissolve = TRUE)
  
  # assign id and zone
  basin_i$gauge_id <- paste0(catchments_pmet_ng[basin,]$gauge_id, "_", basin_ds$table$No.zone)
  catchments_pmet_eb <- rbind(catchments_pmet_eb, basin_i)
  print(basin)
}
  
# catchment: all 
catchments_all  <- vect("GIS/Basins_Patagonia_all.shp")
catchments_all_ng  <- erase(catchments_all, rgi7_hydro) 

## TODO this does not have an impact as the glacier runoff is weighted by the glacier area (??)
to_recover <- catchments_all$gauge_id[!(catchments_all$gauge_id %in% catchments_all_ng$gauge_id)]
to_recover <- catchments_all[catchments_all$gauge_id %in% to_recover]
catchments_all_ng <- rbind(catchments_all_ng, to_recover) # only 2

catchments_all_eb <- vect()
for(basin in 1:length(catchments_all_ng)) {
  
  # get the zones
  basin_ds <- elevationZones(x=catchments_all_ng[basin,], dem=dem, max.zones = 10, min.elevZ = 300)
  basin_i  <- as.polygons(basin_ds$zonesRaster, na.rm = TRUE, dissolve = TRUE)
  
  # assign id and zone
  basin_i$gauge_id <- paste0(catchments_all_ng[basin,]$gauge_id, "_", basin_ds$table$No.zone)
  catchments_all_eb <- rbind(catchments_all_eb, basin_i)
  print(basin)
}


# Baseline climate: PMET -------------------------------------------------------------------------

pp   <- rast(paste0(path_pmet, "Zenodo/v11/PP_PMETsim_1980_2020_v11d.nc"))
pp   <- pp[[as.Date("1987-01-01") <= time(pp) &  time(pp) <= as.Date("2019-12-31")]]

write.csv(ts_extract(pp, catchments_pmet),    paste0(path_disk, "PP_ref_PMET_basins_full.csv"))
write.csv(ts_extract(pp, catchments_pmet_ng), paste0(path_disk, "PP_ref_PMET_basins_ng.csv"))
write.csv(ts_extract(pp, catchments_pmet_eb), paste0(path_disk, "PP_ref_PMET_basins_eb.csv"))
write.csv(ts_extract(pp, catchments_all),     paste0(path_disk, "PP_ref_all_basins_full.csv"))
write.csv(ts_extract(pp, catchments_all_ng),  paste0(path_disk, "PP_ref_all_basins_ng.csv"))
write.csv(ts_extract(pp, catchments_all_eb),  paste0(path_disk, "PP_ref_all_basins_eb.csv"))

t2m  <- rast(paste0(path_pmet, "Zenodo/v11/Tavg_PMETsim_1980_2020_v11d.nc"))
t2m   <- t2m[[as.Date("1987-01-01") <= time(t2m) &  time(t2m) <= as.Date("2019-12-31")]]

write.csv(ts_extract(t2m, catchments_pmet),    paste0(path_disk, "T2M_ref_PMET_basins_full.csv"))
write.csv(ts_extract(t2m, catchments_pmet_ng), paste0(path_disk, "T2M_ref_PMET_basins_ng.csv"))
write.csv(ts_extract(t2m, catchments_pmet_eb), paste0(path_disk, "T2M_ref_PMET_basins_eb.csv"))
write.csv(ts_extract(t2m, catchments_all),     paste0(path_disk, "T2M_ref_all_basins_full.csv"))
write.csv(ts_extract(t2m, catchments_all_ng),  paste0(path_disk, "T2M_ref_all_basins_ng.csv"))
write.csv(ts_extract(t2m, catchments_all_eb),  paste0(path_disk, "T2M_ref_all_basins_eb.csv"))

pet  <- rast(paste0(path_pmet, "Data/Evaporation/Ep_PMETsim_1980_2020d_dev.nc"))
pet  <- pet[[as.Date("1987-01-01") <= time(pet) &  time(pet) <= as.Date("2019-12-31")]]

write.csv(ts_extract(pet, catchments_pmet),    paste0(path_disk, "PET_ref_PMET_basins_full.csv"))
write.csv(ts_extract(pet, catchments_pmet_ng), paste0(path_disk, "PET_ref_PMET_basins_ng.csv"))
write.csv(ts_extract(pet, catchments_pmet_eb), paste0(path_disk, "PET_ref_PMET_basins_eb.csv"))
write.csv(ts_extract(pet, catchments_all),     paste0(path_disk, "PET_ref_all_basins_full.csv"))
write.csv(ts_extract(pet, catchments_all_ng),  paste0(path_disk, "PET_ref_all_basins_ng.csv"))
write.csv(ts_extract(pet, catchments_all_eb),  paste0(path_disk, "PET_ref_all_basins_eb.csv"))


# Climate projections ----------------------------------------------------------------------------------

gcms  <- c("GFDL-ESM4", "IPSL-CM6A-LR", "MIROC6", "MPI-ESM1-2-LR", "MRI-ESM2-0")
gcms  <- c(paste(gcms, "ssp126", sep = "_"), paste(gcms, "ssp585", sep = "_"))

for (gcm in gcms) {
    
  pattern <- paste("PP", gcm, sep = "_")
  pp   <- rast(list.files(path_future, pattern, full.names = T))
  names(pp) <- time(pp)
  write.csv(ts_extract(pp, catchments_all),    paste(paste0(path_disk, "PP"), gcm,  "all_basins_full.csv", sep = "_"))
  write.csv(ts_extract(pp, catchments_all_ng), paste(paste0(path_disk, "PP"), gcm,  "all_basins_ng.csv", sep = "_"))
  write.csv(ts_extract(pp, catchments_all_eb), paste(paste0(path_disk, "PP"), gcm,  "all_basins_eb.csv", sep = "_"))
  
  pattern <- paste("T2M", gcm, sep = "_")
  t2m  <- rast(list.files(path_future, pattern, full.names = T))
  names(t2m) <- time(t2m)
  write.csv(ts_extract(t2m, catchments_all),    paste(paste0(path_disk, "T2M"), gcm,  "all_basins_full.csv", sep = "_"))
  write.csv(ts_extract(t2m, catchments_all_ng), paste(paste0(path_disk, "T2M"), gcm,  "all_basins_ng.csv", sep = "_"))
  write.csv(ts_extract(t2m, catchments_all_eb), paste(paste0(path_disk, "T2M"), gcm,  "all_basins_eb.csv", sep = "_"))
  
  pattern <- paste("PET", gcm, sep = "_")
  pet  <- rast(list.files(path_future, pattern, full.names = T))
  names(pet) <- time(pet)
  write.csv(ts_extract(pet, catchments_all),    paste(paste0(path_disk, "PET"),  "all_basins_full.csv", sep = "_"))
  write.csv(ts_extract(pet, catchments_all_ng), paste(paste0(path_disk, "PET"),  "all_basins_ng.csv", sep = "_"))
  write.csv(ts_extract(pet, catchments_all_eb),   paste(paste0(path_disk, "PET"), gcm,  "all_basins_eb.csv", sep = "_"))
  
  cat(gcm)
}


