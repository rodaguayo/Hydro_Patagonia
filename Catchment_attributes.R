rm(list=ls())
cat("\014")  

library("exactextractr")
library("hydromad")
library("terra")
library("zoo")
library("sf")

terra::gdalCache(40000)

setwd("/home/rooda/OneDrive/DeepHydro/")
path_pmet <- "/home/rooda/OneDrive/PatagoniaMet/"

# Basins to simulate
basin_shp  <- st_read("GIS/Basins_Patagonia_all.shp")

# Reference period
period <- c(as.POSIXct("1989-12-31"), as.POSIXct("2019/12/31"))

dem <- rast(paste0(path_pmet, "GIS/dem_patagonia3f.tif"))
dem <- subst(dem, NA, 0)  # NAs to sea level (= 0)
slope <- terrain(dem, v='aspect', unit='degrees')

# assign values
basin_shp$elev_mean   <- round(exact_extract(dem,   basin_shp, "mean"),   1)
basin_shp$elev_median <- round(exact_extract(dem,   basin_shp, "median"), 1)
basin_shp$slope_mean  <- round(exact_extract(slope, basin_shp, "median"), 1)

# forest cover
forest_cover <- rast(paste0(path_pmet, "GIS/lc_forest_500m.tif"))
forest_cover <- subst(forest_cover, NA, 0)  # NAs to 0 %
basin_shp$forest_cover <- round(exact_extract(forest_cover, basin_shp, "mean"), 1)

# lake cover
lake_cover <- rast(paste0(path_pmet, "GIS/lc_water_500m.tif"))
lake_cover <- subst(lake_cover, NA, 0)  # NAs to 0 %
basin_shp$lake_cover   <- round(exact_extract(lake_cover, basin_shp, "mean"), 1)

# glacier cover (RGI7)
glaciers  <- vect(paste0(path_pmet, "GIS/Glaciers/RGI7.shp"))
glaciers  <- subset(glaciers, glaciers$cenlat <= -40.5)
glaciers  <- rasterize(glaciers, dem, background = 0) * 100
basin_shp$glacier_cover <- round(exact_extract(glaciers,   basin_shp, "mean"), 1)

# Glacier change (elevation in m y-1)
dh_dt <- rast(paste0(path_pmet, "GIS/Glaciers/dhdt_2021.tif"))
dh_dt <- project(dh_dt, crs(basin_shp))
dh_dt <- mask(dh_dt, vect(basin_shp), overwrite=TRUE)
dh_dt <- crop(dh_dt, vect(basin_shp))*1000 # from m to mm
dh_dt <- dh_dt / 1.091 # from ice to water
dh_dt[is.na(dh_dt)] <- 0

basin_shp$glacier_dhdt <- round(exact_extract(dh_dt, basin_shp, "mean"), 1)

# leaf area index (LAI)
lai_data   <- rast(paste0(path_pmet, "GIS/LAI_climatology_1981-2015.nc4"))
lai_data   <- subst(lai_data, NA, 0)  # NAs to 0 %
lai_data   <- round(exact_extract(lai_data, basin_shp, "mean"), 3)
basin_shp$lai_max  <- apply(lai_data, 1, max, na.rm=TRUE)
basin_shp$lai_diff <- apply(lai_data, 1, max, na.rm=TRUE) - apply(lai_data, 1, min, na.rm=TRUE)

# p_mean from PMET v1.1
pp_stack <- rast(paste0(path_pmet, "Zenodo/v11/PP_PMETsim_1980_2020_v11d.nc"))
terra::time(pp_stack) <- as.POSIXct(time(pp_stack), tz= "UTC") 
pp_stack <- subset(pp_stack, which(time(pp_stack) > period[1] & time(pp_stack) <= period[2]))
pp_stack <- mean(tapp(pp_stack, strftime(time(pp_stack), format="%Y"), fun = sum, na.rm = TRUE))
basin_shp$p_mean_PMET <- round(exact_extract(pp_stack, basin_shp, "mean"), 1)

# pet_mean from PMET v1.1 (Hargreaves eq based on temp)
pet_stack <- rast(paste0(path_pmet, "Data/Evaporation/Ep_PMETsim_1980_2020d_dev.nc"))
terra::time(pet_stack) <- as.POSIXct(time(pet_stack), tz= "UTC") 
pet_stack <- subset(pet_stack, which(time(pet_stack) > period[1] & time(pet_stack) <= period[2]))
pet_stack <- mean(tapp(pet_stack, strftime(time(pet_stack),format="%Y"), fun = sum, na.rm = TRUE))
basin_shp$pet_mean_PMET <- round(exact_extract(pet_stack, basin_shp, "mean"), 1)

# aridity using p_mean and pet_mean
basin_shp$aridity_PMET <- round(basin_shp$p_mean_PMET/basin_shp$pet_mean_PMET, 3)

# precipitation data for each basin
pp_stack <- rast(paste0(path_pmet, "Zenodo/v11/PP_PMETsim_1980_2020_v11d.nc"))
terra::time(pp_stack) <- as.POSIXct(time(pp_stack), tz= "UTC") 
pp_stack <- subset(pp_stack, which(time(pp_stack) > period[1] & time(pp_stack) <= period[2]))

time_pp  <- time(pp_stack)
pp_stack <- t(exact_extract(pp_stack, basin_shp, "mean"))
row.names(pp_stack) <- as.character(time_pp)
colnames(pp_stack)  <- basin_shp$gauge_id
pp_stack <- zoo(pp_stack,  order.by = time_pp)

# temperature data for each basin
t2m_stack <- rast(paste0(path_pmet, "Zenodo/v11/Tavg_PMETsim_1980_2020_v11d.nc"))
terra::time(t2m_stack) <- as.POSIXct(time(t2m_stack), tz= "UTC") 
t2m_stack <- subset(t2m_stack, which(time(t2m_stack) > period[1] & time(t2m_stack) <= period[2]))

time_t2m  <- time(t2m_stack)
t2m_stack <- t(exact_extract(t2m_stack, basin_shp, "mean"))
row.names(t2m_stack) <- as.character(time_t2m)
colnames(t2m_stack)  <- basin_shp$gauge_id
t2m_stack <- zoo(t2m_stack,  order.by = time_t2m)

# several precipitation and temperature metrics
for (i in 1:ncol(pp_stack)) {
  metric_i <- eventseq(pp_stack[,i], thresh = 5*mean(pp_stack[,i]))
  metric_i <- eventinfo(pp_stack[,i], metric_i)
  
  # high_prec_freq (≥ 5 times mean daily precipitation) in days yr−1 
  basin_shp$high_prec_freq_PMET[i] <- mean(aggregate(metric_i$Duration, by = list(metric_i$Year), FUN = sum)$x)
  
  # high_prec_dur (number of consecutive days ≥ 5 times mean daily pp)
  basin_shp$high_prec_dur_PMET[i] <- mean(metric_i$Duration)
  
  metric_i <- eventseq(pp_stack[,i], thresh = 1, below = TRUE)
  metric_i <- eventinfo(pp_stack[,i], metric_i)
  
  # low_prec_freq (frequency of dry days (< 1 mm day−1) )
  basin_shp$low_prec_freq_PMET[i] <- mean(aggregate(metric_i$Duration, by = list(metric_i$Year), FUN = sum)$x)
  
  # low_prec_dur (frequency of dry days (< 1 mm day−1) )
  basin_shp$low_prec_dur_PMET[i] <- mean(metric_i$Duration)
  print(i)
}

# p_frac_snow (days colder than 0◦C)
pp_total <- colMeans(aggregate(pp_stack, by = list(substr(row.names(pp_stack),1,4)), FUN = sum))
pp_stack <- as.data.frame(pp_stack)
pp_stack[as.data.frame(t2m_stack) > 0] <- 0
pp_snow  <- aggregate(pp_stack, by = list(substr(row.names(pp_stack),1,4)), FUN = sum)
pp_snow  <- colMeans(pp_snow[,-1])
basin_shp$frac_snow_PMET     <- round(pp_snow/pp_total, 3)

writeVector(vect(basin_shp), "GIS/Basins_Patagonia_all_data.shp", overwrite=TRUE)
write.csv(st_drop_geometry(basin_shp), "Attributes_all_basins.csv", row.names = FALSE)
