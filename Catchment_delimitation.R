rm(list=ls())
cat("\014")  

library("stringr")
library("rgrass") #TODO: This could parallized using https://htmlpreview.github.io/?https://github.com/petrasovaa/FUTURES-CONUS-talk/blob/main/foss4g2022.html#/
library("terra")

setwd("/home/rooda/Dropbox/Patagonia/")

# DEM
dem        <- rast("GIS South/dem_patagonia3f.tif")
dem        <- aggregate(dem, fact=2, fun="mean")
depre      <- is.na(dem) # depression DEM -> water/ocean -> use less RAM

# Station locations 
q_location <- read.csv("Data/Streamflow/Q_PMETobs_v10_metadata.csv")
q_vect     <- vect(q_location, geom=c("gauge_lat", "gauge_lon"), crs="epsg:4326")

# Initializate GRASS and files
initGRASS(gisBase = "/usr/lib/grass82/", home="/home/rooda/Dropbox/Patagonia/", SG = dem, override = TRUE)

write_RAST(dem,    vname = "dem_grass",   flags = c("overwrite", "o"))
write_RAST(depre,  vname = "depre_grass", flags = c("overwrite", "o"))
write_VECT(q_vect, vname = "q_grass",     flags = c("overwrite", "o"))

# based on stream modules: https://grasswiki.osgeo.org/wiki/R.stream.*_modules
execGRASS("r.stream.extract", flags=c("overwrite"), 
          parameters = list(elevation="dem_grass", depression = "depre_grass", threshold = 200, 
                            direction = "fdir", stream_vector="stream_v", stream_raster="stream_r"))

execGRASS("r.stream.basins", flags=c("overwrite", "l"),
          parameters=list(direction="fdir", stream_rast = "stream_r", basins="basins"))


all_basins      <- read_RAST("basins", NODATA = -2147483648)
all_basins      <- as.polygons(all_basins)
all_basins$total_area <- round(expanse(all_basins, unit="km"), 2)
all_basins$lat      <- as.data.frame(geom(centroids(all_basins)))$y
all_basins$lon      <- as.data.frame(geom(centroids(all_basins)))$x

# only the study area (western side)
all_basins <- subset(all_basins, all_basins$lat < -41)
all_basins <- subset(all_basins, all_basins$lon < -68.5)
all_basins <- subset(all_basins, (all_basins$lon < -70) | (all_basins$lat < -54))
all_basins <- subset(all_basins, !((all_basins$lon > -71.5) & (all_basins$lat > -52) & (all_basins$lat < -48))) 
all_basins <- subset(all_basins, !((all_basins$lon > -70) & (all_basins$lat > -54.2) & (all_basins$lat < -53))) 

# small islands outside pmet coverage (n < 10)
all_basins <- subset(all_basins, !((all_basins$lon < -70.4) & (all_basins$lat > -55.2) & (all_basins$lat < -54.95))) 
all_basins <- subset(all_basins, !((all_basins$lon < -74.4) & (all_basins$lat > -53.0) & (all_basins$lat < -52.7))) 
all_basins <- subset(all_basins, !((all_basins$lon < -71.7) & (all_basins$lat > -54.7) & (all_basins$lat < -54.4))) 
all_basins <- subset(all_basins, !((all_basins$lon < -71.1) & (all_basins$lat > -55.0) & (all_basins$lat < -54.7))) 

# subset: area > 5 km2 (4000 -> 3800)
all_basins <- subset(all_basins, all_basins$total_area > 5) # this value could be lower
all_basins

# save file
all_basins$gauge_id <- paste0("Y", sprintf("%08d", seq(1, length(all_basins))))
all_basins <- all_basins[,c("gauge_id", "total_area","lat","lon")]
writeVector(all_basins, "GIS South/Basins_Patagonia_all.shp", overwrite=TRUE)

