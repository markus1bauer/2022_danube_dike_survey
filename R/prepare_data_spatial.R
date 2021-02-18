# Prepare spatial data ####
# Markus Bauer
# Citation: Markus Bauer, Jakob Huber, Johannes Kollmann...

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Preparation ################################################################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

### Packages ###
library(tidyverse)
library(sf)
library(ggmap)

### Start ###
#installr::updateR(browse_news = F, install_R = T, copy_packages = T, copy_Rprofile.site = T, keep_old_packages = T, update_packages = T, start_new_R = F, quit_R = T)
rm(list = ls())
setwd("Z:/Documents/0_Donaudeiche/3_Aufnahmen_und_Ergebnisse/2022_Danube_old_dikes/data/raw")
register_google(key = "AIzaSyB5nQU_dgB_kPsQkk-_cq7pA0g1-2qka4E")


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Load shp files ##########################################################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


## 1 Base map #################################################################################################

ger <- raster::getData('GADM', country = 'DEU', level = 0)
ger <- st_as_sf(ger)
ger <- st_set_crs(ger, 4326)

experiment <- get_map(
  location = c(lon = 12.884, lat = 48.839),
  zoom = 16, 
  scale = 1,
  maptype = "satellite",
  source = "google"
  )
ggmap(experiment)

files <- list.files(pattern = '.gpx$', full.names = T)
allwaypoints <- list()
for (i in 1:length(files)) {
  allwaypoints[[i]] <- plotKML::readGPX(files[i], tracks = F, routes = F)$waypoints[, c('name', 'lon', 'lat')]
}
sites <- do.call('rbind', allwaypoints)
sp::coordinates(sites) <- c("lon","lat")
sites <- st_as_sf(sites)
sites <- st_set_crs(sites, 4326)
sites$dataset <- c(rep("Transects", 40), rep("Blocks", 42))
rm(files, i, allwaypoints)




#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# C Transform data ##############################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#rm(list = setdiff(ls(), c("nsg", "rollfeld", "paths", "parts", "sites", "avp", "avf", "app", "apf")))



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# D Save ##############################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


st_write(ger, layer = "germany.shp", driver = "ESRI Shapefile",
         dsn = "Z:/Documents/0_Donaudeiche/3_Aufnahmen_und_Ergebnisse/2022_Danube_dike_experiment/data/processed/shp_files")
save(experiment, file = "Z:/Documents/0_Donaudeiche/3_Aufnahmen_und_Ergebnisse/2022_Danube_dike_experiment/data/processed/shp_files/map_experiment.rda")
