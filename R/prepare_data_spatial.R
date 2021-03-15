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
register_google(key = "AIzaSyB5nQU_dgB_kPsQkk-_cq7pA0g1-2qka4E")


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Load shp files ##########################################################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

## 1 Sites #################################################################################################

setwd("Z:/Documents/0_Uni/2022_Donaudeiche/3_Aufnahmen_und_Ergebnisse/2022_Danube_old_dikes/data/raw")
sites <- read_csv2("data_raw_sites.csv", col_names = T, na = "na", col_types = 
                     cols(
                       .default = col_double(),
                       id = col_factor(),
                       location = col_factor(),
                       ageCategory = col_factor(),
                       HCl = col_factor(),
                       phosphorousClass = col_factor(),
                       potassiumClass = col_factor(),
                       magnesiumClass = col_character()
                     )        
)
sites <- select(sites, id, location, RW, HW, constructionYear, sand, phosphorous, phosphorousClass)
sites <- st_as_sf(sites, coords = c("RW", "HW"), crs = 31468)
sites <- st_transform(sites, 4326)
coord <- as_tibble(st_coordinates(sites))
sites2 <- st_drop_geometry(sites)
sites2$lon <- coord$X
sites2$lat <- coord$Y
sites2 <- as_tibble(sites2)
blocks <- sites2 %>%
  group_by(location) %>%
  summarise_at(c("lon", "lat", "constructionYear"), mean, na.rm = T) %>%
  rename(lon_cent = lon, lat_cent = lat)
sites2 <- left_join(sites2, blocks, by = "location")
rm(coord)


## 2 Transform shp files #################################################################################################

setwd("Z:/Documents/0_Uni/2022_Donaudeiche/3_Aufnahmen_und_Ergebnisse/2022_Danube_old_dikes/data/raw/shp_files")

data <- st_read("Deich.shp")
data <- st_transform(data, crs = 4326)
dikes <- st_crop(data, ymin = 48.65, ymax = 48.95, xmin = 12.55, xmax = 13.15)
bbox <- st_convex_hull(st_union(dikes))

data <- st_read("beweidung_deiche_wwa_deg.shp")
data <- st_transform(data, crs = 4326)
grazing <- st_intersection(data, bbox)
  
data <- st_read("nsg_epsg31468.shp")
data <- st_transform(data, crs = 4326)
conservation_area <- st_intersection(data, bbox)

data <- st_read("ffh_epsg31468.shp")
data <- st_transform(data, crs = 4326)
ffh_area <- st_intersection(data, bbox)
rm(data, bbox)

## 3 Background map #################################################################################################

ger <- raster::getData('GADM', country = 'DEU', level = 0)
ger <- st_as_sf(ger)
ger <- st_set_crs(ger, 4326)

background_google <- get_map(
  location = c(lon = 12.884, lat = 48.839),
  zoom = 10, 
  scale = 1,
  maptype = "terrain",
  source = "google"
)
ggmap(background_google)
summary(background_google)

background_toner <- get_map(
  location = c(lon = 12.884, lat = 48.839),
  zoom = 10, 
  scale = 1,
  maptype = "toner",
  source = "stamen"
)
ggmap(background_toner)

background_terrain <- get_map(
  location = c(left = 12.55, bottom = 48.65, right = 13.15, top = 48.95),
  zoom = 10, 
  scale = 1,
  maptype = "terrain",
  source = "stamen"
)
ggmap(background_terrain)



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# C Save ##############################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


save(background_google, file = "Z:/Documents/0_Uni/2022_Donaudeiche/3_Aufnahmen_und_Ergebnisse/2022_Danube_old_dikes/data/processed/shp_files/background_google.rda")
save(background_toner, file = "Z:/Documents/0_Uni/2022_Donaudeiche/3_Aufnahmen_und_Ergebnisse/2022_Danube_old_dikes/data/processed/shp_files/background_toner.rda")
save(background_terrain, file = "Z:/Documents/0_Uni/2022_Donaudeiche/3_Aufnahmen_und_Ergebnisse/2022_Danube_old_dikes/data/processed/shp_files/background_terrain.rda")
st_write(ger, layer = "germany.shp", driver = "ESRI Shapefile",
         dsn = "Z:/Documents/0_Uni/2022_Donaudeiche/3_Aufnahmen_und_Ergebnisse/2022_Danube_old_dikes/data/processed/shp_files")
st_write(grazing, layer = "grazing.shp", driver = "ESRI Shapefile",
         dsn = "Z:/Documents/0_Uni/2022_Donaudeiche/3_Aufnahmen_und_Ergebnisse/2022_Danube_old_dikes/data/processed/shp_files")
st_write(dikes, layer = "dikes.shp", driver = "ESRI Shapefile",
         dsn = "Z:/Documents/0_Uni/2022_Donaudeiche/3_Aufnahmen_und_Ergebnisse/2022_Danube_old_dikes/data/processed/shp_files")
st_write(conservation_area, layer = "conservation_area.shp", driver = "ESRI Shapefile",
         dsn = "Z:/Documents/0_Uni/2022_Donaudeiche/3_Aufnahmen_und_Ergebnisse/2022_Danube_old_dikes/data/processed/shp_files")
st_write(ffh_area, layer = "ffh_area.shp", driver = "ESRI Shapefile",
         dsn = "Z:/Documents/0_Uni/2022_Donaudeiche/3_Aufnahmen_und_Ergebnisse/2022_Danube_old_dikes/data/processed/shp_files")
st_write(sites, layer = "sites.shp", driver = "ESRI Shapefile",
         dsn = "Z:/Documents/0_Uni/2022_Donaudeiche/3_Aufnahmen_und_Ergebnisse/2022_Danube_old_dikes/data/processed/shp_files")
setwd("Z:/Documents/0_Uni/2022_Donaudeiche/3_Aufnahmen_und_Ergebnisse/2022_Danube_old_dikes/data/processed/shp_files")
write_csv2(sites2, file = "sites2.csv")
write_csv2(blocks, file = "blocks.csv")
