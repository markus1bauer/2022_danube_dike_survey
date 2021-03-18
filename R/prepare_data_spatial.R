# Prepare spatial data ####
# Markus Bauer


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Preparation ################################################################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

### Packages ###
library(tidyverse)
library(sf)
library(ggmap)
library(mapview)
library(mapedit)

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
                       side = col_factor(),
                       exposition = col_factor(),
                       ageCategory = col_factor(),
                       HCl = col_factor(),
                       humusLevel = col_factor(),
                       cnLevel = col_factor(),
                       phosphorousClass = col_factor(),
                       potassiumClass = col_factor(),
                       magnesiumClass = col_character()
                     )) %>%
  select(id, location, RW, HW, constructionYear, sandPerc, phosphorous, phosphorousClass) %>%
  st_as_sf(coords = c("RW", "HW"), crs = 31468) %>%
  st_transform(4326)

coord <- as_tibble(st_coordinates(sites))
sites2 <- st_drop_geometry(sites) %>%
  mutate(lon = coord$X) %>%
  mutate(lat = coord$Y) %>%
  as_tibble()
rm(coord)

blocks <- sites2 %>%
  group_by(location) %>%
  summarise_at(c("lon", "lat", "constructionYear"), mean, na.rm = T) %>%
  rename(lon_cent = lon, lat_cent = lat)

sites2 <- left_join(sites2, blocks, by = "location")


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


## 3 Digitize shp files #################################################################################################

data <- mapview() %>% editMap()
mapview(data$finished)
danube <- data$finished %>%
  st_as_sf() %>%
  rename(river = X_leaflet_id) %>%
  mutate(river = str_replace(as.character(river), ".", "Danube")) %>%
  mutate(river = str_extract(river, "Danube")) %>%
  st_crop(ymin = 48.65, ymax = 48.95, xmin = 12.55, xmax = 13.15)
plot(st_geometry(danube))
rm(data)


## 4 Background map #################################################################################################

data <- raster::getData('GADM', country = 'DEU', level = 0, download = F)
data <- st_as_sf(data)
germany <- st_set_crs(data, 4326)

data <- rnaturalearth::ne_download(scale = 10, type = 'rivers_lake_centerlines', category = 'physical')
data <- st_as_sf(data)
data <- st_set_crs(data, 4326)
rivers <- st_intersection(data, bbox)

background_google <- get_map(
  location = c(lon = 12.884, lat = 48.839),
  zoom = 10, 
  scale = 1,
  maptype = "terrain",
  source = "google"
)
ggmap(background_google)

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


save(background_google, file = "Z:/Documents/0_Uni/2022_Donaudeiche/3_Aufnahmen_und_Ergebnisse/2022_Danube_old_dikes/data/processed/spatial/background_google.rda")
save(background_toner, file = "Z:/Documents/0_Uni/2022_Donaudeiche/3_Aufnahmen_und_Ergebnisse/2022_Danube_old_dikes/data/processed/spatial/background_toner.rda")
save(background_terrain, file = "Z:/Documents/0_Uni/2022_Donaudeiche/3_Aufnahmen_und_Ergebnisse/2022_Danube_old_dikes/data/processed/spatial/background_terrain.rda")
st_write(germany, layer = "germany.shp", driver = "ESRI Shapefile", delete_layer = T,
         dsn = "Z:/Documents/0_Uni/2022_Donaudeiche/3_Aufnahmen_und_Ergebnisse/2022_Danube_old_dikes/data/processed/spatial")
st_write(rivers, layer = "rivers.shp", driver = "ESRI Shapefile", delete_layer = T,
         dsn = "Z:/Documents/0_Uni/2022_Donaudeiche/3_Aufnahmen_und_Ergebnisse/2022_Danube_old_dikes/data/processed/spatial")
st_write(danube, layer = "danube.shp", driver = "ESRI Shapefile", delete_layer = T,
         dsn = "Z:/Documents/0_Uni/2022_Donaudeiche/3_Aufnahmen_und_Ergebnisse/2022_Danube_old_dikes/data/processed/spatial")
st_write(grazing, layer = "grazing.shp", driver = "ESRI Shapefile", delete_layer = T,
         dsn = "Z:/Documents/0_Uni/2022_Donaudeiche/3_Aufnahmen_und_Ergebnisse/2022_Danube_old_dikes/data/processed/spatial")
st_write(dikes, layer = "dikes.shp", driver = "ESRI Shapefile", delete_layer = T,
         dsn = "Z:/Documents/0_Uni/2022_Donaudeiche/3_Aufnahmen_und_Ergebnisse/2022_Danube_old_dikes/data/processed/spatial")
st_write(conservation_area, layer = "conservation_area.shp", driver = "ESRI Shapefile", delete_layer = T,
         dsn = "Z:/Documents/0_Uni/2022_Donaudeiche/3_Aufnahmen_und_Ergebnisse/2022_Danube_old_dikes/data/processed/spatial")
st_write(ffh_area, layer = "ffh_area.shp", driver = "ESRI Shapefile", delete_layer = T,
         dsn = "Z:/Documents/0_Uni/2022_Donaudeiche/3_Aufnahmen_und_Ergebnisse/2022_Danube_old_dikes/data/processed/spatial")
st_write(sites, layer = "sites.shp", driver = "ESRI Shapefile", delete_layer = T,
         dsn = "Z:/Documents/0_Uni/2022_Donaudeiche/3_Aufnahmen_und_Ergebnisse/2022_Danube_old_dikes/data/processed/spatial")
setwd("Z:/Documents/0_Uni/2022_Donaudeiche/3_Aufnahmen_und_Ergebnisse/2022_Danube_old_dikes/data/processed/spatial")
write_csv2(sites2, file = "sites2.csv")
write_csv2(blocks, file = "blocks.csv")
