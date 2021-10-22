# Prepare spatial data ####
# Markus Bauer



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Preparation ################################################################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


### Packages ###
library(here)
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

sites <- read_csv(here("data/raw/data_raw_sites.csv"), col_names = T, na = "na", col_types = 
                     cols(
                       .default = "?",
                       id = "f",
                       location = "f",
                       side = "f",
                       exposition = "f"
                       )) %>%
  select(id, location, longitude, latitude, constructionYear, sandPerc, phosphorus, phosphorusClass) %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 31468) %>%
  st_transform(4326)

coord <- as_tibble(st_coordinates(sites))
sites2 <- st_drop_geometry(sites) %>%
  mutate(longitude = coord$X) %>%
  mutate(latitude = coord$Y) %>%
  as_tibble()
rm(coord)

blocks <- sites2 %>%
  group_by(location) %>%
  summarise(across(c(longitude, latitude), mean, na.rm = T)) %>%
  rename(longitude_center = longitude, latitude_center = latitude)

sites2 <- left_join(sites2, blocks, by = "location")


## 2 Transform shp files #################################################################################################

setwd(here("data/raw/spatial"))

dikes <- st_read("Deich.shp") %>%
  st_transform(crs = 4326) %>%
  st_crop(ymin = 48.65, ymax = 48.95, xmin = 12.55, xmax = 13.15)

bbox <- st_convex_hull(st_union(dikes))

grazing <- st_read("beweidung_deiche_wwa_deg.shp") %>%
  st_transform(crs = 4326) %>%
  st_intersection(bbox)
  
conservation_area <- st_read("nsg_epsg31468.shp") %>%
  st_transform(crs = 4326) %>% #problem
  st_intersection(bbox)

ffh_area <- st_read("ffh_epsg31468.shp") %>%
  st_transform(crs = 4326) %>% #problem
  st_intersection(bbox)


## 3 Digitize shp files #################################################################################################

#data <- mapview() %>% editMap()
#mapview(data$finished)
#danube <- data$finished %>%
  #st_as_sf() %>%
  #rename(river = X_leaflet_id) %>%
  #mutate(river = str_replace(as.character(river), ".", "Danube")) %>%
  #mutate(river = str_extract(river, "Danube")) %>%
  #st_crop(ymin = 48.65, ymax = 48.95, xmin = 12.55, xmax = 13.15)
#plot(st_geometry(danube))
#rm(data)


## 4 Background map #################################################################################################

germany <- raster::getData('GADM', country = 'DEU', level = 0, download = F) %>%
  st_as_sf() %>%
  st_set_crs(4326)

rivers <- rnaturalearth::ne_download(scale = 10, type = 'rivers_lake_centerlines', category = 'physical') %>% #problem
  st_as_sf() %>%
  st_set_crs(4326) %>%
  st_intersection(bbox)

background_google <- get_map(
  location = c(left = 12.55, bottom = 48.65, right = 13.15, top = 48.95),
  zoom = 10, 
  scale = 1,
  maptype = "terrain",
  source = "google" #problem
)
ggmap(background_google)

background_toner <- get_map(
  location = c(left = 12.55, bottom = 48.65, right = 13.15, top = 48.95),
  zoom = 10, 
  scale = 1,
  maptype = "toner-background", #problem
  source = "stamen"
)
ggmap(background_toner)

background_terrain <- get_map(
  location = c(left = 12.55, bottom = 48.65, right = 13.15, top = 48.95),
  zoom = 10, 
  scale = 1,
  maptype = "terrain-background",
  source = "stamen"
)
ggmap(background_terrain)



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# C Save ##############################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


save(background_google, 
     file = paste0(here("data/processed/spatial"), "/", "background_google.rda"))
save(background_toner, 
     file = paste0(here("data/processed/spatial"), "/", "background_toner.rda"))
save(background_terrain, 
     file = paste0(here("data/processed/spatial"), "/", "background_terrain.rda"))
st_write(germany, layer = "germany_epsg4326.shp", driver = "ESRI Shapefile", delete_layer = T,
         dsn = here("data/processed/spatial"))
st_write(rivers, layer = "rivers_epsg4326.shp", driver = "ESRI Shapefile", delete_layer = T,
         dsn = here("data/processed/spatial"))
st_write(danube, layer = "danubeepsg4326.shp", driver = "ESRI Shapefile", delete_layer = T,
         dsn = here("data/processed/spatial"))
st_write(grazing, layer = "grazing_epsg4326.shp", driver = "ESRI Shapefile", delete_layer = T,
         dsn = here("data/processed/spatial"))
st_write(dikes, layer = "dikes_epsg4326.shp", driver = "ESRI Shapefile", delete_layer = T,
         dsn = here("data/processed/spatial"))
st_write(conservation_area, layer = "conservation_area_epsg4326.shp", driver = "ESRI Shapefile", delete_layer = T,
         dsn = here("data/processed/spatial"))
st_write(ffh_area, layer = "ffh_area_epsg4326.shp", driver = "ESRI Shapefile", delete_layer = T,
         dsn = here("data/processed/spatial"))
st_write(sites, layer = "sites_epsg4326.shp", driver = "ESRI Shapefile", delete_layer = T,
         dsn = here("data/processed/spatial"))
write_csv(sites2, file = here("data/processed/spatial/sites2.csv"))
write_csv(blocks, file = here("data/processed/spatial/blocks.csv"))
setwd(here("data/processed/spatial"))
sites2 %>%
  select(id, lon, lat) %>%
  mutate(id = as.character(id)) %>%
  as.data.frame() %>%
  pgirmess::writeGPX(type = "w", filename = "danube_old_dikes_plots.gpx")
