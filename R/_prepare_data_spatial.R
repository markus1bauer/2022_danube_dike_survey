# Beta diversity on dike grasslands
# Prepare spatial data ####
# Markus Bauer
# 2024-02-07



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Preparation ###############################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



### Packages ###
library(here)
library(sf)
library(tidyverse)
library(ggmap)
library(mapview)
library(mapedit)

### Start ###
rm(list = ls())
register_google(key = "AIzaSyB5nQU_dgB_kPsQkk-_cq7pA0g1-2qka4E")



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Load shp files ############################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



#______________________________________________________________________________
## 1 Sites ####################################################################


sites <- read_csv(here("data", "raw", "data_raw_sites.csv"),
  col_names = TRUE,
  na = "na", col_types =
    cols(
      .default = "?",
      id = "f",
      location = "f"
    )) %>%
  select(id, location, longitude, latitude, construction_year) %>%
  mutate(
    plot = str_sub(id, start = 1, end = 2),
    location_abb = str_sub(location, 1, 3),
    location_abb = str_to_upper(location_abb),
    location_abb = factor(
      location_abb,
      levels = unique(location_abb[order(construction_year)])
      ),
    location_construction_year = str_c(location_abb, construction_year, sep = "-")
  ) %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 31468) %>%
  st_transform(4326)

### Create tbl instead of sf file ###
coord <- as_tibble(st_coordinates(sites))
sites_with_spatial_data <- sites %>%
  st_drop_geometry() %>%
  mutate(longitude = coord$X) %>%
  mutate(latitude = coord$Y) %>%
  as_tibble()
rm(coord)



#______________________________________________________________________________
## 2 Transform shp files ######################################################


dikes <- st_read(here("data", "raw", "spatial", "dikes_epsg31468.shp")) %>%
  st_transform(crs = 4326) %>%
  st_crop(ymin = 48.65, ymax = 48.95, xmin = 12.55, xmax = 13.15)

bbox <- st_convex_hull(st_union(dikes))

rivers <- st_read(here("data", "raw", "spatial", "rivers_epsg25832.shp")) %>%
  st_transform(crs = 4326) %>%
  st_intersection(bbox) %>%
  filter(GEWKZ_S == "1" | GEWKZ_S == "2" | GEWKZ_S == "3" | GEWKZ_S == "4")

grazing <- st_read(here("data", "raw", "spatial", "grazing_epsg31468.shp")) %>%
  st_transform(crs = 4326) %>%
  st_intersection(bbox)

conservation_area <- st_read(
  here("data", "raw", "spatial", "conservation_area_epsg31468.shp")
  ) %>%
  st_transform(crs = 4326)  %>%
  st_make_valid() %>%
  st_intersection(bbox)

ffh_area <- st_read(here("data", "raw", "spatial", "ffh_epsg31468.shp")) %>%
  st_transform(crs = 4326)  %>%
  st_make_valid() %>%
  st_intersection(bbox)

biotope_mapping <- st_read(
  here("data", "raw", "spatial", "bio_fbk_epsg25832_shp.shp")
  ) %>%
  st_transform(crs = 4326)  %>%
  st_intersection(bbox)



#______________________________________________________________________________
## 3 Background map ###########################################################


germany <- geodata::gadm(
  country = "DEU", level = 0, version = "latest", resolution = 2, path = here()
  ) %>%
  st_as_sf() %>%
  st_set_crs(4326)

# problem since 2021
# https://github.com/ropensci/rnaturalearth/issues/29
# https://github.com/nvkelso/natural-earth-vector/issues/528
#rivers <- rnaturalearth::ne_download(
#  scale = 10,
#  type = "rivers_lake_centerlines",
#  category = "physical",
#  load = FALSE,
#  returnclass = "sf"
#) %>%
#  st_set_crs(4326) %>%
#  st_intersection(bbox)

# problem
# no valid key any longer
#background_google <- get_map(
#  location = c(left = 12.55, bottom = 48.65, right = 13.15, top = 48.95),
#  zoom = 10,
#  scale = 1,
#  maptype = "terrain",
#  source = "google" 
#)
#ggmap(background_google)

background_toner <- get_map(
  location = c(left = 12.55, bottom = 48.65, right = 13.15, top = 48.95),
  zoom = 10,
  scale = 1,
  maptype = "toner-background",
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



#______________________________________________________________________________
## 4 Digitize shp files #######################################################


# Was only done once

# data <- mapview() %>% editMap()
# mapview(data$finished)
# danube_isar <- data$finished %>%
# st_as_sf() %>%
# rename(river = X_leaflet_id) %>%
# mutate(river = str_replace(as.character(river), ".", "Danube")) %>%
# mutate(river = str_extract(river, "Danube")) %>%
# st_crop(ymin = 48.65, ymax = 48.95, xmin = 12.55, xmax = 13.15)
# plot(st_geometry(danube_isar))
# rm(data)
### Here the digitized file ###
danube_isar <- st_read(
  here("data", "raw", "spatial", "danube_isar_digitized_epsg4326.shp")
  )



#______________________________________________________________________________
## 5 Calculate new variables ##################################################


### a Calculate center of locations -------------------------------------------

locations <- sites_with_spatial_data %>%
  group_by(location_construction_year) %>%
  summarise(
    across(c(longitude, latitude, construction_year), ~ mean(.x, na.rm = TRUE))
    ) %>%
  rename(longitude_center = longitude, latitude_center = latitude)
sites_with_spatial_data <- left_join(
  sites_with_spatial_data,
  locations %>% select(-construction_year),
  by = "location_construction_year"
)


### b Distance to river -------------------------------------------------------

### Prepare data ###
coordinates_plots <- sites %>%
  st_coordinates()
coordinates_danube_isar <- danube_isar %>%
  st_coordinates() %>%
  as_tibble() %>%
  select(X, Y)
### Calculate distances ###
distance <- geosphere::dist2Line(
  p = coordinates_plots, line = coordinates_danube_isar
  ) %>%
  as_tibble() %>%
  mutate(distance = round(distance, digits = 0)) %>%
  rename(river_distance = distance)
river_distance <- distance$river_distance
sites_with_spatial_data <- sites_with_spatial_data %>%
  add_column(river_distance)
### Plot for proof ###
dist_sf <- st_as_sf(distance, coords = c("lon", "lat")) %>%
  st_set_crs(value = 4326)
ggplot() +
  geom_sf(data = danube_isar, fill = "grey50", color = "grey50") +
  geom_sf(data = sites) +
  geom_sf(data = dist_sf, colour = "grey60")


### c Amount of surrounding habitats ------------------------------------------

### Prepare data ###
buffered_plots <- sites %>%
  st_buffer(dist = 500)
biotope_types <- biotope_mapping %>%
  select(id, datum, titel, haupt_typ, neben_typ, geometry) %>%
  filter(
    str_detect(haupt_typ, "Flachland-Mähwiesen") |
      str_detect(haupt_typ, "Magerrasen") |
      str_detect(haupt_typ, "Artenreiches Extensivgrünland") |
      str_detect(haupt_typ, "Magere Altgrasbestände") |
      str_detect(haupt_typ, "Pfeifengraswiese") |
      str_detect(neben_typ, "Flachland-Mähwiesen") |
      str_detect(neben_typ, "Magerrasen") |
      str_detect(neben_typ, "Artenreiches Extensivgrünland") |
      str_detect(neben_typ, "Pfeifengraswiese")
      )
sites_with_biotopes <- buffered_plots %>%
  st_intersection(biotope_types) %>%
  st_make_valid() %>%
  group_by(id) %>%
  summarise(st_union(geometry))
### Calculate area sizes ###
biotope_area <- sites_with_biotopes %>%
  mutate(biotope_area = st_area(sites_with_biotopes),
         biotope_area = round(biotope_area, digits = 0)) %>%
  st_drop_geometry()
sites_with_spatial_data <- sites_with_spatial_data %>%
  left_join(biotope_area, by = "id")
### Plot for proof ###
ggplot() +
  geom_sf(data = danube_isar, fill = "grey50", color = "grey50") +
  geom_sf(data = buffered_plots, fill = "transparent") +
  geom_sf(data = sites, size = .01) +
  geom_sf(data = sites_with_biotopes, fill = "red", color = "red")


### d Distance to closest habitat ---------------------------------------------

coordinates_biotopes <- biotope_types %>%
 as(Class = "Spatial")
distance <- geosphere::dist2Line(
  p = coordinates_plots, line = coordinates_biotopes
  ) %>%
  as_tibble() %>%
  rename(biotope_distance = distance)
biotope_distance <- distance$biotope_distance
sites_with_spatial_data <- sites_with_spatial_data %>%
  add_column(biotope_distance) %>%
  mutate(biotope_distance = round(biotope_distance, digits = 0))
### Plot for proof ###
dist_sf <- st_as_sf(distance, coords = c("lon", "lat")) %>%
  st_set_crs(value = 4326)
ggplot() +
  geom_sf(data = danube_isar, fill = "grey50", color = "grey50")
  geom_sf(data = buffered_plots, fill = "transparent") +
  geom_sf(data = sites, size = 0.01) +
  geom_sf(data = dist_sf, colour = "grey60")

  

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# C Save #######################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


### Save shp files ###
save(
  background_google,
  file = paste0(here("data", "processed", "spatial"),
                "/", "background_google.rda")
  )
save(
  background_toner,
  file = paste0(here("data", "processed", "spatial"),
                "/", "background_toner.rda")
  )
save(
  background_terrain,
  file = paste0(here("data", "processed", "spatial"),
                "/", "background_terrain.rda")
  )
st_write(
  germany,
  layer = "germany_epsg4326.shp", driver = "ESRI Shapefile",
  delete_layer = TRUE, dsn = here("data", "processed", "spatial")
  )
st_write(
  rivers,
  layer = "rivers_epsg4326.shp", driver = "ESRI Shapefile",
  delete_layer = TRUE, dsn = here("data", "processed", "spatial")
  )
st_write(
  grazing,
  layer = "grazing_epsg4326.shp", driver = "ESRI Shapefile",
  delete_layer = TRUE, dsn = here("data", "processed", "spatial")
  )
st_write(
  dikes,
  layer = "dikes_epsg4326.shp", driver = "ESRI Shapefile",
  delete_layer = TRUE, dsn = here("data", "processed", "spatial")
  )
st_write(
  conservation_area,
  layer = "conservation_area_epsg4326.shp", driver = "ESRI Shapefile",
  delete_layer = TRUE, dsn = here("data", "processed", "spatial")
  )
st_write(
  ffh_area,
  layer = "ffh_area_epsg4326.shp", driver = "ESRI Shapefile",
  delete_layer = TRUE, dsn = here("data", "processed", "spatial")
  )
st_write(
  sites,
  layer = "sites_epsg4326.shp", driver = "ESRI Shapefile",
  delete_layer = TRUE, dsn = here("data", "processed", "spatial")
  )
### River layer (was only once digitized) ###
#st_write(danube_isar, layer = "danube_isar_digitized_epsg4326.shp", driver = "ESRI Shapefile", delete_layer = TRUE, dsn = here("data", "processed", "spatial"))

st_write(
  obj = rivers,
  layer = "rivers_epsg4326.shp", driver = "ESRI Shapefile",
  delete_layer = FALSE, dsn = here("data", "processed", "spatial")
)
### Save data frames ###
write_csv(
  sites_with_spatial_data,
  file = here("data", "processed", "spatial", "sites_with_spatial_data.csv")
  )
write_csv(
  locations,
  file = here("data", "processed", "spatial", "locations.csv")
  )

### Save gpx file ###
sites_with_spatial_data %>%
  select(id, longitude, latitude) %>%
  mutate(id = as.character(id)) %>%
  as.data.frame() %>%
  pgirmess::writeGPX(type = "w", filename = here(
    "data", "processed", "spatial",
    "danube_old_dikes_plots.gpx"
  ))
