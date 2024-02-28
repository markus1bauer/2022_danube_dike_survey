# Beta diversity on dike grasslands
# Plot map of study site ####

# Markus Bauer
# 2022-01-11



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Preparation #########################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



### Packages ###
library(here)
library(tidyverse)
library(sf)
library(ggmap)
library(tmap)
library(ggrepel)
library(RColorBrewer)
library(patchwork)
library(grid)

### Start ###
rm(list = ls())



## 1 Load data ################################################################


### a Load locations ----------------------------------------------------------

filter <- read_csv(
  here("data", "processed", "data_processed_sites_spatial.csv"),
  col_names = TRUE, na = c("na", "NA"), col_types =
    cols(
      .default = "?"
    )
) %>%
  filter(survey_year == 2017)
layer_sites <- st_read(
  here("data", "processed", "spatial", "sites_epsg4326.shp")
  ) %>%
  semi_join(filter, by = "plot") %>%
  mutate(plot = factor(plot))
sites_ggmap <- st_coordinates(layer_sites) %>%
  as_tibble()
locations <- read_csv(
  here("data", "processed", "spatial", "locations.csv"),
  col_names = TRUE, col_types =
    cols(
      location_construction_year = "f"
    )
) %>%
  semi_join(filter, by = "location_construction_year")
rm("filter")


### b Load land use and borders ------------------------------------------------

layer_danube <- st_read(
  here("data", "raw", "spatial", "danube_isar_digitized_epsg4326.shp")
  )
layer_danube$river[2] <- "Isar"
layer_dikes <- st_read(
  here("data", "processed", "spatial", "dikes_epsg4326.shp")
  )
layer_germany <- st_read(
  here("data", "processed", "spatial", "germany_epsg4326.shp")
  )
layer_rivers <- st_read(
  here("data", "processed", "spatial", "rivers_epsg4326.shp")
  )


### c Load background ----------------------------------------------------------

load(here("data", "processed", "spatial", "background_terrain.rda"))



## 2 Load theme ################################################################


theme_mb <- function() {
  theme(
    panel.background = element_rect(fill = NA),
    panel.grid = element_line(colour = NA),
    text = element_text(size = 10, color = "black"),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    legend.key = element_rect(fill = "white"),
    legend.margin = margin(0, 0, 0, 0, "cm"),
    plot.margin = margin(.5, 0, 0, 0, "cm")
  )
}




#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Plot ################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



## 1 Map with ggmap #####################################################


### a Map of project site -----------------------------------------------

(graph_sites <- ggmap(
  background_terrain,
  base_layer = ggplot(sites_ggmap, aes(x = X, y = Y))
  ) +
  geom_point(size = 2, color = "black", pch = 15) +
  geom_text_repel(
    data = locations,
    aes(label = construction_year, x = longitude_center, y = latitude_center),
    min.segment.length = 0
  ) +
  scale_x_continuous(breaks = seq(10, 15, 0.1)) +
  scale_y_continuous(breaks = seq(48, 50, 0.1)) +
  coord_sf(crs = st_crs(4326)) +
  ggspatial::annotation_scale(
    width_hint = 0.4,
    height = unit(0.2, "cm"),
    pad_y = unit(0.6, "cm"),
    pad_x = unit(0.7, "cm")
  ) +
  ggspatial::annotation_north_arrow(
    which_north = "true",
    style = ggspatial::north_arrow_fancy_orienteering(),
    height = unit(1, "cm"),
    width = unit(1, "cm"),
    pad_y = unit(0.9, "cm"),
    pad_x = unit(0.6, "cm")
  ) +
  theme_mb() +
  theme(
    legend.position = c(0.8, 0.8),
    legend.background = element_rect(linetype = "solid", colour = "black")
  )
)

### b Germany -----------------------------------------------------------

graph_germany <- ggplot() +
  geom_sf(data = layer_germany, fill = "transparent", colour = "black") +
  geom_point(aes(x = 12.885, y = 48.839), size = 1) +
  theme_mb() +
  theme(
    plot.background = element_blank()
  )


### c Inset -------------------------------------------------------------

graph_sites + inset_element(graph_germany,
  left = .7,
  bottom = .65,
  right = .99,
  top = .99,
  on_top = TRUE
)


### d Save --------------------------------------------------------------

# ggsave(
# "figure_1_map_ggmap_300dpi_17x11cm.tiff",
# dpi = 300, width = 17, height = 11, units = "cm",
# path = here("outputs", "figures")
# )



## 2 Map with ggplot2 ###################################################


### a Map of project site -----------------------------------------------

set.seed(2)
(graph_sites <- ggplot() +
  geom_sf(data = layer_rivers, color = "black", linewidth = .3) +
  geom_sf(data = layer_danube, colour = "black", linewidth = .8) +
  geom_sf(data = layer_dikes, colour = "grey60", linewidth = .6) +
  geom_label_repel(
    data = locations, aes(
      label = location_construction_year,
      x = longitude_center,
      y = latitude_center
    ),
    min.segment.length = 0, box.padding = .6, fill = "white"
  ) +
  geom_sf(data = layer_sites, colour = "red", size = 2) +
  coord_sf(crs = st_crs(4326)) +
  annotate("text", label = "Danube", x = 12.92, y = 48.92, angle = -30) +
  annotate("text", label = "Isar", x = 12.71, y = 48.695, angle = 22) +
  ggspatial::annotation_scale(
    width_hint = 0.4,
    height = unit(0.2, "cm"),
    pad_y = unit(0.6, "cm"),
    pad_x = unit(0.7, "cm")
  ) +
  ggspatial::annotation_north_arrow(
    which_north = "true",
    style = ggspatial::north_arrow_fancy_orienteering(),
    height = unit(1, "cm"),
    width = unit(1, "cm"),
    pad_y = unit(0.9, "cm"),
    pad_x = unit(0.6, "cm")
  ) +
  theme_mb() +
  theme(
    legend.position = c(0.8, 0.8),
    legend.background = element_rect(linetype = "solid", colour = "black")
  ))


### b Inset -------------------------------------------------------------

set.seed(2)
graph_sites + inset_element(graph_germany,
  left = .72,
  bottom = .65,
  right = .99,
  top = .99,
  on_top = TRUE
)


### c Save --------------------------------------------------------------

# ggsave(
#   "figure_1_map_ggplot_300dpi_17x11cm.tiff",
#   dpi = 300, width = 17, height = 11, units = "cm",
#   path = here("outputs", "figures")
# )



## 3 Map with tmap ######################################################


### a Map of project site -----------------------------------------------

tmap_mode("plot")
# tm_shape(ffh_area) +
# tm_fill(col = "grey40") +
# tm_shape(conservation_area) +
# tm_fill(col = "grey60") +
tmap <- tm_shape(layer_danube) +
  tm_lines(col = "grey40") +
  tm_text("river", ymod = 1.2) +
  tm_shape(layer_dikes) +
  tm_lines() +
  tm_shape(layer_sites) +
  tm_dots(col = "red", size = .2, shape = 16) +
  tm_compass(position = c("left", "bottom"), size = 2) +
  tm_scale_bar(position = c("left", "bottom", with = 0.4)) +
  tm_layout(frame = FALSE)
tmap_ger <- tm_shape(layer_germany) +
  tm_borders(col = "black") +
  tm_layout(frame = FALSE)


### b Save --------------------------------------------------------------

# tmap_save(tmap,
#   insets_tm = tmap_ger,
#   insets_vp = viewport(
#     x = unit(3.1, "cm"),
#     y = unit(4.3, "cm"),
#     width = unit(3, "cm"),
#     height = unit(4, "cm")
#   ),
#   filename = paste0(
#     here("outputs", "figures"), "/", "figure_1_map_tmap_300dpi_8x11cm.tiff"
#     ),
#   dpi = 300
# )
