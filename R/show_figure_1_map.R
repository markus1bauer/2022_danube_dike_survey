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
setwd(here("data", "processed", "spatial"))

### Load data ###
germany <- st_read("germany_epsg4326.shp")
danube <- st_read(here("data", "raw", "spatial",
                       "danube_isar_digitized_epsg4326.shp"))
danube$river[2] <- "Isar"
filter <- read_csv(here("data", "processed",
                        "data_processed_sites_spatial.csv"),
  col_names = TRUE, na = c("na", "NA"), col_types =
    cols(
      .default = "?"
    )
) %>%
  filter(surveyYear == 2017)
sites <- st_read("sites_epsg4326.shp") %>%
  semi_join(filter, by = "plot") %>%
  mutate(plot = factor(plot))
sites_ggmap <- st_coordinates(sites) %>%
  as_tibble()
grazing <- st_read("grazing_epsg4326.shp")
conservation_area <- st_read("conservation_area_epsg4326.shp")
ffh_area <- st_read("ffh_area_epsg4326.shp")
dikes <- st_read("dikes_epsg4326.shp")
locations <- read_csv("locations.csv",
  col_names = TRUE, col_types =
    cols(
      locationYear = "f"
    )
) %>%
  semi_join(filter, by = "locationYear")
load("background_toner.rda")
load("background_terrain.rda")
load("background_google.rda")



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Plot ################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


themeMB <- function() {
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


## 1 Map with ggmap #####################################################

### a Map of project site -----------------------------------------------

(graphSites <- ggmap(background_terrain,
  base_layer = ggplot(sites_ggmap, aes(x = X, y = Y))
) +
  geom_point(size = 2, color = "black", pch = 15) +
  geom_text_repel(
    data = locations, aes(
      label = constructionYear,
      x = longitude_center,
      y = latitude_center
    ),
    min.segment.length = 0
  ) +
  scale_x_continuous(breaks = seq(10, 15, 0.1)) +
  scale_y_continuous(breaks = seq(48, 50, 0.1)) +
  coord_sf(crs = st_crs(4326)) +
  # scale_fill_brewer(palette = "Greens", type = "seq", direction = 1, na.value = "grey", name = "Class of P") +
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
  themeMB() +
  theme(
    legend.position = c(0.8, 0.8),
    legend.background = element_rect(linetype = "solid", colour = "black")
  )
)

### b Germany -----------------------------------------------------------

graphGermany <- ggplot() +
  geom_sf(data = germany, fill = "transparent", colour = "black") +
  geom_point(aes(x = 12.885, y = 48.839), size = 1) +
  themeMB() +
  theme(
    plot.background = element_blank()
  )

### c Inset -------------------------------------------------------------

graphSites + inset_element(graphGermany,
  left = .7,
  bottom = .65,
  right = .99,
  top = .99,
  on_top = TRUE
)


### d Save --------------------------------------------------------------

ggsave("figure_1_map_ggmap_300dpi_17x11cm.tiff",
  dpi = 300, width = 17, height = 11, units = "cm",
  path = here("outputs", "figures")
)


## 2 Map with ggplot2 ###################################################

### a Map of project site -----------------------------------------------
set.seed(2)
(graphSites <- ggplot() +
  # geom_sf(data = ffh_area, fill = "grey50", color = "grey50") +
  geom_sf(data = dikes, colour = "grey60") +
  geom_sf(data = danube, colour = "black", size = 1) +
  geom_label_repel(
    data = locations, aes(
      label = locationYear,
      x = longitude_center,
      y = latitude_center
    ),
    min.segment.length = 0, box.padding = .6, fill = "white"
  ) +
  geom_sf(data = sites, colour = "red", size = 2) +
  coord_sf(crs = st_crs(4326)) +
  annotate("text", x = 13.11, y = 48.72, angle = -54, label = "Danube") +
  annotate("text", x = 12.79, y = 48.72, angle = 32, label = "Isar") +
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
  themeMB() +
  theme(
    legend.position = c(0.8, 0.8),
    legend.background = element_rect(linetype = "solid", colour = "black")
  ))

### b Inset -------------------------------------------------------------

set.seed(2)
graphSites + inset_element(graphGermany,
  left = .72,
  bottom = .65,
  right = .99,
  top = .99,
  on_top = TRUE
)

### c Save --------------------------------------------------------------

ggsave("figure_1_map_ggplot_300dpi_17x11cm.tiff",
  dpi = 300, width = 17, height = 11, units = "cm",
  path = here("outputs", "figures")
)


## 3 Map with tmap ######################################################

### a Map of project site -----------------------------------------------

tmap_mode("plot")
# tm_shape(ffh_area) +
# tm_fill(col = "grey40") +
# tm_shape(conservation_area) +
# tm_fill(col = "grey60") +
tmap <- tm_shape(danube) +
  tm_lines(col = "grey40") +
  tm_text("river", ymod = 1.2) +
  tm_shape(dikes) +
  tm_lines() +
  tm_shape(sites) +
  tm_dots(col = "red", size = .2, shape = 16) +
  tm_compass(position = c("left", "bottom"), size = 2) +
  tm_scale_bar(position = c("left", "bottom", with = 0.4)) +
  tm_layout(frame = FALSE)
tmap_ger <- tm_shape(germany) +
  tm_borders(col = "black") +
  tm_layout(frame = FALSE)

### b Save --------------------------------------------------------------

tmap_save(tmap,
  insets_tm = tmap_ger,
  insets_vp = viewport(
    x = unit(3.1, "cm"),
    y = unit(4.3, "cm"),
    width = unit(3, "cm"),
    height = unit(4, "cm")
  ),
  filename = paste0(here("outputs", "figures"), "/",
                    "figure_1_map_tmap_300dpi_8x11cm.tiff"),
  dpi = 300
)
