# Show map of the Danube dike experiment ####
# Markus Bauer
# Citation: Markus Bauer, Jakob Huber, Katharina Beck, Johannes Kollmann (xxxx)


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Preparation ################################################################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

### Packages ###
library(tidyverse)
library(sf)
#library(sp)
library(patchwork)

### Start ###
rm(list = ls())
setwd("Z:/Documents/0_Donaudeiche/3_Aufnahmen_und_Ergebnisse/2022_Danube_dike_experiment/data/processed/shp_files")

### Load data ###
ger <- st_read("germany.shp")
map <- load("map_experiment.rda")



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Plot ##############################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
themeMB <- function(){
  theme(
    panel.background = element_rect(fill = NA),
    text  = element_text(size = 10, color = "black"),
    axis.line.y = element_line(),
    axis.line.x = element_line(),
    axis.title = element_blank(),
    legend.key = element_rect(fill = "white"),
    legend.position = "none",
    legend.margin = margin(0, 0, 0, 0, "cm"),
    plot.margin = margin(.5, 0, 0, 0, "cm")
  )
}


## 1 Preparation ##############################################################################

### a Garchinger Heide -----------------------------------------------------------------------
(expGraph <- ggmap(experiment) +#, 
                   #base_layer = ggplot(sites, aes(lon, lat))) +
    #geom_point(data = sites, aes(x = lon, y = lat, color = age))
    #geom_sf(data = sites, aes(shape = dataset), size = .7) +
    #scale_shape_manual(values = c(1, 3)) +
    scale_x_continuous(breaks = seq(0, 12.9, 0.005)) +
    scale_y_continuous(breaks = seq(48, 50, 0.005)) +
    #ggspatial::annotation_scale(width_hint = 0.4, height = unit(0.2, "cm")) +
    ggspatial::annotation_north_arrow(which_north = "true", style = ggspatial::north_arrow_fancy_orienteering(), height = unit(1, "cm"), width = unit(1, "cm"), pad_y = unit(0.6, "cm")) +
    themeMB() +
    theme(axis.text.y = element_text(angle = 90, hjust = .5)))
### b Germany -----------------------------------------------------------------------
library(tmap)
tm_shape(ger) +
   tm_borders() + 
   tm_compass()
(gerGraph <- ggplot() +
   geom_sf(data = ger, fill = "transparent") +
   geom_point(aes(x = 12.885, y = 48.839), size = 1) +
   themeMB() +
   theme(
     plot.background = element_blank(),
     axis.ticks = element_blank(), 
     axis.text = element_blank(),
     axis.line.x = element_blank(),
     axis.line.y = element_blank()
   ))
### c Inset -----------------------------------------------------------------------
expGraph + inset_element(gerGraph, .0, .6, .5, .9, on_top = T)


# 2 Save ##############################################################################

ggsave("figure_1_map_(300dpi_8x11cm).tiff", 
       dpi = 300, width = 8, height = 11, units = "cm",
       path = "Z:/Documents/0_Donaudeiche/3_Aufnahmen_und_Ergebnisse/2022_Danube_dike_experiment/outputs/figures")
