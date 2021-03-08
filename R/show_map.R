# Show map of the Danube old dikes ####
# Markus Bauer
# Citation: Markus Bauer, Jakob Huber, Katharina Beck, Johannes Kollmann (xxxx)



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Preparation ################################################################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

### Packages ###
library(tidyverse)
library(sf)
library(ggmap)
library(ggrepel)
library(RColorBrewer)
library(patchwork)

### Start ###
rm(list = ls())
setwd("Z:/Documents/0_Donaudeiche/3_Aufnahmen_und_Ergebnisse/2022_Danube_old_dikes/data/processed/shp_files")
#remotes::install_github("poissonconsulting/poisspatial")
#ps_ggmap_to_raster(background_map)


### Load data ###
ger <- st_read("germany.shp")
sites <- st_read("sites.shp")
#sites <- do.call(rbind, st_geometry(sites)) %>% 
#  as_tibble() %>% 
#  setNames(c("lon","lat"))
sites2 <- read_csv2("sites2.csv", col_names = T, col_types = 
                      cols(
                        id = col_factor(),
                        phosphorousClass = col_factor(levels = c("A","B","C","D","E"))
                        )
                      )
blocks <- read_csv2("blocks.csv", col_names = T, col_types = 
                      cols(location = col_factor())
                    )
dikes_dataload("background_toner.rda")
load("background_terrain.rda")
load("background_google.rda")

plot(dikes)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Plot ##############################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
themeMB <- function(){
  theme(
    panel.background = element_rect(fill = NA),
    panel.grid = element_line(colour = NA),
    text  = element_text(size = 10, color = "black"),
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


## 1 Preparation ##############################################################################

### a Map of project site -----------------------------------------------------------------------
(sitesGraph <- ggmap(background_terrain, 
                      base_layer = ggplot(sites2, aes(x = lon_cent, y = lat_cent))) +
    geom_point(size = 2, color = "black", pch = 15) +
    geom_text_repel(data = blocks, aes(label = constructionYear, x = lon_cent, y = lat_cent),
                    min.segment.length = 0) +
    scale_x_continuous(breaks = seq(10, 15, 0.1)) +
    scale_y_continuous(breaks = seq(48, 50, 0.1)) +
    coord_sf(crs = st_crs(4326)) +
    #scale_fill_brewer(palette = "Greens", type = "seq", direction = 1, na.value = "grey", name = "Class of P") +
    ggspatial::annotation_scale(width_hint = 0.4, height = unit(0.2, "cm"), pad_y = unit(0.6, "cm"), pad_x = unit(0.7, "cm")) +
    ggspatial::annotation_north_arrow(which_north = "true", style = ggspatial::north_arrow_fancy_orienteering(), height = unit(1, "cm"), width = unit(1, "cm"), pad_y = unit(0.9, "cm"), pad_x = unit(0.6, "cm")) +
    themeMB() +
    theme(legend.position = c(0.8, 0.8),
          legend.background = element_rect(linetype = "solid", colour = "black"))
 )

### b Germany -----------------------------------------------------------------------
gerGraph <- ggplot() +
   geom_sf(data = ger, fill = "transparent", colour = "black") +
   geom_point(aes(x = 12.885, y = 48.839), size = 1) +
   themeMB() +
   theme(
     plot.background = element_blank()
   )

### c Inset -----------------------------------------------------------------------
sitesGraph + inset_element(gerGraph, .7, .65, .99, .99, on_top = T)



# 2 Save ##############################################################################

ggsave("figure_1_map_(300dpi_17x11cm).tiff", 
       dpi = 300, width = 17, height = 11, units = "cm",
       path = "Z:/Documents/0_Donaudeiche/3_Aufnahmen_und_Ergebnisse/2022_Danube_old_dikes/outputs/figures")
