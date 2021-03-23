# Show interactive map of the Danube old dikes ####
# Markus Bauer



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Preparation ################################################################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


### Packages ###
library(here)
library(tidyverse)
library(sf)
library(ggmap)
library(leaflet)
library(leaflet.extras)
library(htmltools)
library(mapview)

### Start ###
rm(list = ls())
setwd(here("data/processed/spatial"))

### Load data ###
wms_flood <- "https://www.lfu.bayern.de/gdi/wms/wasser/ueberschwemmungsgebiete?"
sites <- st_read("sites.shp")
ffh_area <- st_read("ffh_area.shp")
dikes <- st_read("dikes.shp")

(map <- sites %>%
  leaflet() %>% 
  addTiles(group = "OSM") %>%
  addProviderTiles("Esri.WorldTopoMap", group = "Esri") %>%
  setView(lng = 12.885, lat = 48.839, zoom = 10) %>% 
  addPolygons(data = ffh_area, color = "grey", opacity = .5, fillColor = "grey", fillOpacity = .6, weight = 0,
              popup = ~paste0("<b>", htmlEscape(NAME), "</b>", "<br/>",
                              htmlEscape(ID)),
              group = "FFH Areas") %>%
  addPolylines(data = dikes, color = "black", opacity = 1, weight = 1,
               group = "Dikes") %>%
  addCircleMarkers(radius = 3, color = "red", opacity = 1,
                   popup = ~paste0("<b>",htmlEscape(locatin), "</b>", "<br/>",
                                   htmlEscape(id)),
                   label = ~paste0(locatin, " (ID: ", id, ")"),
                   group = "Plots") %>%
  addWMSTiles(wms_flood,
              layers = "Deich",
              group = "HQ100",
              options = WMSTileOptions(format = "image/png", transparent = F)) %>%
  addLayersControl(baseGroups = c("OSM", "Esri"),
                   overlayGroups = c("FFH Areas", "Dikes", "Plots", "HQ100"),
                   options = layersControlOptions(collapsed = F)) %>%
  addScaleBar() %>%
  addMeasure(primaryLengthUnit = "kilometers",
             secondaryLengthUnit = F,
             activeColor = "red",
             completedColor = "darkred",
             primaryAreaUnit = "sqkilometers",
             localization = "en") %>%
  addMiniMap(
    toggleDisplay = T,
    tiles = providers$OSM,
    zoomLevelOffset = -6,
    minimized = T
  ))

setwd(here("outputs/figures"))
mapshot(map, url = "figure_1_map_interactive.html")
