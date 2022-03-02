# Beta diversity on dike grasslands
# Figure A5 ####
# Markus Bauer
# 2022-01-11



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Preparation ##########################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


### Packages ###
library(here)
library(tidyverse)
library(sf)
library(ggmap)
library(leaflet)
library(leaflet.extras)
library(leaflet.extras2)
library(htmltools)
library(mapview)

### Start ###
rm(list = ls())
setwd(here("data", "processed", "spatial"))

### Load data ###
wms_flood <- "https://www.lfu.bayern.de/gdi/wms/wasser/ueberschwemmungsgebiete?"
sites <- st_read("sites_epsg4326.shp")
ffh_area <- st_read("ffh_area_epsg4326.shp")
dikes <- st_read("dikes_epsg4326.shp")

(map <- sites %>%
  leaflet() %>%
  addTiles(group = "OSM") %>%
  addProviderTiles("Esri.WorldTopoMap", group = "Esri") %>%
  setView(lng = 12.885, lat = 48.839, zoom = 10) %>%
  addPolygons(
    data = ffh_area,
    color = "grey60", opacity = 0, weight = 0,
    fillColor = "grey60", fillOpacity = .6,
    highlight = highlightOptions(
      weight = 2, color = "blue",
      opacity = .8, bringToFront = FALSE
    ),
    popup = ~ paste0(
      "<b>", htmlEscape(NAME), "</b>", "<br/>",
      "ID: ", htmlEscape(ID)
    ),
    group = "FFH Areas"
  ) %>%
  addWMSTiles(
    baseUrl = wms_flood,
    layers = "hwgf_hq100",
    group = "HQ100",
    options = WMSTileOptions(format = "image/png", transparent = TRUE)
  ) %>%
  addPolylines(
    data = dikes, color = "black", opacity = 1, weight = 2,
    label = ~ paste0(
      "Baujahr: ", htmlEscape(BAUJAHR), ", Sanierung:",
      htmlEscape(SANIERUNG)
    ),
    popup = ~ paste0(
      "Baujahr: ", htmlEscape(BAUJAHR), "<br/>",
      "Sanierung: ", htmlEscape(SANIERUNG)
    ),
    highlight = highlightOptions(weight = 2, color = "blue", opacity = .8),
    group = "Dikes"
  ) %>%
  addCircleMarkers(
    radius = 5, color = "red", weight = 2, opacity = 1, fillOpacity = .5,
    popup = ~ paste0(
      "<b>", htmlEscape(locatin), "</b>", "<br/>",
      "ID: ", htmlEscape(id)
    ),
    label = ~ paste0(locatin, " (ID: ", id, ")"),
    group = "Plots"
  ) %>%
  addLayersControl(
    baseGroups = c("OSM", "Esri"),
    overlayGroups = c("FFH Areas", "Dikes", "Plots", "HQ100"),
    options = layersControlOptions(collapsed = FALSE)
  ) %>%
  addScaleBar() %>%
  addMeasure(
    primaryLengthUnit = "kilometers",
    secondaryLengthUnit = FALSE,
    activeColor = "red",
    completedColor = "darkred",
    primaryAreaUnit = "sqkilometers",
    localization = "en"
  ) %>%
  addMiniMap(
    toggleDisplay = TRUE,
    tiles = providers$OSM,
    zoomLevelOffset = -6,
    minimized = TRUE
  ))

mapshot(map, url = here("outputs", "figures", "figure_a7_map_interactive.html"))


leaflet() %>%
  addTiles() %>%
  setView(lng = 12.885, lat = 48.839, zoom = 11) %>%
  addWMSTiles(
    baseUrl = wms_flood,
    layers = "hwgf_hq100",
    group = "HQ100",
    options = WMSTileOptions(format = "image/png", transparent = FALSE)
  )
