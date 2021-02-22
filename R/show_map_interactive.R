library(leaflet)

wms_administration <- "https://geoservices.bayern.de/wms/v1/ogc_verwaltungsatlas_verwaltungsgrenzen.cgi?"
wms_dGK25 <- "https://www.lfu.bayern.de/gdi/wms/geologie/dgk25?"
wms_uebk25 <- "https://www.lfu.bayern.de/gdi/wms/boden/uebk25?"
wms_flood <- "https://www.lfu.bayern.de/gdi/wms/wasser/ueberschwemmungsgebiete?"
wms_conservation <- "https://www.lfu.bayern.de/gdi/wms/natur/schutzgebiete?"
wms_biotope <- "https://www.lfu.bayern.de/gdi/wms/natur/biotopkartierung?"

leaflet() %>% 
  setView(lng = 12.885, lat = 48.839, zoom = 15) %>% 
  #addTiles(group = "OSM (default)") %>%
  #addProviderTiles(providers$Stamen.Toner, group = "Toner") %>%
  addProviderTiles(providers$Stamen.TonerLite, group = "Toner Lite") %>%
  addWMSTiles(
    wms_flood,
    layers = "Hochwassergefahrenflächen HQ100",
    options = WMSTileOptions(format = "image/png", transparent = TRUE)
  ) %>%
  addLayersControl(
    baseGroups = c("OSM (default)", "Toner", "Toner Lite"),
    #overlayGroups = c("Quakes", "Outline"),
    options = layersControlOptions(collapsed = FALSE)
  )

