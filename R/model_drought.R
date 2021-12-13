library(ncdf4) # package for netcdf manipulation
library(raster) # package for raster manipulation
library(rgdal) # package for geospatial analysis
library(here)
library(tidyverse)

#installr::install.Rtools(check = TRUE, check_r_update = TRUE, GUI = TRUE)
rm(list = ls())
setwd(here("data/raw"))


nc_data <- nc_open('data_drought.nc')
# Save the print(nc) dump to a text file
{
  sink('data_drought_metadata.txt')
  print(nc_data)
  sink()
}

lon <- ncvar_get(nc_data, "lon")
lat <- ncvar_get(nc_data, "lat", verbose = F)
t <- ncvar_get(nc_data, "time")

head(lon) # look at the first few entries in the longitude vector

ndvi.array <- ncvar_get(nc_data, "SMI") # store the data in a 3-dimensional array
dim(ndvi.array) 

fillvalue <- ncatt_get(nc_data, "SMI", "_FillValue")
fillvalue

nc_close(nc_data) 


ndvi.array[ndvi.array == fillvalue$value] <- NA

ndvi.slice <- ndvi.array[, , 1] 

dim(ndvi.slice)

r <- raster(t(ndvi.slice), xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat), 
            crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))

r <- flip(r, direction='y')

plot(r)



r_brick <- brick(ndvi.array, xmn=min(lat), xmx=max(lat), ymn=min(lon), ymx=max(lon), 
                 crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))

toolik_lon <- 12.86298
toolik_lat <- 48.83921
r_brick2 <- flip(t(r_brick), direction='y')
toolik_series <- extract(r_brick2, SpatialPoints(cbind(toolik_lon,toolik_lat)), method='simple')

toolik_df <- data.frame(year= seq(from=1, to=816, by=1), NDVI=t(toolik_series)) %>%
  filter(year > 792)

ggplot(data=toolik_df, aes(x=year, y=NDVI, group=1)) +
  geom_line() + # make this a line plot
  ggtitle("SMI") +     # Set title
  geom_hline(yintercept = c(0.05, 0.1, 0.02)) +
  theme_bw() # use the black and white theme

