library(lidR) # Point cloud classification
library(papeR) # summary tables
library(tidyverse) # tidy essentials (ggplot, purr, tidyr, readr, dplyr)
library(lubridate) # handling dates and time
library(tmap) # map visualization
# library(leaflet) # interactive maps
library(terra) # handling spatial data
library(sf) # handling spatial data
library(rgdal) # export raster objects
library(ncdf4) # export ncdf files
library(janitor) # clean and consistent naming
library(forcats) # handling factor levels
library(raster) # rasterizing vector data
library(knitr) # for pretty tables
library(rayshader) # for pretty DTMs
library(RCSF) # for CSF ground classification
library(RMCC) # for MCC ground classification
library(geometry) # for raserize_canopy function
library(lmom) # for Key structural features of boreal forests
library(purrr) # for map function
library(rjson) # for JSON generation
library(rgl) # for RGL Viewer control functions
library(sabre) # for mapcurves (goodness-of-fit function)
library(fasterize) # for faster rasterization

# Generate rectangular polygon of area of interest
gen_xy_tls <- structure(list(dat = c("AOI TLS", "AOI TLS", 
                                     "AOI TLS", "AOI TLS"),
                             Longitude = c(2575340, 2575340, 
                                           2575480, 2575480),
                             Latitude = c(1178520, 1178780, 
                                          1178780, 1178520)),
                        class = "data.frame", row.names = c(NA,-4L))

bounding_box_tls <- gen_xy_tls %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 2056) %>%
  group_by(dat) %>%
  summarise(geometry = st_combine(geometry)) %>%
  st_cast("POLYGON")


# Read Shapefiles
csf_aoi_shp <- read_sf(dsn = "C:/Daten/math_gubelyve/pcc_standalone/data/dut_filter_230315/dut_filter_230315.shp") 

ggplot() +
  geom_sf(data = bounding_box_tls, aes(fill = dat), alpha = 0.1) +
  geom_sf(data = csf_aoi_shp, aes(fill = as.factor(id)), alpha = 0.5) +
  coord_sf(datum=st_crs(2056)) +
  theme_minimal()

t0_tif <- terra::rast("C:/Daten/math_gubelyve/tiff_data/20202008_Sarine_RGB_ppk_GCP02_transparent_mosaic_group1.tif")
# look at the metadata
t0_tif
# Plot the three layered geotiff
plot(t0_tif)
# Plot the combination
plotRGB(t0_tif)
# Task 5: Adding a background map
# in tmap tm_shape() calls data, and a tm_* object defines how the data is visualized
# tmap_mode("plot")

tm1 <- tm_shape(csf_aoi_shp) +
  tm_borders(col="red", lwd=1) +
  tm_shape(bounding_box_tls) +
  tm_borders(col="black", lwd=1)
  # tm_shape(t0_tif) + 
  # tm_rgb(r=1, g=2, b=3)
  # tm_basemap()
tm1
