library(terra) # handling spatial data
library(sf) # handling spatial data
library(tmap) # map visualization
library(raster) # rasterizing vector data

t0_tif <- terra::rast("C:/Daten/math_gubelyve/tiff_data/20202008_Sarine_RGB_ppk_GCP02_transparent_mosaic_group1.tif")
# look at the metadata
t0_tif
# Plot the three layered geotiff
plot(t0_tif)
# Plot the combination
plotRGB(t0_tif)
# Task 5: Adding a background map
# in tmap tm_shape() calls data, and a tm_* object defines how the data is visualized
tmap_mode("plot")


tm2 <- tm_shape(t0_tif) +
  tm_raster(palette = "viridis") +
  tm_layout(legend.outside = TRUE)
tm2

tm1 <- tm_shape(t0_tif) + 
  tm_rgb(r=1, g=2, b=3)
tm1

oranges <- tmaptools::get_brewer_pal("Oranges", n = 2, contrast = c(0.3, 0.9))
tm_nw_100 <-
  tmap_mode("plot") + # "plot" or "view"
  tm_shape(tm_sed_t0) +
  tm_raster(palette = oranges, title = "Nests", alpha = 1) +
  # tm_shape(wallows_raster_100) +
  # tm_raster(palette = purples, title = "Wallows", alpha = 1) +
  tm_view(control.position = c("right", "top"))

tm_nw_100