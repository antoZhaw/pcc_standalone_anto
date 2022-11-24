# libraries ---------------------------------------------------------------

# install packages

library(lidR) # summary tables
library(papeR) # summary tables
library(tidyverse) # tidy essentials (ggplot, purr, tidyr, readr, dplyr)
library(lubridate) # handling dates and time
# library(tmap) # map visualization
# library(leaflet) # interactive maps
library(terra) # handling spatial data
library(sf) # handling spatial data
library(janitor) # clean and consistent naming
library(forcats) # handling factor levels
library(raster) # rasterizing vector data
library(knitr) # for pretty tables
library(RCSF) # for CSF ground classification
library(geometry) # for raserize_canopy function
library(lmom) # for Key structural features of boreal forests

# Functions------
is.veg <- function(G) {
  is_veg <- G > 120 
  return(list(veg = is_veg))
}

# Read LAS file-----------------------------------------------------------------
# Intensity (i), color information (RGB), number of Returns (r), classification (c)
# of the first point is loaded only to reduce computational time.
las = readLAS(r"(C:\Daten\math_gubelyve\tls_data\saane_20211013_subsample_onlyRGBpts.las)", select = "xyziRGBrc", filter = "-keep_first")

las <- add_attribute(las, 0, "RGBmean")
las$RGBmean <- (las$R + las$G + las$B)/3
# las$veg <- if_else(las$G >= 23051, T, F)

# plot(las, size = 3, color = "RGB", bg = "white")

summary(las$Classification)
typeof(las$Classification)

# Point Metrics calculations (untested, heavy duty)------
# Add attribute on point level
# M <- point_metrics(las, ~is.planar(X,Y,Z), k = 20, filter = ~Classification != LASGROUND)
# Run metrics computation
#  un(M <- point_metrics(las, ~is.veg(G), k = 20)
# merge the output with the point cloud to visualize the result

# Add manual LAS attribute within R---------------------------------------------
# las <- add_attribute(las, FALSE, "veg")
# las$veg <- if_else(las$G >= 23051, T, F)
# plot(las, color = "veg")

# Ground filter with CSF--------------------------------------------------------
# setup csf filter settings
# rigidness: does not seem to have much impact.
# class_threshold and cloth_resolution influence each other. 0.5 x 0.5 is more conservative compared to 0.5 x 1.
mycsf <- csf(TRUE, 0.5, 1, rigidness = 3)
# apply ground classification
las <- classify_ground(las, mycsf)
summary(las$Classification)

# filter non-ground part from classified las------------------------------------
nongnd <- filter_poi(las, Classification %in% c(LASNONCLASSIFIED, LASUNCLASSIFIED))
# plot(nongnd, size = 3, color = "RGB", bg = "white")

# alternative approach (untested)
# poi <- ~Classification == 2
# can <- classify_poi(las, LASHIGHVEGETATION, poi = poi)

# filter ground from classified las---------------------------------------------
gnd <- filter_ground(las)
# plot(gnd, size = 3, color = "RGB", bg = "white")

# plot(can, size = 3, color = "RGB", bg = "white")

# Filter functions--------------------------------------------------------------
# las_sub = filter_poi(las, Classification %in% c(LASGROUND, LASWATER))

# Exploration of Colors---------------------------------------------------------

par(mfrow = c(2,1))
hist(nongnd$RGBmean)
hist(gnd$RGBmean)

hist(nongnd$R)
hist(gnd$R)

hist(nongnd$G)
hist(gnd$G)

hist(nongnd$B)
hist(gnd$B)


# Classification with L-moments-------------------------------------------------
# cloud_metrics(las, func = ~as.list(lmom::samlmu(Z)))



# Maybe worth trying out
# ?classify_poi()


print(las)
summary(las)


# Negation of attributes is also possible (all except intensity and angle)
# las = readLAS(LASfile, select = "* -i -a")

# Prefer filter() befor filter_poi() since it does not read it on C++ level
# show all available filters
readLAS(filter = "-help")
