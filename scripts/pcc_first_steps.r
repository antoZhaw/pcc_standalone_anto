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

summary(las)

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

# filter ground from classified las
gnd <- filter_ground(las)
plot(gnd, size = 3, color = "RGB", bg = "white")

poi <- ~Classification == 2
can <- classify_poi(las, LASHIGHVEGETATION, poi = poi)

plot(can, size = 3, color = "RGB", bg = "white")


# Exploration of Colors---------------------------------------------------------

par(mfrow = c(2,1))
hist(can$R)
hist(gnd$R)

hist(can$G)
hist(gnd$G)

hist(can$B)
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
