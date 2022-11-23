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


# Read LAS file
# Read settings: Intensity (i), color information (RGB) of the first point is loaded only to reduce computational time.
las = readLAS(r"(C:\Daten\math_gubelyve\tls_data\saane_20211013_subsample_onlyRGBpts.las)", select = "xyziRGB", filter = "-keep_first")
# las_check(las)

# setup csf filter settings
# rigidness: does not seem to have much impact.
# class_threshold and cloth_resolution influence each other. 0.5 x 0.5 is more conservative compared to 0.5 x 1.
mycsf <- csf(TRUE, 0.5, 1, rigidness = 3)
# apply ground classification
las <- classify_ground(las, mycsf)
# filter ground from classified las
gnd <- filter_ground(las)
plot(gnd, size = 3, color = "RGB", bg = "white")




plot(las, color ="RGB")


print(las)
summary(las)


# Negation of attributes is also possible (all except intensity and angle)
# las = readLAS(LASfile, select = "* -i -a")

# Prefer filter() befor filter_poi() since it does not read it on C++ level
# show all available filters
readLAS(filter = "-help")
