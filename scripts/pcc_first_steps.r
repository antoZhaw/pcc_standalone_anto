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


# Read LAS file

print(las)
summary(las)


las = readLAS("data/TLS_Saane_2022_WGS84/ScanPos003 - SINGLESCANS - 221011_112204.las", select = "xyzi", filter = "-keep_first")
las_check(las)
plot(las)

# Negation of attributes is also possible (all except intensity and angle)
# las = readLAS(LASfile, select = "* -i -a")

# Prefer filter() befor filter_poi() since it does not read it on C++ level
# show all available filters
readLAS(filter = "-help")
