# libraries --------------------------------------------------------------------

# install packages

library(lidR) # Point cloud classification
library(tidyverse)
library(janitor) # clean and consistent naming
library(lubridate)
library(terra) # handling spatial data
library(sf) # handling spatial data

# Functions---------------------------------------------------------------------
has.lasRGB <- function(las) {
  has_R <- if_else(mean(las$R) != 0,T,F)
  has_G <- if_else(mean(las$G) != 0,T,F)
  has_B <- if_else(mean(las$B) != 0,T,F)
  valid <- has_R & has_G & has_B
  noRGB <- !valid
  return(noRGB)
}

# globals-----------------------------------------------------------------------
user <- Sys.getenv("USERNAME")

dir_repo <- if_else(user == "gubelyve", "C:/Daten/math_gubelyve/pcc_standalone", "C:/code_wc/pcc_standalone")
dir_data <- "C:/Daten/math_gubelyve"
path_input_data <- file.path(dir_data, "tls_data/2022/wholeset/pointcloud_2022_lv95.las")
path_output_data <- file.path(dir_data, "tls_data/2022/wholeset/2022-tls-wholeset-2.las")


# read data---------------------------------------------------------------------
# Using filter argument empty discards essential information such as RGB.
las <- readLAS(path_input_data)
summary(las)

if (has.lasRGB(las)) 
{stop("The read LAS file has no colour information, script stops.")}

# Check whether a Coordinate Reference System is assigned already
st_crs(las)

crs_name <- st_crs(las)$Name
if(is.na(crs_name)){
  print("No CRS is assigned yet!")
}

# Assign Coordinate Reference System
# Go to https://spatialreference.org/ or epsg.io and find your CRS.
# Common EPSG: 
# CH1903+ / LV95 = 2056
# WGS84 = 

EPSG_nr <- 2056
las <- st_set_crs(las, EPSG_nr)
st_crs(las)$Name

writeLAS(las, file = path_output_data)

# Now one can transform the data
st_bbox(las)
tlas <- sf::st_transform(las, sf::st_crs(2056))
st_crs(tlas)$Name
st_bbox(tlas)

writeLAS(tlas, file = path_output_data)

