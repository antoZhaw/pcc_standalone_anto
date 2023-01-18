# libraries --------------------------------------------------------------------

# install packages

library(lidR) # Point cloud classification
library(tidyverse)
library(janitor) # clean and consistent naming
library(lubridate)
library(terra) # handling spatial data
library(sf) # handling spatial data

# globals-----------------------------------------------------------------------
user <- Sys.getenv("USERNAME")

dir_repo <- if_else(user == "gubelyve", "C:/Daten/math_gubelyve/pcc_standalone", "C:/code_wc/pcc_standalone")
dir_data <- "C:/Daten/math_gubelyve"
path_input_data <- file.path(dir_data, "tls_data/2021/wholeset/pointcloud_2021_lv95_reg.las")
path_output_data <- file.path(dir_data, "tls_data/2021/wholeset/2021-tls-wholeset-2.las")

# read data---------------------------------------------------------------------
# Using filter argument empty discards essential information such as RGB.
las <- readLAS(path_input_data)


# Classify outliers as noise
# It computes first the average distance of each point to its k nearest neighbors.
# A higher k leads to a more precise but also higher computational time.
# Then it rejects the points that are farther than the average distance
# plus a number of times the standard deviation (multiplier m).
k <- 150
# m = 5 is quiet aggressive and takes away tips of leafes.
m <- 5
las <- classify_noise(las, sor(k, m))
las_outliers <- filter_poi(las, Classification %in% c(LASNOISE))
plot(las_outliers, size = 1, color = "RGB", bg = "white")
las <- filter_poi(las, Classification != LASNOISE)

writeLAS(las, file = path_output_data)
