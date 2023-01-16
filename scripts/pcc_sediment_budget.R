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
las_filter = "-drop_z_above 580
              -drop_x_below 2575330
              -drop_x_above 2575460
              -drop_y_below 1178450
              -drop_y_above 1178800"
# -drop_z_below 570

dir_repo <- if_else(user == "gubelyve", "C:/Daten/math_gubelyve/pcc_standalone", "C:/code_wc/pcc_standalone")
dir_data <- "C:/Daten/math_gubelyve"
path_las1 <- file.path(dir_data, "tls_data/2021/wholeset/output/2023-01-16-13-18-2021-tls-wholeset-2/2023-01-16-13-18-2021-tls-wholeset-2-sed.las")
path_las2 <- file.path(dir_data, "tls_data/2022/wholeset/output/2023-01-16-10-50-2022-tls-wholeset-2/2023-01-16-10-50-2022-tls-wholeset-2-sed.las")
path_output_data <- file.path(dir_data, "tls_data/2021/wholeset/2021-tls-wholeset-2.las")


# read data---------------------------------------------------------------------
# Using filter argument empty discards essential information such as RGB.
las1 <- readLAS(path_las1, filter = las_filter)
# min z 559.1937, max z 694.1542 

las2 <- readLAS(path_las2, filter = las_filter)
# min z 564.7118, max z 684.0745 

DEM1 <- rasterize_canopy(las1, res = 1, algorithm = p2r())
DEM2 <- rasterize_canopy(las2, res = 1, algorithm = p2r())

rs_DEM2 <- resample(DEM1, DEM2)
DEM3 <- rs_DEM2 - DEM2

col <- height.colors(30)
par(mfrow=c(1,3))

plot(DEM1, col = col, main = "TLS 2021")
plot(DEM2, col = col, main = "TLS 2022")
plot(DEM3, col = col, main = "DEM of difference")


# means that every triangle with an edge longer than 8 will be discarded from the triangulation.
# chm <- rasterize_canopy(las2, res = 0.5, algorithm = dsmtin(max_edge = 8))
col <- height.colors(25)
plot(chm, col = col)

readLAS(filter = "-help")