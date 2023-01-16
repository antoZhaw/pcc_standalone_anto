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
path_input_data <- file.path(dir_data, "tls_data/2021/wholeset/pointcloud_2021_lv95_reg.las")
path_output_data <- file.path(dir_data, "tls_data/2021/wholeset/2021-tls-wholeset-2.las")


# read data---------------------------------------------------------------------
# Using filter argument empty discards essential information such as RGB.
las <- readLAS(path_input_data)
summary(las)

if (has.lasRGB(las)) 
{stop("The read LAS file has no colour information, script stops.")}
