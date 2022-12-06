# libraries ---------------------------------------------------------------

# Fragen an Nils
# 1) Eher attribute erstellen und einfache Formulas oder umgekehrt?
# 2) Operatoren && und &
# 3) Bug in poi_darksky

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

# Functions---------------------------------------------------------------------
is.veg <- function(G) {
  is_veg <- G > 120 
  return(list(veg = is_veg))
}

to.LAScolor <- function(small_RGB) {
  # According to LAS Speciﬁcation 1.4 - R14 of 
  # the American Society for Photogrammetry and Remote Sensing (ASPRS)
  # A normalization of each pixel channel to a two byte integer is recommended.
  LAScolor <- small_RGB * 256 
  return(as.integer(LAScolor))
}

# Globals-----------------------------------------------------------------------

sky_upper_RGB <- as.integer(c("150", "175", "250")) %>% to.LAScolor()
# 1st Quantile for Blue = 170, 120 takes away more sediment.
sky_lower_RGB <- as.integer(c("25", "50", "120")) %>% to.LAScolor()

# currently unused
darksky_upper_RGB <- as.integer(c("67", "94", "108")) %>% to.LAScolor()
darksky_lower_RGB <- as.integer(c("25", "43", "54")) %>% to.LAScolor()



# Formulas----------------------------------------------------------------------

poi_whitenoise <- ~if_else(las$RGBmean >= 45000, T, F)

# poi_sky <- ~if_else(las$R <= sky_upper_RGB[1] & las$R >= sky_lower_RGB[1] &
#                       las$G <= sky_upper_RGB[2] & las$G >= sky_lower_RGB[2] &
#                       las$B <= sky_upper_RGB[3] & las$B >= sky_lower_RGB[3]
#                     , T, F)

# Factor does not work for some reason.
# poi_darksky <- ~if_else(las$R <= darksky_upper_RGB[1] & las$R >= darksky_lower_RGB[1] &
                    #   las$G <= darksky_upper_RGB[2] & las$G >= darksky_lower_RGB[2] &
                    #   las$B <= darksky_upper_RGB[3] & las$B >= darksky_lower_RGB[3] &
                    #     as.numeric(las$B/las$R) >= 1.40
                    # , T, F)

# RtoB überschneidet sich mit rock.
poi_sky <- ~if_else(las$RtoB <= 0.722, T, F)

# Maximales RtoB: 1.21 in cliff_bright
# Minimales RtoB: 0.722 in cliff_dark kann aber zu sky zugeordnet werden.
poi_rock <- ~if_else(las$RtoB >= 0.75 & las$RtoB <= 1.22, T, F)

# Minimales RtoB: 1.265, minimales GtoB: 1.374 (Stand 6.12.22)
poi_veg <- ~if_else(las$RtoB >= 1.25 & las$GtoB >= 1.3, T, F)

# Read LAS file-----------------------------------------------------------------
# Intensity (i), color information (RGB), number of Returns (r), classification (c)
# of the first point is loaded only to reduce computational time.
las_origin = readLAS(r"(C:\Daten\math_gubelyve\tls_data\saane_20211013_subsample_onlyRGBpts.las)", select = "xyzRGBci", filter = "-keep_first")
las <- las_origin

# Data exploration--------------------------------------------------------------
# plot(las, size = 1, color = "Intensity", bg = "black")

# Classify noise----------------------------------------------------------------
# Add RGBmean attribute
las <- add_attribute(las, 0, "RGBmean")
las$RGBmean <- (las$R + las$G + las$B)/3

# Add Ratio R to B
las <- add_attribute(las, 0, "RtoB")
las$RtoB <- (las$R/las$B)

# Add Ratio G to B
las <- add_attribute(las, 0, "GtoB")
las$GtoB <- (las$G/las$B)

# summary(las$RGBmean)
# hist(las$RGBmean)
# max(las$RGBmean)

# Classify white noise 
# good thresholds for white noise filter between 40000...(45000)...48000
las <- classify_poi(las, class = LASNOISE, poi = poi_whitenoise)
# las <- filter_poi(las, Classification != LASNOISE)
# las <- classify_noise(las, las$RGBmean <= 2000)
# plot(las, size = 1, color = "RGB", bg = "black")

# Classify blue points (classified as LASRAIL here)
las <- classify_poi(las, class = LASRAIL, poi = poi_sky)
# las_blue <- filter_poi(las, Classification == LASRAIL)
# plot(las_blue, size = 1, color = "RGB", bg = "black")

# Classify dark blue points (classified as LASWIREGUARD here), currently unused
# las <- classify_poi(las, class = LASWIREGUARD, poi = poi_darksky)
# las_darksky <- filter_poi(las, Classification == LASWIREGUARD)
# plot(las_darksky, size = 1, color = "RGB", bg = "white")

# Plot denoised LAS
las <- filter_poi(las, Classification != LASRAIL)
las <- filter_poi(las, Classification != LASNOISE)
# las <- filter_poi(las, Classification != LASWIREGUARD)
plot(las, size = 1, color = "RGB", bg = "black")

# Classify vegetation
las <- classify_poi(las, class = LASMEDIUMVEGETATION, poi = poi_veg)
# las_veg <- filter_poi(las, Classification == LASMEDIUMVEGETATION)
# plot(las_veg, size = 1, color = "RGB", bg = "black")

# Classify sediment
las <- classify_poi(las, class = LASROADSURFACE, poi = poi_rock)
# las_rock <- filter_poi(las, Classification == LASROADSURFACE)
# plot(las_rock, size = 1, color = "RGB", bg = "black")


# Point Metrics calculations (untested, heavy duty)-----------------------------
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
