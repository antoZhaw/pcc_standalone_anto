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
# good thresholds for white noise filter between 40000...(45000)...48000
poi_whitenoise <- ~if_else(las$RGBmean >= 45000, T, F)

poi_veg_GLI <- ~if_else(las$GLI >= 0.10 & 
                          las$Classification == LASNONCLASSIFIED, T, F)
# GLI filters a broad range from greyish and yellowish parts.
# GLI = 0.15 is conservative, no sediment and cliff is affected
# GPI = 0.11 is ideal.
# GPI = 0.09 boarder of watercourse is affected entirely.
# GPI = 0.07 is already the lower limit, cliff and sediment points are affected.


poi_veg_GPI <- ~if_else(las$GPI >= 0.37 & 
                          las$Classification == LASNONCLASSIFIED, T, F)
# GPI filters a broad range from greyish and yellowish parts.
# GPI = 0.4 is conservative, no sediment and cliff is affected
# GPI = 0.37 is ideal.
# GPI = 0.36 boarder of watercourse is affected entirely.
# GPI = 0.35 is already the lower limit, cliff and sediment points are affected.

poi_veg_ExG <- ~if_else(las$ExG >= 7000 & 
                          las$Classification == LASNONCLASSIFIED, T, F)
# ExG filters vegetation in general
# ExG = 7500, border of watercourse is affected partly.
# ExG = 6000, border of watercourse is affected almost entirely
# ExG = 5500, border of watercourse and some sediment is affected.

poi_veg_ExGR <- ~if_else(las$ExGR >= 14000 & 
                           las$Classification == LASNONCLASSIFIED, T, F)
# ExGR filters specially bright green and blue points, yellowish points not
# ExGR = 14000, Point of cliff are affected
# ExGR = 10000, Parts of cliff are affected

poi_sky_ExB <- ~if_else(las$ExB >= 13500 & 
                          las$Classification == LASNONCLASSIFIED, T, F)
# ExB = 18500, already takes away some cliff parts.
# ExB = 14500, doubles cliff part but no sediment.
# ExB = 13500, sediment is affected.

poi_sky_BPI <- ~if_else(las$BPI >= 0.391 & 
                          las$Classification == LASNONCLASSIFIED, T, F)
#BPI = 0.391 represents 3rd Qu. and already takes away sediment


poi_rock_wide <- ~if_else(las$RtoB >= 0.8221 & las$RtoB <= 1.00 &
                           las$GtoB >= 0.9428 & las$GtoB <= 1.0507 &
                           las$ground == F &
                           las$Classification == LASNONCLASSIFIED, T, F)

# Maximales RtoB: 1.21 in cliff_bright
# Minimales RtoB: 0.722 in cliff_dark kann aber zu sky zugeordnet werden.

poi_rock_nar <- ~if_else(las$RBtimesGB >= 0.9370 & las$RBtimesGB <= 1.0433 &
                           las$ground == F &
                           las$Classification == LASNONCLASSIFIED, T, F)

poi_sed_wide <- ~if_else(las$RtoB >= 0.8221 & las$RtoB <= 1.00 &
                           las$GtoB >= 0.9428 & las$GtoB <= 1.0507 &
                           las$ground == T &
                           las$Classification == LASNONCLASSIFIED, T, F)

# Maximales RtoB: 1.21 in cliff_bright
# Minimales RtoB: 0.722 in cliff_dark kann aber zu sky zugeordnet werden.

poi_sed_nar <- ~if_else(las$RBtimesGB >= 0.9370 & las$RBtimesGB <= 1.0433 &
                           las$ground == T &
                           las$Classification == LASNONCLASSIFIED, T, F)


# Unused formulas---------------------------------------------------------------

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

# Read LAS file-----------------------------------------------------------------
# Intensity (i), color information (RGB), number of Returns (r), classification (c)
# of the first point is loaded only to reduce computational time.
las_origin = readLAS(r"(C:\Daten\math_gubelyve\tls_data\saane_20211013_subsample_onlyRGBpts.las)", select = "xyzRGBc", filter = "-keep_first")

# Create copy of read LAS to omit loading procedure.
las <- las_origin

# Display loaded classes (supposed to be zero)
# factor(las$Classification)

# Data exploration--------------------------------------------------------------
# plot(las, size = 1, color = "Intensity", bg = "black")


# Create attributes for classification------------------------------------------
# Add RGBmean attribute
las <- add_attribute(las, 0, "RGBmean")
las$RGBmean <- (las$R + las$G + las$B)/3

# Add Ratio R to B
las <- add_attribute(las, 0, "RtoB")
las$RtoB <- (las$R/las$B)

# Add Ratio G to B
las <- add_attribute(las, 0, "GtoB")
las$GtoB <- (las$G/las$B)

# Add Factor RtoB times GtoB
las <- add_attribute(las, 0, "RBtimesGB")
las$RBtimesGB <- (las$RtoB*las$GtoB)

# Add Green Percentage Index GPI
# las <- add_attribute(las, 0, "GPI")
# las$GPI <- (las$G/(las$R+las$G+las$B))

# Add Excess Green Index ExG
# las <- add_attribute(las, 0, "ExG")
# las$ExG <- (2*las$G-las$R-las$B)

# Add Blue Percentage Index BPI
# las <- add_attribute(las, 0, "BPI")
# las$BPI <- (las$B/(las$R+las$G+las$B))

# Add Excess Blue Index ExB
las <- add_attribute(las, 0, "ExB")
las$ExB <- (2*las$B-las$R-las$G)

# Add Excess Green minus Excess Red Index ExGR
# las <- add_attribute(las, 0, "ExGR")
# las$ExGR <- ((2*las$G-las$R-las$B)-(2*las$R-las$G-las$B))

# Add Green Leaf Index GLI
las <- add_attribute(las, 0, "GLI")
las$GLI <- (((las$G-las$R)+(las$G-las$B))/(2*las$G+las$R+las$B))

# Buggy forumulas---------------------------------------------------------------

# Add RGB Vegetation Index RGBVI
# Produces NAs only
# las <- add_attribute(las, 0, "RGBVI")
# las$RGBVI <- ((2*las$G-(las$B*las$R))/(2*las$G+(las$B*las$R)))

# Add Visible Atmospherically Resistant Index VARI
# Ranges from -Inf to Inf
# las <- add_attribute(las, 0, "VARI")
# las$VARI <- ((las$G-las$R)/(las$G+las$R-las$B))

# Add Normalised Green Red Difference Index VARI
# Throws an error
# las <- add_attribute(las, 0, "NGRDI")
# las$NGDRI <- ((las$G-las$R)/(las$G+las$R))

# Segment Ground with Cloth Simulation Filter-----------------------------------
# setup csf filter settings
# rigidness: does not seem to have much impact.
# class_threshold and cloth_resolution influence each other. 0.5 x 0.5 is more conservative compared to 0.5 x 1.
mycsf <- csf(TRUE, 0.5, 1, rigidness = 3)
# apply ground classification
las <- classify_ground(las, mycsf)
# gnd <- filter_ground(las)
las <- add_attribute(las, FALSE, "ground")
las$ground <- if_else(las$Classification == LASGROUND, T, F)

# filter non-ground part from classified las
# nongnd <- filter_poi(las, Classification %in% c(LASNONCLASSIFIED, LASUNCLASSIFIED))
# plot(nongnd, size = 3, color = "RGB", bg = "white")

# Reset class LASGROUND for further procedure
las$Classification <- LASNONCLASSIFIED

# Check whether classes are reset (Levels are supposed to be zero)
# factor(las$Classification)

# Classify noise----------------------------------------------------------------

# Classify white noise 
las <- classify_poi(las, class = LASNOISE, poi = poi_whitenoise)
las <- filter_poi(las, Classification != LASNOISE)
# las_noise <- filter_poi(las, Classification == LASNOISE)
# plot(las_noise, size = 1, color = "RGB", bg = "black")

# Classify sky------------------------------------------------------------------
# Vegetation filter priority: ExB, BPI (some might be deactivated)

las <- classify_poi(las, class = LASWIREGUARD, poi = poi_sky_ExB)
# las <- filter_poi(las, Classification != LASWIREGUARD)
# las_sky_nar <- filter_poi(las, Classification == LASWIREGUARD)
# plot(las_sky_nar, size = 1, color = "RGB", bg = "black")


# Testing-----------------------------------------------------------------------
# summary(las$GLI)

# las <- las_origin
# las <- classify_poi(las, class = LASHIGHVEGETATION, poi = poi_veg_GLI)
# las <- filter_poi(las, Classification != LASHIGHVEGETATION)
# las_test <- filter_poi(las, Classification == LASHIGHVEGETATION)
# plot(las_test, size = 1, color = "RGB", bg = "black")

# plot(las, size = 1, color = "RGB", bg = "black")
# plot(las, size = 1, color = "ExG", bg = "black")
# plot(las, size = 1, color = "GPI", bg = "black")
# plot(las, size = 1, color = "GLI", bg = "black")
# plot(las, size = 1, color = "ExGR", bg = "black")
 
# Explore them
# summary(las$RGBmean)
# summary(las$GtoB)
# hist(las$RGBmean)
# max(las$RGBmean)


# Classify vegetation-----------------------------------------------------------
# Vegetation filter priority: GLI, GPI, ExG, ExGR (some might be deactivated)

las <- classify_poi(las, class = LASLOWVEGETATION, poi = poi_veg_GLI)
# las <- filter_poi(las, Classification != LASLOWVEGETATION)
# las_veg_nar <- filter_poi(las, Classification == LASLOWVEGETATION)
# plot(las_veg_nar, size = 1, color = "RGB", bg = "black")

# Classify sediment-------------------------------------------------------------

# wide approach
# las <- classify_poi(las, class = LASKEYPOINT, poi = poi_sed_wide)
# las_sed_wide <- filter_poi(las, Classification == LASKEYPOINT)
# plot(las_sed_wide, size = 1, color = "RGB", bg = "black")

# narrow approach
las <- classify_poi(las, class = LASLOWPOINT, poi = poi_sed_nar)
# las_sed_nar <- filter_poi(las, Classification == LASLOWPOINT)
# plot(las_sed_nar, size = 1, color = "RGB", bg = "black")

# Classify rocks and cliffs-----------------------------------------------------

# wide approach
# las <- classify_poi(las, class = LASWIRECONDUCTOR, poi = poi_rock_wide)
# las_rock_wide <- filter_poi(las, Classification == LASWIRECONDUCTOR)
# plot(las_rock_wide, size = 1, color = "RGB", bg = "black")

# narrow approach
las <- classify_poi(las, class = LASROADSURFACE, poi = poi_rock_nar)
# las_rock_nar <- filter_poi(las, Classification == LASROADSURFACE)
# plot(las_rock_nar, size = 1, color = "RGB", bg = "black")


# Plot classified point cloud---------------------------------------------------

# Filter out noise and unclassified points for reduced computational time.
las <- filter_poi(las, Classification != LASNOISE)
las <- filter_poi(las, Classification != LASNONCLASSIFIED)
las <- filter_poi(las, Classification != LASUNCLASSIFIED)
# Sediment currently disabled
las <- filter_poi(las, Classification != LASROADSURFACE)
plot(las, size = 1, color = "Classification", bg = "black", legend = T)


# Outdated stuff----------------------------------------------------------------
# Point Metrics calculations (untested, heavy duty)
# Add attribute on point level
# M <- point_metrics(las, ~is.planar(X,Y,Z), k = 20, filter = ~Classification != LASGROUND)
# Run metrics computation
#  un(M <- point_metrics(las, ~is.veg(G), k = 20)
# merge the output with the point cloud to visualize the result

# Alternative approach for Ground filtration (untested)
# poi <- ~Classification == 2
# can <- classify_poi(las, LASHIGHVEGETATION, poi = poi)

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
