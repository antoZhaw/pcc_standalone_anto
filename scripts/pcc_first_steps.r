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
library(purrr)


# Functions---------------------------------------------------------------------
is.veg <- function(G) {
  is_veg <- G > 120 
  return(list(veg = is_veg))
}

is.LAScorrupt <- function(las) {
  has_R <- if_else(mean(las$R) != 0,T,F)
  has_G <- if_else(mean(las$G) != 0,T,F)
  has_B <- if_else(mean(las$B) != 0,T,F)
  valid <- has_R & has_G & has_B
  corrupt <- !valid
  return(corrupt)
}

to.LAScolor <- function(small_RGB) {
  # According to LAS Speciﬁcation 1.4 - R14 of 
  # the American Society for Photogrammetry and Remote Sensing (ASPRS)
  # A normalization of each pixel channel to a two byte integer is recommended.
  LAScolor <- small_RGB * 256 
  return(as.integer(LAScolor))
}

gen.attribute.plot <- function(input_attr, attr_name, plot_title, sub_title, post) {
  # receive attribute name, uncleaned "$" might cause errors.
  suffix <- if_else(post == T, "_post", "_pre")
  filename <- paste("export/", plot_title, "_", attr_name, suffix, ".png", sep = "")
  ggplot(las@data) +
    aes(x = input_attr, fill = as.factor(las$Classification)) + 
    geom_density(alpha = 0.5) + 
    labs(title = plot_title, 
         subtitle = sub_title,
         x = attr_name) +
    theme_minimal() +
    theme(legend.position = c(.9, .90),
          legend.title = element_blank())
  ggsave(filename = filename, bg = "white", units = "mm")
  # alternative print function, which is faster but throws an error.
  # dev.print(file=filename, device=png, width=800)
}

# Globals-----------------------------------------------------------------------
start <- lubridate::now()
date <- as.Date(start)
hour <- hour(start)
minute <- minute(start)

timestamp <- as.character(paste(date, hour, minute, sep = "-"))

user <- Sys.getenv("USERNAME")

wholeset <- F
year <- "2021"
perspective <- "tls"
settype <- if_else(wholeset == T, "wholeset", "subset")

if(user == "gubelyve"){
  dir_repo <- "C:/Daten/math_gubelyve/pcc_standalone"
  dir_data <- "C:/Daten/math_gubelyve/tls_data"
  
} else{
  dir_repo <- "C:/code_wc/"
  dir_data <- "C:/Daten/math_gubelyve/tls_data"
}

if(wholeset){
  dataset <- "wholeset_221011.las"
} else{
  dataset <- "saane_20211013_subsample_onlyRGBpts.las"
}

output_id <- as.character(paste(timestamp, perspective, settype, year, sep = "-"))
output_las_name <- as.character(paste(output_id, ".las", sep = ""))
output_las_path <- file.path(dir_data, year, settype, output_las_name, fsep="/")

data_path <- file.path(dir_data, year, settype, dataset)

output_las_path

sky_upper_RGB <- as.integer(c("150", "175", "250")) %>% to.LAScolor()
# 1st Quantile for Blue = 170, 120 takes away more sediment.
sky_lower_RGB <- as.integer(c("25", "50", "120")) %>% to.LAScolor()

# currently unused
darksky_upper_RGB <- as.integer(c("67", "94", "108")) %>% to.LAScolor()
darksky_lower_RGB <- as.integer(c("25", "43", "54")) %>% to.LAScolor()

# Formulas----------------------------------------------------------------------
poi_whitenoise <- ~if_else(las$RGBmean >= whitenoise_tresh, T, F)

poi_blacknoise <- ~if_else(las$RGBmean <= blacknoise_tresh, T, F)

poi_sky_ExB <- ~if_else(las$ExB >= ExB_tresh & las$ground == F &
                          las$Classification == LASNONCLASSIFIED, T, F)

poi_sky_ExB_band <- ~if_else(las$ExB >= ExB_tresh_min & las$ExB <= ExB_tresh_max &
                          las$Classification == LASNONCLASSIFIED, T, F)

poi_sky_BPI <- ~if_else(las$BPI >= BPI_tresh & las$ground == F &
                          las$Classification == LASNONCLASSIFIED, T, F)

poi_red_ExR <- ~if_else(las$ExR >= ExR_tresh &
                          las$Classification == LASNONCLASSIFIED, T, F)

poi_sed_ExR_band <- ~if_else(las$ExR >= ExR_tresh_min & las$ExR <= ExR_tresh_max &
                          las$ground == T &
                          las$Classification == LASNONCLASSIFIED, T, F)

poi_red_RPI <- ~if_else(las$RPI >= RPI_tresh & 
                          las$Classification == LASNONCLASSIFIED, T, F)

poi_veg_GLI <- ~if_else(las$GLI >= GLI_tresh & 
                          las$Classification == LASNONCLASSIFIED, T, F)

poi_veg_GPI <- ~if_else(las$GPI >= GPI_tresh & 
                          las$Classification == LASNONCLASSIFIED, T, F)

poi_veg_ExG <- ~if_else(las$ExG >= ExG_tresh & 
                          las$Classification == LASNONCLASSIFIED, T, F)

poi_veg_ExGR <- ~if_else(las$ExGR >= ExGR_tresh & 
                           las$Classification == LASNONCLASSIFIED, T, F)

poi_rock_ratios <- ~if_else(las$RtoB >= RtoBmin & las$RtoB <= RtoBmax &
                              las$GtoB >= GtoBmin & las$GtoB <= GtoBmax &
                              las$ground == F &
                              las$Classification == LASNONCLASSIFIED, T, F)

poi_sed_ratios <- ~if_else(las$RtoB >= RtoBmin & las$RtoB <= RtoBmax &
                             las$GtoB >= GtoBmin & las$GtoB <= GtoBmax &
                             las$ground == T &
                             las$Classification == LASNONCLASSIFIED, T, F)

poi_rock_times <- ~if_else(las$RBtimesGB >= RBtimesGB_min & las$RBtimesGB <= RBtimesGB_max &
                           las$ground == F &
                           las$Classification == LASNONCLASSIFIED, T, F)

# Maximales RtoB: 1.21 in cliff_bright
# Minimales RtoB: 0.722 in cliff_dark kann aber zu sky zugeordnet werden.

poi_sed_times <- ~if_else(las$RBtimesGB >= RBtimesGB_min & las$RBtimesGB <= RBtimesGB_max &
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
las <- readLAS(data_path, select = "xyzRGBci", filter = "-keep_first")
# -keep_xy 4343860 542298 4344110 542457
if (is.LAScorrupt(las)) {stop("The read LAS file has no colour information, script stops.")}

# Create copy of read LAS to omit loading procedure.
# las <- las_origin

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
las <- add_attribute(las, 0, "GPI")
las$GPI <- (las$G/(las$R+las$G+las$B))

# Add Excess Green Index ExG
las <- add_attribute(las, 0, "ExG")
las$ExG <- (2*las$G-las$R-las$B)

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

# Add Red Percentage Index RPI
las <- add_attribute(las, 0, "RPI")
las$RPI <- (las$R/(las$R+las$G+las$B))

# Add Excess Red Index ExR
las <- add_attribute(las, 0, "ExR")
las$ExR <- (2*las$R-las$G-las$B)

# Buggy Indices-----------------------------------------------------------------

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


# Generate attribute plots before classification--------------------------------
# tbd - Starting point for a for loop, piping might be better

active_attr <- names(las)
active_attr <- active_attr[! active_attr %in% c("X","Y","Z","Classification", "RtoB", "RGtoB", "RBtimesGB", "Intensity" )]
 
static_subtitle <- "(0) No class, (3) Veg., (6) Sky, (8) Sediment, (10) Rocks."
las_post <- F

map(active_attr, function(x){
  gen.attribute.plot(las[[x]], x, output_id, static_subtitle, las_post)
})

# Segment Ground with Cloth Simulation Filter-----------------------------------
# setup csf filter settings
# rigidness: does not seem to have much impact.
# class_threshold and cloth_resolution influence each other. 0.5 x 0.5 is more conservative compared to 0.5 x 1.

class_tresh <- if_else(wholeset == T, 0.5, 0.5)
cloth_res <- if_else(wholeset == T, 0.5, 1)
rigid <- if_else(wholeset == T, 1, 3)

mycsf <- csf(TRUE, class_threshold = class_tresh, cloth_res, rigidness = rigid)
# apply ground classification
las <- classify_ground(las, mycsf)
# gnd <- filter_ground(las)
las <- add_attribute(las, FALSE, "ground")
las$ground <- if_else(las$Classification == LASGROUND, T, F)

las_gnd <- filter_poi(las, Classification == LASGROUND)
plot(las_gnd, size = 1, color = "RGB", bg = "white")


# filter non-ground part from classified las
# nongnd <- filter_poi(las, Classification %in% c(LASNONCLASSIFIED, LASUNCLASSIFIED))
# plot(nongnd, size = 3, color = "RGB", bg = "white")

# Reset class LASGROUND for further procedure
las$Classification <- LASNONCLASSIFIED

# Check whether classes are reset (Levels are supposed to be zero)
# factor(las$Classification)

# Classify noise----------------------------------------------------------------

# Classify white noise 
whitenoise_tresh <- 38000
# good thresholds for white noise filter between 40000...(45000)...48000

las <- classify_poi(las, class = LASNOISE, poi = poi_whitenoise)
# las <- filter_poi(las, Classification != LASNOISE)

# Classify black noise
blacknoise_tresh <- 6000
# good thresholds for black noise filter between 4000...(6000)...8000

las <- classify_poi(las, class = LASNOISE, poi = poi_blacknoise)

# Plot filtered noise
# las_noise <- filter_poi(las, Classification == LASNOISE)
# plot(las_noise, size = 1, color = "RGB", bg = "white")

las <- filter_poi(las, Classification != LASNOISE)

# Classify sky------------------------------------------------------------------
# Vegetation filter priority: ExB, BPI (some might be deactivated)

ExB_tresh <- 3500
# ExB = 18500, already takes away some cliff parts.
# ExB = 14500, doubles cliff part but no sediment.
# ExB = 13500, cliff wall is affected.
# ExB = 3500, cliff wall is affected on a wide range but still works.

BPI_tresh <- 0.391
#BPI = 0.391 represents 3rd Qu. and already takes away sediment

# las_origin <- las
# las <- las_origin

las <- classify_poi(las, class = LASBUILDING, poi = poi_sky_ExB)
# las <- filter_poi(las, Classification != LASBUILDING)
# las_sky <- filter_poi(las, Classification == LASBUILDING)
# plot(las_sky, size = 1, color = "RGB", bg = "black")

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
summary(las$ExR)
# hist(las$RGBmean)
# max(las$RGBmean)

# Classify green parts of vegetation--------------------------------------------
# Vegetation filter priority: GLI, ExG or GPI, ExGR (some might be deactivated)

GLI_tresh <- 0.09
# GLI filters a broad range from greyish and yellowish parts.
# GLI = 0.15 is conservative, no sediment and cliff is affected
# GPI = 0.11 is ideal.
# GPI = 0.09 boarder of watercourse is affected entirely.
# GPI = 0.07 cliff and sediment points are affected.
# GPI = 0.04 bush and canopy around sediment areas are included, except brown.

GPI_tresh <- 0.35
# GPI filters a broad range from greyish and yellowish parts.
# GPI = 0.4 is conservative, no sediment and cliff is affected
# GPI = 0.37 is ideal.
# GPI = 0.36 boarder of watercourse is affected entirely.
# GPI = 0.35 is already the lower limit, cliff and sediment points are affected.

ExG_tresh <- 4500
# ExG filters vegetation in general, neglects rather dark points
# ExG = 7500, border of watercourse is affected partly.
# ExG = 6000, border of watercourse is affected almost entirely
# ExG = 5500, border of watercourse and some sediment is affected.
# ExG = 4500, lower limit.

ExGR_tresh <- 14000
# ExGR filters specially bright green and blue points, yellowish points not
# ExGR = 14000, Point of cliff are affected
# ExGR = 10000, Parts of cliff are affected

las <- classify_poi(las, class = LASLOWVEGETATION, poi = poi_veg_GLI)
# las <- filter_poi(las, Classification != LASLOWVEGETATION)

# Classify Red parts of vegetation----------------------------------------------
# Red filter priority: ExR, RPI (some might be deactivated)
ExR_tresh <- 20000
# ExR = 30000, broad vegetation is affected, some sediment parts too.
# ExR = 20000, first lines of cliff relief is affected.

RPI_tresh <- 0.9

las <- classify_poi(las, class = LASLOWVEGETATION, poi = poi_red_ExR)
# las <- filter_poi(las, Classification != LASBRIGDE)
# las_red <- filter_poi(las, Classification == LASBRIGDE)
# plot(las_red, size = 1, color = "RGB", bg = "black")

# Plot vegetation
# las_veg <- filter_poi(las, Classification == LASLOWVEGETATION)
# plot(las_veg, size = 1, color = "RGB", bg = "black")

# Classify sediment-------------------------------------------------------------

# Maximales RtoB: 1.21 in cliff_bright
# Minimales RtoB: 0.722 in cliff_dark kann aber zu sky zugeordnet werden.

#Min: 0.7222 and Max. 1.1290 derived from cliff_dark
RtoBmin <- 0.7222
RtoBmax <- 1.1290

#Min: 0.8333 and Max. 1.1613 derived from cliff_dark
#Min: 0.7928 from cliff_blue
GtoBmin <- 0.7928
GtoBmax <- 1.1613

# Approach with ratios RtoB and GtoB
las <- classify_poi(las, class = LASKEYPOINT, poi = poi_sed_ratios)

# Classify band of negative excess red parts------------------------------------
# Negative ExR values appear to be sediment, ground criteria is set true.
ExR_tresh_max <- 1000
ExR_tresh_min <- -53000
# ExR = -53000 ... 1000 seems to affect sediment points more than others
# No big difference between -20000 and -53000

las <- classify_poi(las, class = LASKEYPOINT, poi = poi_sed_ExR_band)

# For an exclusive plot of negative ExR values change to LASBRIDGE as class.
# las_sed_ratios <- filter_poi(las, Classification == LASBRIGDE)
# plot(las_sed_ratios, size = 1, color = "RGB", bg = "white")

# Approach with RtoB times GtoB
# Set limits for sediment.
RBtimesGB_min <- 0.9370
RBtimesGB_max <- 1.0433
# las <- classify_poi(las, class = LASLOWPOINT, poi = poi_sed_times)
# las_sed_times <- filter_poi(las, Classification == LASLOWPOINT)
# plot(las_sed_nar, size = 1, color = "RGB", bg = "black")


# Classify rocks and cliffs-----------------------------------------------------

# Set narrow thresholds to filter only very grey points. Only active in poi_rock_ratios
# RtoBmin <- 0.96
# RtoBmax <- 1.04
# GtoBmin <- 0.96
# GtoBmax <- 1.04

# Approach with ratios RtoB and GtoB
las <- classify_poi(las, class = LASRAIL, poi = poi_rock_ratios)
# las_rock_ratios <- filter_poi(las, Classification == LASRAIL)
# plot(las_rock_ratios, size = 1, color = "RGB", bg = "black")


# Approach with RtoB times GtoB
# Set limits again for rock. If not set, limits of sediment filter is applied.
RBtimesGB_min <- 0.98
RBtimesGB_max <- 1.02
# las <- classify_poi(las, class = LASRAIL, poi = poi_rock_times)
# las_rock_times <- filter_poi(las, Classification == LASRAIL)
# plot(las_rock_times, size = 1, color = "RGB", bg = "black")

# Test - Classify band of negative excess blue parts----------------------------
# Negative ExB values appear to be sediment, ground criteria is set true.
# ExB_tresh_max <- 3499
# ExB_tresh_min <- -8000
# ExB = -8000 ... 3499 are greyish points but no distinct classes
# ExB = -60000 ... -8000 are greenish, greyish points but no distinct classes
# No big difference between -20000 and -53000

# las <- classify_poi(las, class = LASBRIGDE, poi = poi_sky_ExB_band)
# las_exb_neg <- filter_poi(las, Classification == LASBRIGDE)
# plot(las_exb_neg, size = 1, color = "RGB", bg = "white")
# summary(las$ExB)


# Generate attribute plots after classification--------------------------------
las_post = T

map(active_attr, function(x){
  gen.attribute.plot(las[[x]], x, output_id, static_subtitle, las_post)
})

# Save generated output---------------------------------------------------------
writeLAS(las, file = output_las_path)

# Plot classified point cloud---------------------------------------------------

# Show the unclassified---------------------------------------------------------
las_foreveralone <- filter_poi(las, Classification == LASNONCLASSIFIED)
plot(las_foreveralone, size = 1, color = "RGB", bg = "black")

# Filter out noise and unclassified points for reduced computational time.
las <- filter_poi(las, Classification != LASNOISE)
las <- filter_poi(las, Classification != LASNONCLASSIFIED)
las <- filter_poi(las, Classification != LASUNCLASSIFIED)

# disable rocks since it does not work properly
# las <- filter_poi(las, Classification != LASRAIL)
# plot(las, size = 1, color = "Classification", bg = "black")

# Plot separated classes
# las_sky <- filter_poi(las, Classification == LASBUILDING)
# plot(las_sky, size = 1, color = "RGB", bg = "black")

las_veg <- filter_poi(las, Classification == LASLOWVEGETATION)
plot(las_veg, size = 1, color = "RGB", bg = "black")

las_sed_ratios <- filter_poi(las, Classification == LASKEYPOINT)
plot(las_sed_ratios, size = 1, color = "RGB", bg = "white")

end <- lubridate::now()
end

diff <- end - start
diff

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

# Alternative filter functions--------------------------------------------------------------
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


# Help lines
# Negation of attributes is also possible (all except intensity and angle)
# las = readLAS(LASfile, select = "* -i -a")

# Prefer filter() befor filter_poi() since it does not read it on C++ level
# show all available filters
readLAS(filter = "-help")
