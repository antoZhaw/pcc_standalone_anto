# libraries ---------------------------------------------------------------
# install packages

library(lidR) # Point cloud classification
library(papeR) # summary tables
library(tidyverse) # tidy essentials (ggplot, purr, tidyr, readr, dplyr)
library(lubridate) # handling dates and time
# library(tmap) # map visualization
# library(leaflet) # interactive maps
library(terra) # handling spatial data
library(sf) # handling spatial data
library(rgdal) # export raster objects
library(ncdf4) # export ncdf files
library(janitor) # clean and consistent naming
library(forcats) # handling factor levels
library(raster) # rasterizing vector data
library(knitr) # for pretty tables
library(rayshader) # for pretty DTMs
library(RCSF) # for CSF ground classification
library(RMCC) # for MCC ground classification
library(geometry) # for raserize_canopy function
library(lmom) # for Key structural features of boreal forests
library(purrr) # for map function
library(rjson) # for JSON generation

# Functions---------------------------------------------------------------------

has.lasClassification <- function(las) {
  classified <- if_else(mean(las$Classification) != 0,T,F)
  return(classified)
}

is.lasCRScompliant <- function(las, target_epsg) {
  las_epsg <- st_crs(las)$epsg
  compliant <- if_else(las_epsg == target_epsg,T,F)
  return(compliant)
}

to.LAScolor <- function(small_RGB) {
  # According to LAS Speciﬁcation 1.4 - R14 of 
  # the American Society for Photogrammetry and Remote Sensing (ASPRS)
  # A normalization of each pixel channel to a two byte integer is recommended.
  LAScolor <- small_RGB * 256 
  return(as.integer(LAScolor))
}

gen.attribute.plot <- function(input_attr, attr_name, plot_title, sub_title, post, file_path) {
  # receive attribute name, uncleaned "$" might cause errors.
  suffix <- if_else(post == T, "_post", "_pre")
  filename <- paste(file_path, "/",  plot_title, "_", attr_name, suffix, ".png", sep = "")
  ggplot(las@data) +
    aes(x = input_attr, fill = fct_recode(as.factor(Classification), 
                                          "Nonclassified" = "0",
                                          "Unclassified" = "1",
                                          "Ground" = "2",
                                          "Low vegetation" = "3",
                                          "Mid vegetation" = "4",
                                          "High vegetation" = "5",
                                          "Building" = "6",
                                          "Low Point" = "7",
                                          "Sediment" = "8",
                                          "Water" = "9",
                                          "Rail" = "10",
                                          "Road" = "11",
                                          "Wireguard" = "12",
                                          "Wireconductor" = "13",
                                          "Tower" = "14",
                                          "Bridge" = "15",
                                          "Noise" = "16")) + 
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


classify_gnd <- function(las, class_thres, cloth_res, rigid) {
  mycsf <- csf(TRUE, class_thres, cloth_res, rigid)
  las <- classify_ground(las, mycsf)
  las_gnd <- filter_poi(las, Classification == LASGROUND)
  plot(las_gnd, size = 1, color = "RGB", bg = "white")
  return(T)
}

# Globals for Configuration-----------------------------------------------------
# Specify dataset
dataset_id <- "1"
wholeset <- T
year <- "2021"
perspective <- "uav"
settype <- if_else(wholeset == T, "wholeset", "subset")

# Internal globals such as paths and IDs----------------------------------------
# Record start date and time
start <- as_datetime(lubridate::now())
date <- as.Date(start)
hour <- hour(start)
minute <- minute(start)
# Generate timestamp for this run.
timestamp <- as.character(paste(date, hour, minute, sep = "-"))

# Build unique datasetname for reporting purposes.
datasetname <- as.character(paste(year, perspective, settype, dataset_id, sep = "-"))
dataset <- paste(datasetname, ".las", sep = "")

# Load environment dependent paths.
user <- Sys.getenv("USERNAME")
if(user == "gubelyve"){
  dir_repo <- "C:/Daten/math_gubelyve/pcc_standalone"
  dir_data <- "C:/Daten/math_gubelyve"
} else{
  dir_repo <- "C:/code_wc/"
  dir_data <- "C:/Daten/math_gubelyve"
}

dir_persp <- if_else(perspective == "tls", "tls_data", "uav_data")
dir_config <-  file.path(dir_repo, "config", fsep="/")

config_id <- as.character(paste(year, perspective, settype, dataset_id, sep = "-"))
output_id <- as.character(paste(timestamp, config_id, sep = "-"))
output_path <- file.path(dir_data, dir_persp, year, settype, "output", output_id, fsep="/")

config_json_name <- as.character(paste(config_id, ".json", sep = ""))
config_json_path <- file.path(dir_config, config_json_name, fsep="/")

aoi_path <- file.path(dir_repo, "data/area_of_interest_final", "AOI_final.shp", fsep="/")
aoi_dir <- file.path(dir_repo, "data/area_of_interest_final", fsep="/")

# load dataset specific parameter set
cfg <- fromJSON(file = config_json_path)

# Create run specific output folder
dir.create(output_path)

output_json_name <- as.character(paste(output_id, ".json", sep = ""))
output_json_path <- file.path(output_path, output_json_name, fsep="/")

output_las_inp_name <- as.character(paste(output_id, "-inp.txt", sep = ""))
output_las_inp_path <- file.path(output_path, output_las_inp_name, fsep="/")

output_las_aoi_name <- as.character(paste(output_id, "-aoi.txt", sep = ""))
output_las_aoi_path <- file.path(output_path, output_las_aoi_name, fsep="/")

output_las_class_name <- as.character(paste(output_id, "-classified.txt", sep = ""))
output_las_class_path <- file.path(output_path, output_las_class_name, fsep="/")

output_asc_name <- as.character(paste(output_id, ".asc", sep = ""))
output_asc_path <- file.path(output_path, output_asc_name, fsep="/")

output_ncdf_name <- as.character(paste(output_id, ".nc", sep = ""))
output_ncdf_path <- file.path(output_path, output_ncdf_name, fsep="/")

output_las_sed_name <- as.character(paste(output_id, "-sed.las", sep = ""))
output_las_sed_path <- file.path(output_path, output_las_sed_name, fsep="/")

output_las_gnd_name <- as.character(paste(output_id, "-gnd.las", sep = ""))
output_las_gnd_path <- file.path(output_path, output_las_gnd_name, fsep="/")

output_las_all_name <- as.character(paste(output_id, "-all.las", sep = ""))
output_las_all_path <- file.path(output_path, output_las_all_name, fsep="/")

data_path <- file.path(dir_data, dir_persp, year, settype, dataset)

# Formulas----------------------------------------------------------------------
poi_whitenoise <- ~if_else(las$RGBmean >= whitenoise_thresh, T, F)

poi_blacknoise <- ~if_else(las$RGBmean <= blacknoise_thresh, T, F)

poi_sky_ExB <- ~if_else(las$ExB >= ExB_thresh & las$ground == F &
                          las$Classification == LASNONCLASSIFIED, T, F)

poi_sky_ExB_band <- ~if_else(las$ExB >= ExB_thresh_min & las$ExB <= ExB_thresh_max &
                          las$Classification == LASNONCLASSIFIED, T, F)

poi_sky_BPI <- ~if_else(las$BPI >= BPI_thresh & las$ground == F &
                          las$Classification == LASNONCLASSIFIED, T, F)

poi_red_ExR <- ~if_else(las$ExR >= ExR_thresh &
                          las$Classification == LASNONCLASSIFIED, T, F)

poi_sed_ExR_band <- ~if_else(las$ExR >= ExR_thresh_min & las$ExR <= ExR_thresh_max &
                          las$ground == T &
                          las$Classification == LASNONCLASSIFIED, T, F)

poi_red_RPI <- ~if_else(las$RPI >= RPI_thresh & 
                          las$Classification == LASNONCLASSIFIED, T, F)

poi_veg_GLI <- ~if_else(las$GLI >= GLI_thresh & 
                          las$Classification == LASNONCLASSIFIED, T, F)

poi_veg_GPI <- ~if_else(las$GPI >= GPI_thresh & 
                          las$Classification == LASNONCLASSIFIED, T, F)

poi_veg_ExG <- ~if_else(las$ExG >= ExG_thresh & 
                          las$Classification == LASNONCLASSIFIED, T, F)

poi_veg_ExGR <- ~if_else(las$ExGR >= ExGR_thresh & 
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

poi_sed_times <- ~if_else(las$RBtimesGB >= RBtimesGB_min & las$RBtimesGB <= RBtimesGB_max &
                           las$ground == T &
                           las$Classification == LASNONCLASSIFIED, T, F)


# Read files--------------------------------------------------------------------

# empty warnings if existing.
if(length(warnings())!=0){
  assign("last.warning", NULL, envir = baseenv())
}

# Read LAS
# Intensity (i), color information (RGB), number of Returns (r), classification (c)
# of the first point is loaded only to reduce computational time.
las <- readLAS(data_path, select = "xyzRGBc", filter = cfg$las_filter)

aoi_shp <- read_sf(dsn = aoi_path)

# Filter points which are not within area of interest---------------------------
las <- classify_poi(las, class = LASNOISE, roi = aoi_shp, inverse_roi = T)
las <- filter_poi(las, Classification != LASNOISE)
# plot(las, size = 1, color = "RGB", bg = "white")

# Reset class LASNOISE for further procedure
las$Classification <- LASNONCLASSIFIED

# Create copy of read LAS to omit loading procedure.
# las_origin <- las
# las <- las_origin

# Check LAS whether it complies with the required
if(has.lasClassification(las)){
  print("LAS file is already classified. Are you sure to continue?")}
# if (length(warnings())>=1) {stop("The read LAS file throws warnings, script stops.")}

# Data exploration--------------------------------------------------------------
data_path

# las_origin <- las
# las <- las_origin

# Create attributes for classification------------------------------------------


# Generate list of active attributes 
active_attr <- names(las)
# Exclude some irrelevant attributes manually
active_attr <- active_attr[! active_attr %in% c("X","Y","Z","Classification", "RtoB", "RGtoB", "RBtimesGB", "Intensity" )]
active_attr

static_subtitle <- "Derivat aus Klassifikation"
las_post <- F

# map(active_attr, function(x){
#   gen.attribute.plot(las[[x]], x, output_id, static_subtitle, las_post, output_path)
# })


las_origin <- las
# las <- las_origin

# Segment Ground with Cloth Simulation Filter-----------------------------------
# setup csf filter settings
# rigidness: does not seem to have much impact.
# class_threshold and cloth_resolution influence each other. 0.5 x 0.5 is more conservative compared to 0.5 x 1.
# If sharper watercourse desired: increase threshold to 0.7 (leads to more canopy in ground points)


# las_origin <- las

# las <- las_origin

rigid_i <- 1:3
class_thres_i <- c(0.1,0.3,0.5,0.7,0.9)
cloth_res_i <- c(0.1,0.3,0.5,0.7,0.9)

for (i in class_thres_i) {
  for (j in cloth_res_i) {
    classify_gnd(las, i, j, 2)
  }
}


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

# Classify sky------------------------------------------------------------------
# Vegetation filter priority: ExB, BPI (some might be deactivated)

ExB_thresh <- cfg$sky_ExB_threshold
# tls: ExB = 18500, already takes away some cliff parts.
# tls: ExB = 14500, doubles cliff part but no sediment.
# tls: ExB = 13500, cliff wall is affected.
# tls: ExB = 3500, cliff wall is affected on a wide range but still works.
# uav: ExB between 3500 ... 6000, takes away sediment which is not ground = T

# BPI_thresh <- 0.391
# BPI = 0.391 represents 3rd Qu. and already takes away sediment

# las_origin <- las
# las <- las_origin

las <- classify_poi(las, class = LASBUILDING, poi = poi_sky_ExB)
# las <- filter_poi(las, Classification != LASBUILDING)
las_sky <- filter_poi(las, Classification == LASBUILDING)
# plot(las_sky, size = 1, color = "RGB", bg = "black")

# Classify green parts of vegetation--------------------------------------------
# Vegetation filter priority: GLI, ExG or GPI, ExGR (some might be deactivated)

GLI_thresh <- cfg$veg_GLI_threshold

las <- classify_poi(las, class = LASLOWVEGETATION, poi = poi_veg_GLI)
# las <- filter_poi(las, Classification != LASLOWVEGETATION)
las_veg <- filter_poi(las, Classification == LASLOWVEGETATION)
plot(las_veg, size = 1, color = "RGB", bg = "black")

# tls: GLI filters a broad range from greyish and yellowish parts.
# tls: GLI = 0.15 is conservative, no sediment and cliff is affected
# tls: GLI = 0.11 is ideal.
# tls: GLI = 0.09 boarder of watercourse is affected entirely.
# tls: GLI = 0.07 cliff and sediment points are affected.
# tls: GLI = 0.04 bush and canopy around sediment areas are included, except brown.
# uav: GLI between 0.06 ...(0.08)... 0.09 is ideal.

# GPI_thresh <- 0.35
# GPI filters a broad range from greyish and yellowish parts.
# GPI = 0.4 is conservative, no sediment and cliff is affected
# GPI = 0.37 is ideal.
# GPI = 0.36 boarder of watercourse is affected entirely.
# GPI = 0.35 is already the lower limit, cliff and sediment points are affected.

# ExG_thresh <- 4500
# ExG filters vegetation in general, neglects rather dark points
# ExG = 7500, border of watercourse is affected partly.
# ExG = 6000, border of watercourse is affected almost entirely
# ExG = 5500, border of watercourse and some sediment is affected.
# ExG = 4500, lower limit.

# ExGR_thresh <- 14000
# ExGR filters specially bright green and blue points, yellowish points not
# ExGR = 14000, Point of cliff are affected
# ExGR = 10000, Parts of cliff are affected


# Classify Red parts of vegetation----------------------------------------------
# Red filter priority: ExR, RPI (some might be deactivated)

ExR_thresh <- cfg$veg_ExR_threshold
# tls & uav: ExR = 30000, broad vegetation is affected, some sediment parts too.
# tls & uav: ExR = 20000, first lines of cliff relief is affected.

# RPI_thresh <- 0.9

las <- classify_poi(las, class = LASLOWVEGETATION, poi = poi_red_ExR)
# las <- filter_poi(las, Classification != LASBRIGDE)
# las_red <- filter_poi(las, Classification == LASBRIGDE)
# plot(las_red, size = 1, color = "RGB", bg = "black")

# Plot vegetation
# las_veg <- filter_poi(las, Classification == LASLOWVEGETATION)
# plot(las_veg, size = 1, color = "RGB", bg = "black")

# Classify sediment-------------------------------------------------------------
# las_origin <- las
# las <- las_origin

# Maximales RtoB: 1.21 in cliff_bright
# Minimales RtoB: 0.722 in cliff_dark kann aber zu sky zugeordnet werden.

#Min: 0.7222 and Max. 1.1290 derived from cliff_dark
RtoBmin <- cfg$sed_RtoBmin
RtoBmax <- cfg$sed_RtoBmax

#Min: 0.8333 and Max. 1.1613 derived from cliff_dark
#Min: 0.7928 from cliff_blue
GtoBmin <- cfg$sed_GtoBmin
GtoBmax <- cfg$sed_GtoBmax

# Approach with ratios RtoB and GtoB
las <- classify_poi(las, class = LASKEYPOINT, poi = poi_sed_ratios)
las_sed_ratios <- filter_poi(las, Classification == LASKEYPOINT)
# plot(las_sed_ratios, size = 1, color = "RGB", bg = "black")

# Classify band of negative excess red parts------------------------------------
# Negative ExR values appear to be sediment, ground criteria is set true.
ExR_thresh_max <- cfg$sed_ExR_max
ExR_thresh_min <- cfg$sed_ExR_min
# ExR = -53000 ... 1000 seems to affect sediment points more than others
# No big difference between -20000 and -53000

las <- classify_poi(las, class = LASKEYPOINT, poi = poi_sed_ExR_band)
las_sed <- filter_poi(las, Classification == LASKEYPOINT)

# For an exclusive plot of negative ExR values change to LASBRIDGE as class.
# las_sed_ratios <- filter_poi(las, Classification == LASBRIGDE)
plot(las_sed, size = 1, color = "RGB", bg = "white")

# Approach with RtoB times GtoB
# Set limits for sediment.
# RBtimesGB_min <- 0.9370
# RBtimesGB_max <- 1.0433
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
# RBtimesGB_min <- 0.98
# RBtimesGB_max <- 1.02
# las <- classify_poi(las, class = LASRAIL, poi = poi_rock_times)
# las_rock_times <- filter_poi(las, Classification == LASRAIL)
# plot(las_rock_times, size = 1, color = "RGB", bg = "black")

# Test - Classify band of negative excess blue parts----------------------------
# Negative ExB values appear to be sediment, ground criteria is set true.
# ExB_thresh_max <- 3499
# ExB_thresh_min <- -8000
# ExB = -8000 ... 3499 are greyish points but no distinct classes
# ExB = -60000 ... -8000 are greenish, greyish points but no distinct classes
# No big difference between -20000 and -53000

# las <- classify_poi(las, class = LASBRIGDE, poi = poi_sky_ExB_band)
# las_exb_neg <- filter_poi(las, Classification == LASBRIGDE)
# plot(las_exb_neg, size = 1, color = "RGB", bg = "white")
# summary(las$ExB)


# Generate subsets before filtering
las_foreveralone <- filter_poi(las, Classification == LASNONCLASSIFIED)

# Generate report of classification
sink(output_las_class_path)
print("Classified LAS in total:")
summary(las)
print("Classified black and white noise:")
summary(las_noise)
print("Classified sky points:")
summary(las_sky)
print("Classified ground points:")
summary(las_gnd)
print("Classified sediment points:")
summary(las_sed)
print("Classified vegetation points:")
summary(las_veg)
print("Remaining unclassified points:")
summary(las_foreveralone)

sink(append = T)

# Save generated output---------------------------------------------------------

writeLAS(las_sed, file = output_las_sed_path)
writeLAS(las_gnd, file = output_las_gnd_path)

# Filter out noise and unclassified points for a clean output file.
las <- filter_poi(las, Classification != LASNOISE)
las <- filter_poi(las, Classification != LASNONCLASSIFIED)
las <- filter_poi(las, Classification != LASUNCLASSIFIED)
writeLAS(las, file = output_las_all_path)

# Generate attribute plots after classification---------------------------------
las_post = T

map(active_attr, function(x){
  gen.attribute.plot(las[[x]], x, output_id, static_subtitle, las_post, output_path)
})

# Generate JSON report----------------------------------------------------------
end <- as_datetime(lubridate::now())
timespan <- interval(start, end)

run_time <- end - start
run_time

report <- cfg
report$timestamp <- timestamp
report$run_time_minutes <- as.numeric(timespan, "minutes")
report$data <- data_path
report$config_json <- config_json_path

json_report <- toJSON(report, indent = 1)
write(json_report, output_json_path, append = F)

# Generate DTM------------------------------------------------------------------
# TIN: fast and efficient, robust to empty regions. weak at edges.
# IDW: fast, not very realistic but good at edges. Compromise between TIN and Kriging.
# Kriging: very slow, not recommended for large areas.

# Class Nr. 8: LASKEYPOINT, here sediment. use "sfc" in shape for specific polygon boundaries.
# tin_sed <- rasterize_terrain(las_sed, res = cfg$DTM_resolution, algorithm = tin(), use_class = 8, shape = "convex")

# Generate DTM of ground points for comparison.
# tin_gnd <- rasterize_terrain(las_gnd, res = 0.45, algorithm = tin(), use_class = 2, shape = "convex")
 
# writeCDF(tin_sed, output_ncdf_path, overwrite = T)
# writeRaster(tin_sed, output_asc_path, overwrite = T)

# writeCDF(tin_gnd, output_ncdf_path, overwrite = T)
# writeRaster(tin_gnd, output_asc_path, overwrite = T)


# Plot classified point cloud---------------------------------------------------

# Show the unclassified---------------------------------------------------------
plot_dtm3d(tin_sed, bg = "white")
# plot_dtm3d(tin_gnd, bg = "white")

plot(las_foreveralone, size = 1, color = "RGB", bg = "black")

# disable rocks since it does not work properly
# las <- filter_poi(las, Classification != LASRAIL)
# plot(las, size = 1, color = "Classification", bg = "black")

# Plot separated classes
plot(las_sky, size = 1, color = "RGB", bg = "black")
plot(las_sed, size = 1, color = "RGB", bg = "white")

las_veg <- filter_poi(las, Classification == LASLOWVEGETATION)
plot(las_veg, size = 1, color = "RGB", bg = "black")

# Outdated stuff----------------------------------------------------------------

# Help lines
# Negation of attributes is also possible (all except intensity and angle)
# las = readLAS(LASfile, select = "* -i -a")
# Alternative filter functions--------------------------------------------------
# las_sub = filter_poi(las, Classification %in% c(LASGROUND, LASWATER))

# Prefer filter() befor filter_poi() since it does not read it on C++ level
# show all available filters
readLAS(filter = "-help")
