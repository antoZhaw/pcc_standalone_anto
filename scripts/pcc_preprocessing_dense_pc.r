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
library(rgl) # for RGL Viewer control functions

# Functions---------------------------------------------------------------------

has.lasRGB <- function(las) {
  has_R <- if_else(mean(las$R) != 0,T,F)
  has_G <- if_else(mean(las$G) != 0,T,F)
  has_B <- if_else(mean(las$B) != 0,T,F)
  valid <- has_R & has_G & has_B
  noRGB <- !valid
  return(noRGB)
}

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

set.RGLtopview <- function(x_scale = 800, y_scale = 800) {
  view3d(theta = 0, phi = 0, zoom = 0.6)
  par3d(windowRect = c(30, 30, x_scale, y_scale))
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

# Globals for Configuration-----------------------------------------------------
# Specify dataset
dataset_id <- "3"
wholeset <- T
year <- "2022"
perspective <- "tls"
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
if(user == "gubelyve"| user == "xgby"){
  dir_repo <- "C:/Daten/math_gubelyve/pcc_standalone"
  dir_data <- "C:/Daten/math_gubelyve"
} else{
  dir_repo <- "C:/code_wc/pcc_standalone"
  dir_data <- "C:/Daten/math_gubelyve"
}

dir_persp <- if_else(perspective == "tls", "tls_data", "uav_data")
dir_config <-  file.path(dir_repo, "config", fsep="/")

config_id <- as.character(paste(year, perspective, settype, dataset_id, sep = "-"))
output_id <- as.character(paste(timestamp, config_id, sep = "-"))
output_path <- file.path(dir_data, dir_persp, year, settype, "output", output_id, fsep="/")

config_json_name <- as.character(paste(config_id, ".json", sep = ""))
config_json_path <- file.path(dir_config, config_json_name, fsep="/")

# load dataset specific parameter set
cfg <- fromJSON(file = config_json_path)

# Create run specific output folder
dir.create(output_path)

aoi_shp <- as.character(paste(cfg$area_of_interest, ".shp", sep = ""))
aoi_path <- file.path(dir_repo, "data", cfg$area_of_interest, aoi_shp, fsep="/")

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

# Read files--------------------------------------------------------------------

# empty warnings if existing.
if(length(warnings())!=0){
  assign("last.warning", NULL, envir = baseenv())
}

# Read LAS
# Intensity (i), color information (RGB), number of Returns (r), classification (c)
# of the first point is loaded only to reduce computational time.
las <- readLAS(data_path, select = "xyzRGBc", filter = cfg$las_filter)
# aoi_shp <- read_sf(dsn = aoi_path)

# Save relevant information about the whole dataset.
sink(output_las_inp_path)
summary(las)
sink(append = T)

# Reset class LASNOISE for further procedure
las$Classification <- LASNONCLASSIFIED

# Create copy of read LAS to omit loading procedure.
# las_origin <- las
# las <- las_origin

# Check LAS whether it complies with the required
if (has.lasRGB(las)) 
  {stop("The read LAS file has no colour information, script stops.")}
if (!is.lasCRScompliant(las, cfg$crs_epsg))
  {stop("The read LAS file does not comply with the expected coordinate reference system, script stops.")}
if(has.lasClassification(las)){
  print("LAS file is already classified. Are you sure to continue?")}
# if (length(warnings())>=1) {stop("The read LAS file throws warnings, script stops.")}

# Display loaded classes (supposed to be zero)
# levels <- factor(las$Classification)

# Data exploration--------------------------------------------------------------
# plot(las, size = 1, color = "Intensity", bg = "black")
# plot(las, size = 1, color = "RGB", bg = "black")
# summary(las)
data_path

# las_origin <- las
# las <- las_origin

# Create attributes for classification------------------------------------------

# Add RGBmean attribute
las <- add_attribute(las, 0, "RGBmean")
las$RGBmean <- (las$R + las$G + las$B)/3

# Add Excess Blue Index ExB
las <- add_attribute(las, 0, "ExB")
las$ExB <- (2*las$B-las$R-las$G)

# Generate attribute plots before classification--------------------------------

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

# Classify noise in tls data----------------------------------------------------
if(perspective=="tls"){
# Classify outliers as noise. Can be skipped since it has high computational time.
  # if(cfg$outlier_already_filtered==F){
    # las <- classify_noise(las, sor(50,5))
    # las_outliers <- filter_poi(las, Classification %in% c(LASNOISE))
    # plot(las_outliers, size = 1, color = "RGB", bg = "white")
    # las <- filter_poi(las, Classification != LASNOISE)
  # }

# Classify white noise 
whitenoise_thresh <- cfg$whitenoise_threshold
# tls: good thresholds for white noise filter between 40000...(45000)...48000
# uav: good thresholds for white noise filter between 62000...65000

las <- classify_poi(las, class = LASNOISE, poi = poi_whitenoise)

# Classify black noise
blacknoise_thresh <- cfg$blacknoise_threshold
# tls: good thresholds for black noise filter between 4000...(6000)...8000
# uav: good thresholds for black noise filter between 8000...(10000)...12000

las <- classify_poi(las, class = LASNOISE, poi = poi_blacknoise)

# Classify sky------------------------------------------------------------------
# Vegetation filter priority: ExB, BPI (some might be deactivated)

ExB_thresh <- cfg$sky_ExB_threshold
# tls: ExB = 18500, already takes away some cliff parts.
# tls: ExB = 14500, doubles cliff part but no sediment.
# tls: ExB = 13500, cliff wall is affected.
# tls: ExB = 3500, cliff wall is affected on a wide range but still works.
# uav: ExB between 3500 ... 6000, takes away sediment which is not ground = T

las <- classify_poi(las, class = LASBUILDING, poi = poi_sky_ExB)
# las <- filter_poi(las, Classification != LASBUILDING)
las_sky <- filter_poi(las, Classification == LASBUILDING)
# plot(las_sky, size = 1, color = "RGB", bg = "black")

# Plot filtered noise
las_noise <- filter_poi(las, Classification == LASNOISE)
# plot(las_noise, size = 1, color = "RGB", bg = "black")

}
# las_origin <- las
# las <- las_origin

# Filter points which are not within area of interest---------------------------
# las <- classify_poi(las, class = LASNOISE, roi = aoi_shp, inverse_roi = T)
# las <- filter_poi(las, Classification != LASNOISE)
# plot(las, size = 1, color = "RGB", bg = "white")
# set.RGLtopview()

# Save relevant information about the area of interest.
sink(output_las_aoi_path)
summary(las)
sink(append = T)

# Generate report of classification
sink(output_las_class_path)
print("Remaining LAS for further processing:")
summary(las)
print("Classified black and white noise:")
summary(las_noise)
print("Classified sky points:")
summary(las_sky)
sink(append = T)

# Save generated output---------------------------------------------------------

# Filter out noise and unclassified points for a clean output file.
las <- filter_poi(las, Classification != LASNOISE)
las <- filter_poi(las, Classification != LASBUILDING)
writeLAS(las, file = output_las_all_path)

# Generate attribute plots after classification---------------------------------
# las_post = T
# 
# map(active_attr, function(x){
#   gen.attribute.plot(las[[x]], x, output_id, static_subtitle, las_post, output_path)
# })

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
