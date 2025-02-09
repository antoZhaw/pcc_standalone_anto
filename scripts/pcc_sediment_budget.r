#### Quantitative assessment of terrestrial sediment budgets in a ####
# hydropower-regulated floodplain using point cloud classification
# Documentation:  
# https://github.com/gubely/math_documentation/blob/main/main_math_documentation.pdf

# Author: Yves Gubelmann
# Contact: yves.gubelmann@gmx.ch

# Description: This script calculates a terrrestrial sediment budget based on two
# provided point clouds. First, terrestrial sediment is classified with cloth
# simulation filter (CSF) based on given .json-parameter. Second, an 2D habitat
# change map is calculated. Finally, assessed uncertainties are discarded from
# the calculated raster of difference (RoD) and displayed in an elevation change
# distribution (ECD). The calculated terrestrial volumes and areas are also saved
# as .csv-file for further analysis. Beforehand, the best CSF parameter were
# searched with a genetic algorithm. For detailed information, see:
# https://github.com/gubely/pcc_standalone/tree/master

# libraries ####
library(lidR) # Point cloud classification
library(tidyverse) # tidy essentials (ggplot, purr, tidyr, readr, dplyr)
library(lubridate)
library(tmap) # map visualization
# library(leaflet) # interactive maps
library(terra) # handling spatial data
library(sf) # handling spatial data
library(rgdal) # export raster objects
library(janitor) # clean and consistent naming
library(forcats) # handling factor levels
library(raster) # rasterizing vector data
library(fasterize) # for faster rasterization
library(geometry) # for raserize_canopy function
library(purrr) # for map function
library(rjson) # for JSON generation
library(rgl) # for RGL Viewer control functions
library(viridis) # for pretty color schemes
library(cowplot) # for multiple plot export

# Functions ####

# Classifies ground points based on given parameter
classify.gnd <- function(las, steep_slopes = F, class_thres, cloth_res, rigid) {
  mycsf <- csf(steep_slopes, class_thres, cloth_res, rigid)
  las <- classify_ground(las, mycsf)
  las_gnd <- filter_poi(las, Classification == LASGROUND)
  return(las_gnd)
}

# Set current RGL device to top view with a given scale
set.RGLtopview <- function(x_scale = 800, y_scale = 800) {
  view3d(theta = 0, phi = 0, zoom = 0.8)
  par3d(windowRect = c(30, 30, x_scale, y_scale))
}

# Creates a raster mask to subtract water surface from sediment points
mask.raster.layer <- function(raster_layer) {
  raster_layer[raster_layer != 0] <- 0
  raster_layer[is.na(raster_layer)] <- 1
  raster_layer
}

# Sets all values unlike NA to a normalised value
normalise.raster <- function(raster_layer, normal_value = 1) {
  raster_layer[raster_layer != 0] <- normal_value
  raster_layer[is.na(raster_layer)] <- 0
  raster_layer
}

# Filters a certain value and sets other values to zero
filter.raster <- function(raster_layer, filter_value = 1) {
  raster_layer[raster_layer != filter_value] <- 0
  raster_layer[raster_layer == filter_value] <- 1
  raster_layer[is.na(raster_layer)] <- 0
  raster_layer
}

# Sets a certain value to NA
value.to.na.raster <- function(raster_layer, na_value = 0) {
  raster_layer[raster_layer == na_value] <- NA
  raster_layer
}

# Sets all raster cells under the level of detection to NA
discard.uncertain.raster <- function(raster_layer, z_level_of_detection) {
  raster_layer[abs(raster_layer) < z_level_of_detection] <- NA
  raster_layer
}

# Gathers all uncertain raster cells under the value 10
gather.uncertain.raster <- function(raster_layer, z_level_of_detection) {
  raster_layer[abs(raster_layer) >= z_level_of_detection] <- NA
  raster_layer[!is.na(raster_layer)] <- 10
  raster_layer
}

# Create a summary of normalized raster based on a given resolution
summary.norm.raster <- function(raw_raster, res_m, summary_name) {
  summary <- data.frame(values(raw_raster)) %>% 
  filter(!is.na(values.raw_raster.)) %>%
  filter(values.raw_raster. != 0) %>%
  summarize(timestamp, summary_name, n=n(), area=n*res_m^2)
  summary
}

# Generate a plot of classified result and reference including aoi and layout variables
plot.csf.result.vs.target <- function(raster_bin, target_shp, aoi, plot_title, spec_layout, persp) {
  if(persp == "uav"){ # set customized plot settings for uav
    cp_position_csf <- c("right", "top")
    sb_breaks_csf <- c(0, 0.05, 0.1)
    sb_position_csf <- c(0.0, 0.2)
    sb_textsize_csf <- 0.8
  }else{ # set customized plot settings for tls
    cp_position_csf <- c(0, 0.73)
    sb_textsize_csf <- 1
    sb_position_csf <- c(0,0.65)
    sb_breaks_csf <- c(0, 0.01, 0.02)
  }
  bbox_aoi <- st_bbox(aoi) # crop with given bounding box
  palcsf <- c("#FFFFFF", "#2c7bb6")
  tmap_mode("plot") + # "plot" or "view"
  tm_shape(raster_bin, bbox = bbox_aoi) +
  tm_raster(title="", alpha = 1, palette = palcsf, style = "cat", 
            labels = c("Unclassified", "Classified area")) +
  tm_shape(target_shp) +
  tm_polygons(alpha = 0.65, lwd = 0.8, col = "#fdae61") +
  tm_shape(aoi) +
  tm_polygons(alpha = 0.0, lwd = 0.8, border.col = "#000000") +
  tm_view(control.position = c("right", "top")) +
  tm_compass(type = "arrow", size = 2, position = cp_position_csf) +
  tm_scale_bar(breaks = sb_breaks_csf, text.size = sb_textsize_csf, position = sb_position_csf) +
  tm_layout(main.title = plot_title) +
  tm_add_legend('fill', 
                  col = "#fdae61",  alpha = 0.6,
                  labels = c('Reference data')) +
  tm_add_legend('fill', 
                  border.col = "#000000",
                  col = "#ffffff",
                  labels = c('Area of interest')) +
  spec_layout
}

# Create classes of a common sediment budget based on a classified raster and lod_critical
create.budget.classes <- function(raw_raster, lod_critical, raster_res) {
  dt <- data.frame(values(raw_raster)) %>% 
    mutate(discarded = if_else(abs(values.raw_raster.) <= lod_critical, T, F),
           class = case_when(
             discarded == T & values.raw_raster. <= 0 ~ "Discarded Erosion",
             discarded == T & values.raw_raster. > 0 ~ "Discarded Deposition",
             discarded == F & values.raw_raster. <= 0 ~ "Erosion",
             discarded == F & values.raw_raster. > 0 ~ "Deposition",
             TRUE ~ NA # Default case
           ),
           cell_vol = raster_res^2*values.raw_raster.) # calculate volume per cell
  dt
}

# Globals for Configuration ####
# Record start date and time
start <- as_datetime(lubridate::now())
date <- as.Date(start)
hour <- hour(start)
minute <- minute(start)
# Generate timestamp for this run.
timestamp <- as.character(paste(date, hour, minute, sep = "-"))

# Settings which are unique for t0 and t1
# Therefore, such settings are stored here instead of a single-timestep-json
wholeset <- T
settype <- if_else(wholeset == T, "wholeset", "subset")
comment <- "narrow breaks, full saturation"
narrow_breaks <- c(-1, -0.5, 0.5, 1)
wide_breaks <- c(-2, -1, 1, 2)
global_breaks <- narrow_breaks # choose breaks for ECD
sat_basemap <- 1 # saturation suggestion: 1 or 0
alpha_basemap <- 0.3 # alpha suggestion: 0.3 or 0.35

# Settings t0 and t1
# uav 2022-2021
# perspective <- "uav"
# flood_startdate <- "31.05.2022"
# flood_prefix <- "310522"
# tif_path <- "C:/Daten/math_gubelyve/tiff_data/310522_bg.tif"
# aggr_factor <- 18
# t0_dataset_id <- "1"
# t0_year <- "2021"
# t1_dataset_id <- "1"
# t1_year <- "2022"
# raster_res <- 0.4

# uav 2021-2020
# perspective <- "uav"
# flood_startdate <- "11.07.2021"
# flood_prefix <- "110721"
# tif_path <- "C:/Daten/math_gubelyve/tiff_data/110721_bg.tif"
# aggr_factor <- 18
# t0_dataset_id <- "1"
# t0_year <- "2020"
# t1_dataset_id <- "1"
# t1_year <- "2021"
# raster_res <- 0.4

# uav 2020-2020
# perspective <- "uav"
# flood_startdate <- "22.10.2020"
# flood_prefix <- "221020"
# tif_path <- "C:/Daten/math_gubelyve/tiff_data/221020_bg.tif"
# aggr_factor <- 5
# t0_dataset_id <- "2"
# t0_year <- "2020"
# t1_dataset_id <- "1"
# t1_year <- "2020"
# raster_res <- 0.4

# uav overall
# perspective <- "uav"
# flood_startdate <- "NA"
# flood_prefix <- "overall"
# tif_path <- "C:/Daten/math_gubelyve/tiff_data/310522_bg.tif"
# aggr_factor <- 18
# t0_dataset_id <- "2"
# t0_year <- "2020"
# t1_dataset_id <- "1"
# t1_year <- "2022"
# raster_res <- 0.4

# tls 2022-2021
flood_startdate <- "31.05.2022"
flood_prefix <- "310522"
tif_path <- "C:/Daten/math_gubelyve/tiff_data/310522_bg.tif"
aggr_factor <- 18
perspective <- "tls"
t0_dataset_id <- "4"
t0_year <- "2021"
t1_dataset_id <- "4"
t1_year <- "2022"
raster_res <- 0.2

# Generate static tif as backgroud
e <- extent(2575009, 2575489, 1178385, 1178900) # set manual extent for tif file
tif <- terra::rast(x=tif_path)
tot_aoi <- raster(crs=2056, ext=e, resolution=0.2, vals=NULL) # create raster template
tif_aoi <- crop(tif, tot_aoi) # crop with aoi
tif_crop <- terra::aggregate(tif_aoi, aggr_factor) # downsample tif to reasonable size

# Load environment dependent paths.
user <- Sys.getenv("USERNAME")
if(user == "gubelyve"| user == "anto"){
  dir_repo <- "C:/Daten/math_gubelyve/pcc_standalone"
  dir_data <- "C:/Daten/math_gubelyve"
} else{
  dir_repo <- "C:/code_wc/pcc_standalone"
  dir_data <- "C:/code_wc/math_gubelyve"
}
dir_persp <- if_else(perspective == "tls", "tls_data", "uav_data")

# Generate dynamic paths for time step t0
t0_datasetname <- as.character(paste(t0_year, perspective, settype, t0_dataset_id, sep = "-"))
t0_dataset <- paste(t0_datasetname, ".las", sep = "")
t0_dir_config <-  file.path(dir_repo, "config", fsep="/")
t0_config_id <- as.character(paste(t0_year, perspective, settype, t0_dataset_id, sep = "-"))
t0_output_id <- as.character(paste(timestamp, t0_config_id, sep = "-"))
t0_output_path <- file.path(dir_data, dir_persp, t0_year, settype, "output", t0_output_id, fsep="/")
t0_config_json_name <- as.character(paste(t0_config_id, ".json", sep = ""))
t0_config_json_path <- file.path(t0_dir_config, t0_config_json_name, fsep="/")
t0_data_path <- file.path(dir_data, dir_persp, t0_year, settype, t0_dataset)
# load dataset specific parameter set
t0_cfg <- fromJSON(file = t0_config_json_path)
t0_mcdut_shp_name <- as.character(paste(t0_cfg$mapcurve_dut_shp, sep = ""))
t0_mcdut_path <- file.path(dir_repo, t0_mcdut_shp_name, fsep="/")
t0_mctar_shp_name <- as.character(paste(t0_cfg$mapcurve_target_shp, sep = ""))
t0_mctar_path <- file.path(dir_repo, t0_mctar_shp_name, fsep="/")

# Generate dynamic paths for time step t1
t1_datasetname <- as.character(paste(t1_year, perspective, settype, t1_dataset_id, sep = "-"))
t1_dataset <- paste(t1_datasetname, ".las", sep = "")
t1_dir_config <-  file.path(dir_repo, "config", fsep="/")
t1_config_id <- as.character(paste(t1_year, perspective, settype, t1_dataset_id, sep = "-"))
t1_output_id <- as.character(paste(timestamp, t1_config_id, sep = "-"))
t1_output_path <- file.path(dir_data, dir_persp, t1_year, settype, "output", t1_output_id, fsep="/")
t1_config_json_name <- as.character(paste(t1_config_id, ".json", sep = ""))
t1_config_json_path <- file.path(t1_dir_config, t1_config_json_name, fsep="/")
t1_data_path <- file.path(dir_data, dir_persp, t1_year, settype, t1_dataset)
# load dataset specific parameter set
t1_cfg <- fromJSON(file = t1_config_json_path)
t1_mcdut_shp_name <- as.character(paste(t1_cfg$mapcurve_dut_shp, sep = ""))
t1_mcdut_path <- file.path(dir_repo, t1_mcdut_shp_name, fsep="/")
t1_mctar_shp_name <- as.character(paste(t1_cfg$mapcurve_target_shp, sep = ""))
t1_mctar_path <- file.path(dir_repo, t1_mctar_shp_name, fsep="/")

# Load shapefiles
aoi_path <- file.path(dir_repo, "data/aoi_uav_230227", "aoi_uav_230227.shp", fsep="/")
aoi_dir <- file.path(dir_repo, "data/aoi_uav_230227", fsep="/")
aoi_path <- file.path(dir_repo, "data/dut_filter_230315", "dut_filter_230315.shp", fsep="/")
aoi_dir <- file.path(dir_repo, "data/dut_filter_230315", fsep="/")

# Create run specific output folder
bud_output_path <- file.path(dir_data, dir_persp, "budget", timestamp, fsep="/")
dir.create(bud_output_path)
dir.create(t0_output_path)
dir.create(t1_output_path)
output_bud_report_name <- as.character(paste(flood_prefix, "-budget-report.txt", sep = ""))
output_bud_report_path <- file.path(bud_output_path, output_bud_report_name, fsep="/")

# Read point cloud and shapefiles ####
# empty warnings if existing.
if(length(warnings())!=0){
  assign("last.warning", NULL, envir = baseenv())
}

# Generate rectangular polygon of area of interest
if(perspective == "uav"){ # generate uav specific bounding box
gen_xy <- structure(list(dat = c("BURR-1-1-1", "BURR-1-1-1", 
                                 "BURR-1-1-1", "BURR-1-1-1"),
                         Longitude = c(2575011, 2575011, 
                                      2575517.4, 2575517.4),
                         Latitude = c(1178366.7, 1178993, 
                                       1178993, 1178366.7)),
                    class = "data.frame", row.names = c(NA,-4L))
sb_breaks <- c(0, 0.05, 0.1)
sb_textsize <- 0.9
sb_position <- c(0.0, 0.2)
}else{ # generate tls specific bounding box
gen_xy <- structure(list(dat = c("AOI TLS", "AOI TLS", 
                                     "AOI TLS", "AOI TLS"),
                             Longitude = c(2575340, 2575340, 
                                           2575480, 2575480),
                             Latitude = c(1178520, 1178780, 
                                          1178780, 1178520)),
                        class = "data.frame", row.names = c(NA,-4L))
sb_breaks <- c(0, 0.02, 0.04)
sb_textsize <- 0.7
sb_position <- c(0.0, 0.48)
}
bounding_box <- gen_xy %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 2056) %>%
  group_by(dat) %>%
  summarise(geometry = st_combine(geometry)) %>%
  st_cast("POLYGON")

# Read shapefiles
t0_csf_aoi_shp <- read_sf(dsn = t0_mcdut_path)
t1_csf_aoi_shp <- read_sf(dsn = t1_mcdut_path)
#Mask tif with aoi
tif_masked <- terra::mask(tif_crop, t1_csf_aoi_shp, inverse = F)

if(perspective == "tls"){
  t0_csf_aoi_shp <- st_intersection(t0_csf_aoi_shp, bounding_box)
  t1_csf_aoi_shp <- st_intersection(t0_csf_aoi_shp, bounding_box)
}

# Read and generate targets
t0_mctar_all <- read_sf(dsn = t0_mctar_path)
t1_mctar_all <- read_sf(dsn = t1_mctar_path)

if(t1_year == "2022"){
  t1_mctar_all <- t1_mctar_all %>% 
    mutate(
      Id = case_when(
        Class_name == "Water"~1,
        Class_name == "Sediment"~2,
        Class_name == "Vegetation"~3,
        Class_name == "Other"~4,
        Class_name == "Bedrock/Boulders"~5,
        TRUE~99 #Default case
      )
    )
}

# Intersect target with area of interest
t0_targets_aoi_shp <- st_intersection(t0_mctar_all, bounding_box)
t1_targets_aoi_shp <- st_intersection(t1_mctar_all, bounding_box)

t0_target_sed <- t0_targets_aoi_shp %>%  filter(Id == 2)
t1_target_sed <- t1_targets_aoi_shp %>%  filter(Id == 2)

if(perspective == "uav"){
t1_target_wat <- t1_targets_aoi_shp %>%  filter(Id == 1)
t0_target_wat <- t0_targets_aoi_shp %>%  filter(Id == 1)
}

# Read LAS
# Intensity (i), color information (RGB), number of Returns (r), classification (c)
# of the first point is loaded only to reduce computational time.
t0_las <- readLAS(t0_data_path, select = "xyzRGBc", filter = t0_cfg$las_filter)
t1_las <- readLAS(t1_data_path, select = "xyzRGBc", filter = t1_cfg$las_filter)
# generate a text-file with a summary of the loaded point cloud before classification
sink(output_bud_report_path)
print("Summary of LAS t0:")
print(t0_config_id)
summary(t0_las)
print("Summary of LAS t1:")
print(t1_config_id)
summary(t1_las)
sink(append = T)

# Reset class flag for classification
t0_las$Classification <- LASNONCLASSIFIED
t1_las$Classification <- LASNONCLASSIFIED

# Classify Water with Cloth Simulation Filter ####
if(perspective == "uav"){ # water surface is only apparent for uav data
# load best csf parameter settings from config file (.json)
t0_rigid_n_water <- t0_cfg$csf_water_rigidness
t1_rigid_n_water <- t1_cfg$csf_water_rigidness
t0_class_thres_water <- t0_cfg$csf_water_class_threshold
t1_class_thres_water <- t1_cfg$csf_water_class_threshold
t0_cloth_res_water <- t0_cfg$csf_water_cloth_resolution
t1_cloth_res_water <- t1_cfg$csf_water_cloth_resolution
t0_steep_slope_water <- if_else(t0_rigid_n_water == 3, F, T)
t1_steep_slope_water <- if_else(t1_rigid_n_water == 3, F, T)
# Generate filename for classified result
t0_status_water <- as.character(paste("00_water", t0_year, "_rigid", t0_rigid_n_water, 
                                      "_clthres", t0_class_thres_water, "clothres", t0_cloth_res_water, sep = ""))
t1_status_water <- as.character(paste("00_water", t1_year, "_rigid", t1_rigid_n_water,
                                      "_clthres", t1_class_thres_water, "clothres", t1_cloth_res_water, sep = ""))
# Classify water surface
t0_las_water <- classify.gnd(t0_las, t0_steep_slope_water, t0_class_thres_water, t0_cloth_res_water, t0_rigid_n_water)
t1_las_water <- classify.gnd(t1_las, t1_steep_slope_water, t1_class_thres_water, t1_cloth_res_water, t1_rigid_n_water)
# Add classification to additional attribute
t0_las_water <- add_attribute(t0_las_water, FALSE, "water")
t1_las_water <- add_attribute(t1_las_water, FALSE, "water")
t0_las_water$water <- if_else(t0_las_water$Classification == LASGROUND, T, F)
t1_las_water$water <- if_else(t1_las_water$Classification == LASGROUND, T, F)
# Create a subset for visualisation purposes
t0_las_water <- filter_poi(t0_las_water, Classification == LASGROUND)
t1_las_water <- filter_poi(t1_las_water, Classification == LASGROUND)
# plot(las_water, size = 1, color = "RGB", bg = "white", axis = F)

# Filter points which are not within area of interest
t0_las_water <- classify_poi(t0_las_water, class = LASNOISE, roi = t0_csf_aoi_shp, inverse_roi = T)
t0_las_water <- filter_poi(t0_las_water, Classification != LASNOISE)
# plot(t0_las_water, size = 1, color = "RGB", bg = "black", axis = F)
# set.RGLtopview()
# t0_output_water_png_name <- as.character(paste(t0_status_water, ".png", sep = ""))
# t0_output_water_png_path <- file.path(bud_output_path, t0_output_water_png_name, fsep="/")
# rgl.snapshot(t0_output_water_png_path)
# rgl.close()

t1_las_water <- classify_poi(t1_las_water, class = LASNOISE, roi = t1_csf_aoi_shp, inverse_roi = T)
t1_las_water <- filter_poi(t1_las_water, Classification != LASNOISE)
# plot(t1_las_water, size = 1, color = "RGB", bg = "black", axis = F)
# set.RGLtopview()
# t1_output_water_png_name <- as.character(paste(t1_status_water, ".png", sep = ""))
# t1_output_water_png_path <- file.path(bud_output_path, t1_output_water_png_name, fsep="/")
# rgl.snapshot(t1_output_water_png_path)
# rgl.close()

# Rasterize water point cloud ####
t0_DEM_water <- rasterize_canopy(t0_las_water, res = raster_res, p2r(), pkg = "raster")
t1_DEM_water <- rasterize_canopy(t1_las_water, res = raster_res, p2r(), pkg = "raster")
t0_mask_water <- mask.raster.layer(t0_DEM_water)
t1_mask_water <- mask.raster.layer(t1_DEM_water)
}

# Classify sediment with Cloth Simulation Filter ####
# load best csf parameter settings from config file (.json)
t0_rigid_n_sed <- t0_cfg$csf_gnd_rigidness
t1_rigid_n_sed <- t1_cfg$csf_gnd_rigidness
t0_class_thres_sed <- t0_cfg$csf_gnd_class_threshold
t1_class_thres_sed <- t1_cfg$csf_gnd_class_threshold
t0_cloth_res_sed <- t0_cfg$csf_gnd_cloth_resolution
t1_cloth_res_sed <- t1_cfg$csf_gnd_cloth_resolution
t0_steep_slope_sed <- if_else(t0_rigid_n_sed == 3, F, T)
t1_steep_slope_sed <- if_else(t1_rigid_n_sed == 3, F, T)
# Classify ground
t0_las_sed <- classify.gnd(t0_las, t0_steep_slope_sed, t0_class_thres_sed, t0_cloth_res_sed, t0_rigid_n_sed)
t0_las_sed <- classify_poi(t0_las_sed, class = LASNOISE, roi = t0_csf_aoi_shp, inverse_roi = T)
t0_las_sed <- filter_poi(t0_las_sed, Classification != LASNOISE)
t1_las_sed <- classify.gnd(t1_las, t1_steep_slope_sed, t1_class_thres_sed, t1_cloth_res_sed, t1_rigid_n_sed)
t1_las_sed <- classify_poi(t1_las_sed, class = LASNOISE, roi = t0_csf_aoi_shp, inverse_roi = T)
t1_las_sed <- filter_poi(t1_las_sed, Classification != LASNOISE)

# Rasterize sediment point cloud ####
t0_DEM_sed <- rasterize_canopy(t0_las_sed, res = raster_res, p2r(), pkg = "raster")
t1_DEM_sed <- rasterize_canopy(t1_las_sed, res = raster_res, p2r(), pkg = "raster")

# Subtract water surface. this is only required for uav data
if(perspective == "uav"){
  t0_sed <- t0_DEM_sed * t0_mask_water
  t1_sed <- t1_DEM_sed * t1_mask_water
  t0_tm_wat <- normalise.raster(t0_DEM_water)
  t1_tm_wat <- normalise.raster(t1_DEM_water)
}else{
  t0_sed <- t0_DEM_sed
  t1_sed <- t1_DEM_sed
}
t0_tm_sed <- normalise.raster(t0_sed)
t1_tm_sed <- normalise.raster(t1_sed)

# Plot comparison between target and classified raster ####
output_gof_t0_sed_name <- as.character(paste(flood_prefix, "_gof_t0_sed.png", sep = ""))
output_gof_t0_sed_path <- file.path(bud_output_path, output_gof_t0_sed_name, fsep="/")
output_gof_t1_sed_name <- as.character(paste(flood_prefix, "_gof_t1_sed.png", sep = ""))
output_gof_t1_sed_path <- file.path(bud_output_path, output_gof_t1_sed_name, fsep="/")
if(perspective == "uav"){
  output_gof_t0_wat_name <- as.character(paste(flood_prefix, "_gof_t0_wat.png", sep = ""))
  output_gof_t0_wat_path <- file.path(bud_output_path, output_gof_t0_wat_name, fsep="/")
  output_gof_t1_wat_name <- as.character(paste(flood_prefix, "_gof_t1_wat.png", sep = ""))
  output_gof_t1_wat_path <- file.path(bud_output_path, output_gof_t1_wat_name, fsep="/")
}
# Generate paths for output files
output_hab_change_name <- as.character(paste(flood_prefix, "_habitate_change.png", sep = ""))
output_hab_change_path <- file.path(bud_output_path, output_hab_change_name, fsep="/")
output_hab_change_bg_name <- as.character(paste(flood_prefix, "_habitate_change_bg.png", sep = ""))
output_hab_change_bg_path <- file.path(bud_output_path, output_hab_change_bg_name, fsep="/")
output_hab_change_comp_bg_name <- as.character(paste(flood_prefix, "_habitate_change_comp_bg.png", sep = ""))
output_hab_change_comp_bg_path <- file.path(bud_output_path, output_hab_change_comp_bg_name, fsep="/")
output_elev_name <- as.character(paste(flood_prefix, "_elevation_change.png", sep = ""))
output_elev_path <- file.path(bud_output_path, output_elev_name, fsep="/")
output_elev_bg_name <- as.character(paste(flood_prefix, "_elevation_change_bg.png", sep = ""))
output_elev_bg_path <- file.path(bud_output_path, output_elev_bg_name, fsep="/")
output_elev_uncert_name <- as.character(paste(flood_prefix, "_elevation_change_uncert.png", sep = ""))
output_elev_uncert_path <- file.path(bud_output_path, output_elev_uncert_name, fsep="/")
output_elev_uncert_bg_name <- as.character(paste(flood_prefix, "_elevation_change_uncert_bg.png", sep = ""))
output_elev_uncert_bg_path <- file.path(bud_output_path, output_elev_uncert_bg_name, fsep="/")
output_lod_hist_name <- as.character(paste(flood_prefix, "_lod_histogram.png", sep = ""))
output_lod_hist_path <- file.path(bud_output_path, output_lod_hist_name, fsep="/")
output_lod_hist500_name <- as.character(paste(flood_prefix, "_lod_histogram_y500.png", sep = ""))
output_lod_hist500_path <- file.path(bud_output_path, output_lod_hist500_name, fsep="/")
output_lod_bar_name <- as.character(paste(flood_prefix, "_lod_barplot.png", sep = ""))
output_lod_bar_path <- file.path(bud_output_path, output_lod_bar_name, fsep="/")
output_budget_name <- as.character(paste(flood_prefix, "_budget-overview.png", sep = ""))
output_budget_path <- file.path(bud_output_path, output_budget_name, fsep="/")
output_budget500_name <- as.character(paste(flood_prefix, "_budget500-overview.png", sep = ""))
output_budget500_path <- file.path(bud_output_path, output_budget500_name, fsep="/")
# Set specific layout options
if(perspective == "uav"){
  gof_layout <- tm_layout(frame = F, legend.text.size = 1.0, legend.title.size = 1.3, legend.outside = F, legend.position = c("left", "center"),
                        main.title.position = "center", main.title.size = 1.3)
}else{
  gof_layout <- tm_layout(frame = F, legend.text.size = 0.65, legend.title.size = 0.9, legend.outside = F, legend.position = c(0.9, 0.0),
                          main.title.position = "center", main.title.size = 0.9)
}
# Save map of classified sediment results versus reference data
t0_title_sed <- paste("Sediment Classification", t0_cfg$survey_date_pret, sep = " ")
t0_tm_sed_result <- plot.csf.result.vs.target(t0_tm_sed, t0_target_sed, t0_csf_aoi_shp, t0_title_sed, gof_layout, perspective)
tmap_save(tm = t0_tm_sed_result, output_gof_t0_sed_path, width = 1920, height = 1920)
t1_title_sed <- paste("Sediment Classification", t1_cfg$survey_date_pret, sep = " ")  
t1_tm_sed_result <- plot.csf.result.vs.target(t1_tm_sed, t1_target_sed, t1_csf_aoi_shp, t1_title_sed, gof_layout, perspective)
tmap_save(tm = t1_tm_sed_result, output_gof_t1_sed_path, width = 1920, height = 1920)
# Save map of classified water results versus reference data
if(perspective == "uav"){
  t0_title_wat <- paste("Water Classification", t0_cfg$survey_date_pret, sep = " ")  
  t0_tm_wat_result <- plot.csf.result.vs.target(t0_tm_wat, t0_target_wat, t0_csf_aoi_shp, t0_title_wat, gof_layout, perspective)
  tmap_save(tm = t0_tm_wat_result, output_gof_t0_wat_path, width = 1920, height = 1920)
  t1_title_wat <- paste("Water Classification", t1_cfg$survey_date_pret, sep = " ")  
  t1_tm_wat_result <- plot.csf.result.vs.target(t1_tm_wat, t1_target_wat, t1_csf_aoi_shp, t1_title_wat, gof_layout, perspective)
  tmap_save(tm = t1_tm_wat_result, output_gof_t1_wat_path, width = 1920, height = 1920)
}

# Determine habitate change ####
t0_tm_hab <- normalise.raster(t0_sed, 10) # set sediment cells to value 10
t1_tm_hab <- t1_tm_sed
# Generate 2D habitat change map 
tm_habitate <- t0_tm_hab + t1_tm_hab # addition creates pseudo factors (0, 1, 10, 11)
# Generate summary of 2D habitat change map as subproduct (Areas only)
if(perspective == "uav"){
  t0_tm_sedsum <- summary.norm.raster(t0_tm_sed, raster_res, "t0_tm_sed")
  t1_tm_sedsum <- summary.norm.raster(t1_tm_sed, raster_res, "t1_tm_sed")
  t0_tm_watsum <- summary.norm.raster(t0_tm_wat, raster_res, "t0_tm_wat")
  t1_tm_watsum <- summary.norm.raster(t1_tm_wat, raster_res, "t1_tm_wat")
  single_area_results <- rbind(t0_tm_sedsum, t1_tm_sedsum, t0_tm_watsum, t1_tm_watsum)
  bbox_aoi <- st_bbox(t0_csf_aoi_shp)
  tm_default_layout <- tm_layout(frame = F, 
                                 legend.title.size = 1.3, legend.text.size = 1.0, 
                                 legend.outside = F, legend.position = c("left", "center"),
                                 main.title.position = "center", main.title.size = 1.3)
}else{
  t0_tm_sedsum <- summary.norm.raster(t0_tm_sed, raster_res, "t0_tm_sed")
  t1_tm_sedsum <- summary.norm.raster(t1_tm_sed, raster_res, "t1_tm_sed")
  single_area_results <- rbind(t0_tm_sedsum, t1_tm_sedsum)
  gen_xy_tls <- structure(list(dat = c("AOI TLS", "AOI TLS", 
                                       "AOI TLS", "AOI TLS"),
                               Longitude = c(2575310, 2575310, 
                                             2575480, 2575480),
                               Latitude = c(1178520, 1178780, 
                                            1178780, 1178520)),
                          class = "data.frame", row.names = c(NA,-4L))
  bbox_aoi <- gen_xy_tls %>%
    st_as_sf(coords = c("Longitude", "Latitude"), crs = 2056) %>%
    group_by(dat) %>%
    summarise(geometry = st_combine(geometry)) %>%
    st_cast("POLYGON")
  tm_default_layout <- tm_layout(frame = F, 
                                 legend.title.size = 1.0, legend.text.size = 0.8, 
                                 legend.outside = F, legend.position = c(0.0, 0.55),
                                 main.title.position = "center", main.title.size = 1.0)
}
# Save generated summary as csv
write.table(single_area_results, file = "C:/Daten/math_gubelyve/pcc_standalone/export/budget_single_area_results.csv",
            append = T, sep = ";", row.names = F, col.names = F)
# Plot and save 2D habitat change map
pal4div <- c("#FFFFFF", "#440154", "#fde725", "#31688e")
hab_title <- paste(t1_cfg$survey_date_pret, "-", t0_cfg$survey_date_pret, sep = " ")
tm_hab <- tmap_mode("plot") + # "plot" or "view"
  tm_shape(tm_habitate, bbox = bbox_aoi) +
  tm_raster(title="", palette = pal4div, alpha = 1, style = "cat", breaks = c(0, 2, 10.5),
            labels = c("Unchanged", "New deposition", "New erosion ", "Change in elevation")) +
  tm_shape(t0_csf_aoi_shp) +
  tm_polygons(alpha = 0.0, lwd =0.8, border.col = "#000000") +
  # tm_layout(main.title = hab_title) +
  tm_add_legend(type='fill', border.col = "#000000", col = "#ffffff", labels = c('Area of interest')) +
  tm_compass(type = "arrow", size = 2, position = c("right", "top")) +
  tm_scale_bar(breaks = sb_breaks, text.size = sb_textsize, position = sb_position) +
  tm_default_layout
tmap_save(tm = tm_hab, output_hab_change_path, width = 1920, height = 1920)
# Add basemap to 2D habitat change map
tm_hab_bg <- tm_hab +
  tm_shape(tif_masked) +
  tm_rgb(r=1, g=2, b=3, alpha = alpha_basemap, saturation = sat_basemap)
tmap_save(tm = tm_hab_bg, output_hab_change_bg_path, width = 1920, height = 1920)
# Add plot for comparison with other samples
tm_hab_bg_comp <- tm_hab +
  tm_shape(t0_target_sed) +
  tm_polygons(alpha = 0.5, lwd = 0.8, col = "#d73027") +
  tm_shape(t1_target_sed) +
  tm_polygons(alpha = 0.5, lwd = 0.8, col = "#a50026") +
  tm_shape(tif_masked) +
  tm_rgb(r=1, g=2, b=3, alpha = alpha_basemap, saturation = sat_basemap) +
  tm_add_legend('fill', col = "#a50026",  alpha = 0.5, labels = c('Sediment 2022')) +
  tm_add_legend('fill', border.col = "#000000", col = "#d73027", alpha = 0.5,
                labels = c('Sediment 2021'))
tm_hab_bg_comp
tmap_save(tm = tm_hab_bg_comp, output_hab_change_comp_bg_path, width = 1920, height = 1920)

# Calculate Raster of Difference (RoD) ####
# Generate mask for cells which show a change in elevation (pick value 11)
tm_elevation_mask <- filter.raster(tm_habitate, 11)
# Set zero values to na 
tm_elevation_mask <- value.to.na.raster(tm_elevation_mask)
# Apply mask on delta z
delta_z_all <- tm_elevation_mask*(t1_sed - t0_sed)

# Plot elevation change without uncertainty assessment
elev_title <- paste(t1_cfg$survey_date_pret, "-", t0_cfg$survey_date_pret, sep = " ")
tm_elev <- tmap_mode("plot") + # "plot" or "view"
  tm_shape(delta_z_all, bbox = bbox_aoi) +
  tm_raster(title = "Elevation change [m]", alpha = 1, style = "cont", palette = "RdBu", breaks = global_breaks) + 
  tm_shape(t0_csf_aoi_shp) +
  tm_polygons(alpha = 0.0, lwd = 0.8, border.col = "#000000") +
  # tm_layout(main.title = elev_title) +
  tm_add_legend('fill', border.col = "#000000", col = "#ffffff", labels = c('Area of interest')) +
  tm_compass(type = "arrow", size = 2, position = c("right", "top")) +
  tm_scale_bar(breaks = sb_breaks, text.size = sb_textsize, position = sb_position) +
  tm_default_layout
tmap_save(tm = tm_elev, output_elev_path, width = 1920, height = 1920)
# Add basemap to RoD
tm_elev_bg <- tm_elev +
  tm_shape(tif_masked) +
  tm_rgb(r=1, g=2, b=3, alpha = alpha_basemap, saturation = sat_basemap)
tmap_save(tm = tm_elev_bg, output_elev_bg_path, width = 1920, height = 1920)

# Calculate critical level of detection ####
lod_crit <- 1.96*sqrt(t1_cfg$z_sigma_estimated^2 + t0_cfg$z_sigma_estimated^2)
# Propagate assessed uncertainties into RoD 
delta_z_cleaned <- discard.uncertain.raster(delta_z_all, lod_crit)
delta_z_noise <- gather.uncertain.raster(delta_z_all, lod_crit)

# Plot elevation change with uncertainty assessment
paldisc <- c("#000000")
elev_uncert_title <- paste(t1_cfg$survey_date_pret, "-", t0_cfg$survey_date_pret, sep = " ")
tm_elev_uncert <- tmap_mode("plot") + # "plot" or "view"
  tm_shape(delta_z_cleaned, bbox = bbox_aoi) +
  tm_raster(title = "Elevation change [m]", alpha = 1, style = "cont", palette = "RdBu", breaks = global_breaks) + 
  tm_shape(delta_z_noise) +
  tm_raster(title = "", palette = paldisc, alpha = 1, style = "cont", labels = c("Discarded cells")) +
  tm_shape(t0_csf_aoi_shp) +
  tm_polygons(alpha = 0.0, lwd = 0.8, border.col = "#000000") +
  # tm_layout(main.title = elev_uncert_title) +
  tm_add_legend('fill', border.col = "#000000", col = "#ffffff", labels = c('Area of interest')) +
  tm_compass(type = "arrow", size = 2, position = c("right", "top")) +
  tm_scale_bar(breaks = sb_breaks, text.size = sb_textsize, position = sb_position) +
  tm_default_layout
tm_elev_uncert
tmap_save(tm = tm_elev_uncert, output_elev_uncert_path, width = 1920, height = 1920)
# Add basemap to RoD with propagated uncertainties
tm_elev_uncert_bg <- tm_elev_uncert +
  tm_shape(tif_masked) +
  tm_rgb(r=1, g=2, b=3, alpha = alpha_basemap, saturation = sat_basemap)
tmap_save(tm = tm_elev_uncert_bg, output_elev_uncert_bg_path, width = 1920, height = 1920)

# Distill volumes of classified RoD
dist <- create.budget.classes(delta_z_all, lod_crit, yres(delta_z_all)) %>% 
  filter(!is.na(cell_vol)) %>% 
  mutate(cell_status = as.factor(if_else(discarded == T, "Discarded", "Valid")))

# Plot and save elevation change distribution (ECD)
p_hist <- ggplot(dist, aes(values.raw_raster., fill = cell_status)) +
  geom_histogram(binwidth = 0.01) +
  labs(x = "Elevation change [m]", y = "Count", fill = "Cell status") + 
  theme_bw() +
  scale_fill_manual(values = c("#000000", "#35b779")) +
  scale_x_continuous(limits = c(-1, 1)) +
  # scale_y_continuous(limits = c(0, 500)) +
  theme(legend.position = c(0.15, 0.88), legend.text = element_text(size=15), legend.title = element_text(size=15)) +
  theme(axis.line = element_line(color='black'),
        panel.grid.minor = element_blank(),
        panel.border = element_blank())
p_hist
ggsave(output_lod_hist_path, plot = p_hist, height=1800, width=2200, units ="px")

# Plot and save elevation change distribution with limited y-scale (ECD)
p_hist500 <- ggplot(dist, aes(values.raw_raster., fill = cell_status)) +
  geom_histogram(binwidth = 0.01) +
  labs(x = "Elevation change [m]", y = "Count", fill = "Cell status") + 
  theme_bw() +
  # scale_fill_manual(values = c("#fee090", "#74add1")) +
  # scale_fill_manual(values = c("#35b779", "#31688e")) +
  scale_fill_manual(values = c("#000000", "#35b779")) +
  scale_x_continuous(limits = c(-1, 1)) +
  scale_y_continuous(limits = c(0, 500)) +
  theme(legend.position = c(0.15, 0.88), legend.text = element_text(size=15), legend.title = element_text(size=15)) +
  theme(axis.line = element_line(color='black'),
        panel.grid.minor = element_blank(),
        panel.border = element_blank())
p_hist500
ggsave(output_lod_hist500_path, plot = p_hist500, height=1800, width=2200, units ="px")
# Summarise distilled terrestrial volumes
dist_sum <- dist %>%
  filter(!is.na(values.raw_raster.)) %>% 
  group_by(class) %>%
  summarize(n=n(), vol=sum(cell_vol), area=n*raster_res^2) 
.rowNamesDF(dist_sum, make.names=FALSE) <- dist_sum$class
# Add variables to summary
interval <- paste(t0_cfg$survey_date_pret, t1_cfg$survey_date_pret, sep =" - ")
export_detailed <- dist_sum %>% 
  mutate(interval,
         flood_startdate,
         lod_crit_m = lod_crit,
         perspective = perspective,
         raster_res_m = raster_res,
         reported_date = timestamp,
         comment)
# Save detailed summary
write.table(export_detailed, file = "C:/Daten/math_gubelyve/pcc_standalone/export/budget_results_detailed.csv",
            append = T, sep = ";", row.names = F, col.names = F)
# Save shortened summary
export <- data.frame(interval) %>% 
  mutate(flood_startdate,
         lod_crit_m = lod_crit,
         perspective = perspective,
         Ero_validvol_m3 = dist_sum$vol[dist_sum$class=="Erosion"],
         Ero_discvol_m3 = dist_sum$vol[dist_sum$class=="Discarded Erosion"],
         Ero_totvol_m3 = Ero_validvol_m3 + Ero_discvol_m3,
         Ero_lossvol_percent = 100*Ero_discvol_m3/Ero_totvol_m3,
         Ero_validarea_m2 = dist_sum$area[dist_sum$class=="Erosion"],
         Ero_discarea_m2 = dist_sum$area[dist_sum$class=="Discarded Erosion"],
         Ero_totarea_m2 = Ero_validarea_m2 + Ero_discarea_m2,
         Ero_validzoneavg_m_per_cell = Ero_validvol_m3/Ero_validarea_m2,
         Depo_validvol_m3 = dist_sum$vol[dist_sum$class=="Deposition"],
         Depo_discvol_m3 = dist_sum$vol[dist_sum$class=="Discarded Deposition"],
         Depo_totvol_m3 = Depo_validvol_m3 + Depo_discvol_m3,
         Depo_lossvol_percent = 100*Depo_discvol_m3/Depo_totvol_m3,
         Depo_validarea_m2 = dist_sum$area[dist_sum$class=="Deposition"],
         Depo_discarea_m2 = dist_sum$area[dist_sum$class=="Discarded Deposition"],
         Depo_totarea_m2 = Depo_validarea_m2 + Depo_discarea_m2,
         Depo_validzoneavg_m_per_cell = Depo_validvol_m3/Depo_validarea_m2,
         raster_res_m = raster_res,
         reported_date = timestamp,
         comment,
         min_delta_z = minValue(delta_z_all),
         max_delta_z = maxValue(delta_z_all))
write.table(export, file = "C:/Daten/math_gubelyve/pcc_standalone/export/budget_results.csv",
            append = T, sep = ";", row.names = F, col.names = F)

# Additional barplot of volume distribution of calculated budget
budget_title <- paste(t1_cfg$survey_date_pret, "-", t0_cfg$survey_date_pret, sep = " ")
p_bud <- ggplot(dist_sum, aes(fill=class, x=1, y=vol)) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_viridis(discrete = T) +
  ggtitle(budget_title) +
  theme_minimal() +
  theme(axis.line = element_line(color='black'),
        panel.grid.minor = element_blank(),
        panel.border = element_blank()) +
  ylab("Volume [m²]") +
  xlab("")
ggsave(output_lod_bar_path, plot = p_bud, height=1800, width=2200, units ="px")
# Save RoD with correspondent ECD as cowplot
p_hist_grob <- p_hist + theme(legend.position = c(0.17, 0.85), legend.text = element_text(size=12), legend.title = element_text(size=12))
p_elev_unvert_grob <- tm_elev_uncert +
  tm_layout(frame = F, 
            legend.title.size = 1.1, legend.text.size = 0.8, 
            legend.outside = F, legend.position = c("left", "center"),
            main.title.position = "center", main.title.size = 1.3)

p_elev_uncert <- tmap_grob(p_elev_unvert_grob)
plot_grid(p_elev_uncert, p_hist_grob, nrow = 1, labels = c('1', '2'), label_size = 12)
ggsave(output_budget_path, height=1400, width=2800, units ="px")

p_hist500_grob <- p_hist500 + theme(legend.position = c(0.17, 0.85), legend.text = element_text(size=12), legend.title = element_text(size=12))
plot_grid(p_elev_uncert, p_hist500_grob, nrow = 1, labels = c('1', '2'), label_size = 12)
ggsave(output_budget500_path, height=1400, width=2800, units ="px")
