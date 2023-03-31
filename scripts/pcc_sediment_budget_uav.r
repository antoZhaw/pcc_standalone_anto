# libraries ---------------------------------------------------------------
# install packages

library(lidR) # Point cloud classification
library(papeR) # summary tables
library(tidyverse) # tidy essentials (ggplot, purr, tidyr, readr, dplyr)
library(lubridate) # handling dates and time
library(tmap) # map visualization
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
library(sabre) # for mapcurves (goodness-of-fit function)
library(fasterize) # for faster rasterization
library(psych) # for cohen's kappa

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

classify.gnd <- function(las, steep_slopes = F, class_thres, cloth_res, rigid) {
  mycsf <- csf(steep_slopes, class_thres, cloth_res, rigid)
  las <- classify_ground(las, mycsf)
  las_gnd <- filter_poi(las, Classification == LASGROUND)
  # plot(las_gnd, size = 1, color = "RGB", bg = "white")
  return(las_gnd)
}

set.RGLtopview <- function(x_scale = 800, y_scale = 800) {
  view3d(theta = 0, phi = 0, zoom = 0.8)
  par3d(windowRect = c(30, 30, x_scale, y_scale))
}

mask.raster.layer <- function(raster_layer) {
  raster_layer[raster_layer != 0] <- 0
  raster_layer[is.na(raster_layer)] <- 1
  raster_layer
}

normalise.raster <- function(raster_layer, normal_value = 1) {
  raster_layer[raster_layer != 0] <- normal_value
  raster_layer[is.na(raster_layer)] <- 0
  raster_layer
}

filter.raster <- function(raster_layer, filter_value = 1) {
  raster_layer[raster_layer != filter_value] <- 0
  raster_layer[raster_layer == filter_value] <- 1
  raster_layer[is.na(raster_layer)] <- 0
  raster_layer
}

value.to.na.raster <- function(raster_layer, na_value = 0) {
  raster_layer[raster_layer == na_value] <- NA
  raster_layer
}

discard.uncertain.raster <- function(raster_layer, z_level_of_detection) {
  raster_layer[abs(raster_layer) <= z_level_of_detection] <- NA
  raster_layer
}

plot.csf.result.vs.target <- function(raster_bin, target_shp, plot_title) {
  # oranges <- tmaptools::get_brewer_pal("Oranges", n = 2, contrast = c(0.3, 0.9))
  tmap_mode("plot") + # "plot" or "view"
  tm_shape(raster_bin) +
  tm_raster(title = plot_title, alpha = 1) +
  # tm_raster(palette = oranges, title = plot_title, alpha = 1) +
  tm_shape(target_shp) +
  tm_polygons(alpha = 0.5) +
  # tm_shape(wallows_raster_100) +
  # tm_raster(palette = purples, title = "Wallows", alpha = 1) +
  tm_view(control.position = c("right", "top"))
}

# Globals for Configuration-----------------------------------------------------
# Record start date and time
start <- as_datetime(lubridate::now())
date <- as.Date(start)
hour <- hour(start)
minute <- minute(start)
# Generate timestamp for this run.
timestamp <- as.character(paste(date, hour, minute, sep = "-"))

# Settings which apply for t0 and t1.
wholeset <- T
perspective <- "uav"
settype <- if_else(wholeset == T, "wholeset", "subset")
raster_res <- 0.4


# Load environment dependent paths.
user <- Sys.getenv("USERNAME")
if(user == "gubelyve"| user == "xgby"){
  dir_repo <- "C:/Daten/math_gubelyve/pcc_standalone"
  dir_data <- "C:/Daten/math_gubelyve"
} else{
  dir_repo <- "C:/code_wc/pcc_standalone"
  dir_data <- "C:/code_wc/math_gubelyve"
}
dir_persp <- if_else(perspective == "tls", "tls_data", "uav_data")

# Settings t0
t0_dataset_id <- "2"
t0_year <- "2020"
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


# Settings t1
t1_dataset_id <- "1"
t1_year <- "2020"
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

aoi_path <- file.path(dir_repo, "data/area_of_interest_final", "AOI_final.shp", fsep="/")
aoi_dir <- file.path(dir_repo, "data/area_of_interest_final", fsep="/")
bud_output_path <- file.path(dir_data, dir_persp, "budget", timestamp, fsep="/")

# Create run specific output folder
dir.create(bud_output_path)
dir.create(t0_output_path)
dir.create(t1_output_path)

# 
# output_json_name <- as.character(paste(output_id, ".json", sep = ""))
# output_json_path <- file.path(output_path, output_json_name, fsep="/")
# 
# output_target_name <- as.character(paste(output_id, ".png", sep = ""))
# output_target_path <- file.path(output_path, output_target_name, fsep="/")
# 
# output_csv_name <- as.character(paste(output_id, ".csv", sep = ""))
# output_csv_path <- file.path(output_path, output_csv_name, fsep="/")
# 
# output_las_inp_name <- as.character(paste(output_id, "-inp.txt", sep = ""))
# output_las_inp_path <- file.path(output_path, output_las_inp_name, fsep="/")
# 
# output_las_aoi_name <- as.character(paste(output_id, "-aoi.txt", sep = ""))
# output_las_aoi_path <- file.path(output_path, output_las_aoi_name, fsep="/")
# 
# output_las_class_name <- as.character(paste(output_id, "-classified.txt", sep = ""))
# output_las_class_path <- file.path(output_path, output_las_class_name, fsep="/")
# 
# output_asc_name <- as.character(paste(output_id, ".asc", sep = ""))
# output_asc_path <- file.path(output_path, output_asc_name, fsep="/")
# 
# output_ncdf_name <- as.character(paste(output_id, ".nc", sep = ""))
# output_ncdf_path <- file.path(output_path, output_ncdf_name, fsep="/")
# 
# output_las_sed_name <- as.character(paste(output_id, "-sed.las", sep = ""))
# output_las_sed_path <- file.path(output_path, output_las_sed_name, fsep="/")
# 
# output_las_gnd_name <- as.character(paste(output_id, "-gnd.las", sep = ""))
# output_las_gnd_path <- file.path(output_path, output_las_gnd_name, fsep="/")
# 
# output_las_all_name <- as.character(paste(output_id, "-all.las", sep = ""))
# output_las_all_path <- file.path(output_path, output_las_all_name, fsep="/")
# 
# output_target_rast_name <- as.character(paste(output_id, "-target.png", sep = ""))
# output_target_rast_path <- file.path(output_path, output_target_rast_name, fsep="/")

# Read files--------------------------------------------------------------------

# empty warnings if existing.
if(length(warnings())!=0){
  assign("last.warning", NULL, envir = baseenv())
}

# Generate rectangular polygon of area of interest
gen_xy <- structure(list(dat = c("BURR-1-1-1", "BURR-1-1-1", 
                                 "BURR-1-1-1", "BURR-1-1-1"),
                         Longitude = c(2575011, 2575011, 
                                      2575517.4, 2575517.4),
                         Latitude = c(1178366.7, 1178993, 
                                       1178993, 1178366.7)),
                    class = "data.frame", row.names = c(NA,-4L))

bounding_box <- gen_xy %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 2056) %>%
  group_by(dat) %>%
  summarise(geometry = st_combine(geometry)) %>%
  st_cast("POLYGON")

# Read Shapefiles
# mctar_all <- read_sf(dsn = mctar_path) 
t0_csf_aoi_shp <- read_sf(dsn = t0_mcdut_path)
t1_csf_aoi_shp <- read_sf(dsn = t1_mcdut_path)
t0_mctar_all <- read_sf(dsn = t0_mctar_path)
t1_mctar_all <- read_sf(dsn = t1_mctar_path)

aoi_shp <- read_sf(dsn = aoi_path)
# Intersect target with area of interest
t0_targets_aoi_shp <- st_intersection(t0_mctar_all, bounding_box)
t1_targets_aoi_shp <- st_intersection(t1_mctar_all, bounding_box)


t0_target_wat <- t0_targets_aoi_shp %>%  filter(Id == 1)
t0_target_sed <- t0_targets_aoi_shp %>%  filter(Id == 2)

t1_target_wat <- t1_targets_aoi_shp %>%  filter(Id == 1)
t1_target_sed <- t1_targets_aoi_shp %>%  filter(Id == 2)

# Intersect target with area of interest
# mctar_bb <- st_intersection(mctar_all, bounding_box)

# Separate targets filtered by class
# mctar_water <- mctar_bb %>%  filter(Id == 1)
# mctar_sed <- mctar_bb %>%  filter(Id == 2)
# mctar_veg <- mctar_bb %>%  filter(Id == 3)
# mctar_rock <- mctar_bb %>%  filter(Id == 5)

# targeted_class <- mctar_sed

# Generate target beforehand with raster----------------------------------------
# Calculate number of rows and columns of target, based on extent
# This calculation results in a resolution of 0.5
# tar_ncol <- 2*(xmax(extent(targeted_class)) - xmin(extent(targeted_class)))
# tar_nrow <- 2*(ymax(extent(targeted_class)) - ymin(extent(targeted_class)))
# Generate target. This is done beforehand since it is required only once.
# tar_raw <- raster(as(targeted_class, "Spatial"), ncols = tar_ncol, nrows = tar_nrow)
# target <- rasterize(as(targeted_class, "Spatial"), tar_raw, getCover = TRUE, progress = "text")

# Information about rasterized target
# target$layer
# st_crs(target)

# Plot features on one plot
# ggplot() + 
#   geom_sf(data = AOI_xy, mapping = aes()) +
#   geom_sf(data = mctar_sed, mapping = aes()) +
#   coord_sf(crs = st_crs(2056))

# test <- st_intersection(mctar_sed, AOI_xy)

# ggplot() +
#   geom_sf(data = mcdut, mapping = aes()) +
#   coord_sf(crs = st_crs(2056))

# Read LAS
# Intensity (i), color information (RGB), number of Returns (r), classification (c)
# of the first point is loaded only to reduce computational time.
t0_las <- readLAS(t0_data_path, select = "xyzRGBc", filter = t0_cfg$las_filter)
t1_las <- readLAS(t1_data_path, select = "xyzRGBc", filter = t1_cfg$las_filter)

# Filter points which are not within area of interest---------------------------
# las <- classify_poi(las, class = LASNOISE, roi = aoi_shp, inverse_roi = T)
# las <- filter_poi(las, Classification != LASNOISE)
# plot(las, size = 1, color = "RGB", bg = "white", axis = F)

# Reset class LASNOISE for further procedure
t0_las$Classification <- LASNONCLASSIFIED
t1_las$Classification <- LASNONCLASSIFIED

# Create copy of read LAS to omit loading procedure.
# las_origin <- las
# las <- las_origin

# Segment Water with Cloth Simulation Filter------------------------------------
t0_rigid_n_water <- t0_cfg$csf_water_rigidness
t1_rigid_n_water <- t1_cfg$csf_water_rigidness
t0_class_thres_water <- t0_cfg$csf_water_class_threshold
t1_class_thres_water <- t1_cfg$csf_water_class_threshold
t0_cloth_res_water <- t0_cfg$csf_water_cloth_resolution
t1_cloth_res_water <- t1_cfg$csf_water_cloth_resolution
t0_steep_slope_water <- if_else(t0_rigid_n_water == 3, F, T)
t1_steep_slope_water <- if_else(t1_rigid_n_water == 3, F, T)

t0_status_water <- as.character(paste("00_water", t0_year, "_rigid", t0_rigid_n_water,
                                   "_clthres", t0_class_thres_water,
                                   "clothres", t0_cloth_res_water, sep = ""))
t1_status_water <- as.character(paste("00_water", t1_year, "_rigid", t1_rigid_n_water,
                                      "_clthres", t1_class_thres_water,
                                      "clothres", t1_cloth_res_water, sep = ""))

t0_las_water <- classify.gnd(t0_las, t0_steep_slope_water, t0_class_thres_water, t0_cloth_res_water, t0_rigid_n_water)
t1_las_water <- classify.gnd(t1_las, t1_steep_slope_water, t1_class_thres_water, t1_cloth_res_water, t1_rigid_n_water)

t0_las_water <- add_attribute(t0_las_water, FALSE, "water")
t1_las_water <- add_attribute(t1_las_water, FALSE, "water")
t0_las_water$water <- if_else(t0_las_water$Classification == LASGROUND, T, F)
t1_las_water$water <- if_else(t1_las_water$Classification == LASGROUND, T, F)

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

# Rasterize water point cloud    
t0_DEM_water <- rasterize_canopy(t0_las_water, res = raster_res, p2r(), pkg = "raster")
t1_DEM_water <- rasterize_canopy(t1_las_water, res = raster_res, p2r(), pkg = "raster")

# Determine convex hull of both DEM
# tbd: unclear, whether this has an effect on alignment..
# raster_ext <- extent(min(c(xmin(t0_DEM_water), xmin(t1_DEM_water))),
#                     max(c(xmax(t0_DEM_water), xmax(t1_DEM_water))),
#                     min(c(ymin(t0_DEM_water), ymin(t1_DEM_water))),
#                     max(c(ymax(t0_DEM_water), ymax(t1_DEM_water))))
# raster_ext

# t0_raster_water <- raster(nrows=nrow(t0_DEM_water), ncols=ncols(t0_DEM_water), crs=2056,
#                   ext=raster_ext, resolution=raster_res, vals=NULL)
# t1_raster_water <- raster(nrows=nrow(t1_DEM_water), ncols=ncols(t1_DEM_water), crs=2056,
#                        ext=raster_ext, resolution=raster_res, vals=NULL)

# plot(t0_DEM_water)
# plot(t1_DEM_water)

t0_mask_water <- mask.raster.layer(t0_DEM_water)
t1_mask_water <- mask.raster.layer(t1_DEM_water)

# plot(t0_mask_water)
# plot(t1_mask_water)


# Segment Ground with Cloth Simulation Filter-----------------------------------
# classify ground
t0_rigid_n_sed <- t0_cfg$csf_gnd_rigidness
t1_rigid_n_sed <- t1_cfg$csf_gnd_rigidness
t0_class_thres_sed <- t0_cfg$csf_gnd_class_threshold
t1_class_thres_sed <- t1_cfg$csf_gnd_class_threshold
t0_cloth_res_sed <- t0_cfg$csf_gnd_cloth_resolution
t1_cloth_res_sed <- t1_cfg$csf_gnd_cloth_resolution
t0_steep_slope_sed <- if_else(t0_rigid_n_sed == 3, F, T)
t1_steep_slope_sed <- if_else(t1_rigid_n_sed == 3, F, T)

t0_las_sed <- classify.gnd(t0_las, t0_steep_slope_sed, t0_class_thres_sed, t0_cloth_res_sed, t0_rigid_n_sed)
t0_las_sed <- classify_poi(t0_las_sed, class = LASNOISE, roi = t0_csf_aoi_shp, inverse_roi = T)
t0_las_sed <- filter_poi(t0_las_sed, Classification != LASNOISE)

# classify ground
t1_las_sed <- classify.gnd(t1_las, t1_steep_slope_sed, t1_class_thres_sed, t1_cloth_res_sed, t1_rigid_n_sed)
t1_las_sed <- classify_poi(t1_las_sed, class = LASNOISE, roi = t0_csf_aoi_shp, inverse_roi = T)
t1_las_sed <- filter_poi(t1_las_sed, Classification != LASNOISE)

# Rasterize sediment
t0_DEM_sed <- rasterize_canopy(t0_las_sed, res = raster_res, p2r(), pkg = "raster")
t1_DEM_sed <- rasterize_canopy(t1_las_sed, res = raster_res, p2r(), pkg = "raster")

# Subtract water
t0_sed <- t0_DEM_sed * t0_mask_water
t1_sed <- t1_DEM_sed * t1_mask_water

delta_sed <- t0_sed - t1_sed

col <- height.colors(30)
par(mfrow=c(1,3))

# plot(t0_sed, col = col, main = "UAV t0")
# plot(t1_sed, col = col, main = "UAV t1")
# plot(delta_sed, col = col, main = "DEM of difference")


t0_tm_sed <- normalise.raster(t0_sed)
t1_tm_sed <- normalise.raster(t1_sed)
t0_tm_wat <- normalise.raster(t0_DEM_water)
t1_tm_wat <- normalise.raster(t1_DEM_water)

# Plot comparison between target and classified raster
t0_tm_sed_result <- plot.csf.result.vs.target(t0_tm_sed, t0_target_sed, "GOF: Sediment t0")
t0_tm_sed_result

t1_tm_sed_result <- plot.csf.result.vs.target(t1_tm_sed, t1_target_sed, "GOF: Sediment t1")
t1_tm_sed_result

t0_tm_wat_result <- plot.csf.result.vs.target(t0_tm_wat, t0_target_wat, "GOF: Water t0")
t0_tm_wat_result

t1_tm_wat_result <- plot.csf.result.vs.target(t1_tm_wat, t1_target_wat, "GOF: Water t1")
t1_tm_wat_result

# Determine habitate change
t0_tm_hab <- normalise.raster(t0_sed, 10)
t1_tm_hab <- t1_tm_sed
# Generate raster of pseudo factors (with values 0, 1, 10, 11)
tm_habitate <- t0_tm_hab + t1_tm_hab


pal4div <- c("#FFFFFF", "#4daf4a", "#e41a1c", "#377eb8")
tmap_mode("plot") + # "plot" or "view"
  tm_shape(tm_habitate) +
  # tm_raster(title = "2D Habitat change of Sediment", style = "cat",  ) +
  tm_raster(title = "2D Habitat change of Sediment", palette = pal4div, alpha = 1, style = "cat", breaks = c(0, 2, 10.5),
            labels = c("no change", "deposition", "erosion", "change in elevation") ) +
  # tm_raster(palette = oranges, title = plot_title, alpha = 1) +
  # tm_shape(target_shp) +
  # tm_polygons(alpha = 0.5) +
  # tm_shape(wallows_raster_100) +
  # tm_raster(palette = purples, title = "Wallows", alpha = 1) +
  tm_view(control.position = c("right", "top"))

# Generate mask for cells which show a change in elevation (pick value 11)
tm_elevation_mask <- filter.raster(tm_habitate, 11)
# Set zero values to na 
tm_elevation_mask <- value.to.na.raster(tm_elevation_mask)
# Apply mask on delta z
delta_z_all <- tm_elevation_mask*(t1_sed - t0_sed)

# Plot elevation change without uncertainty assessment
tmap_mode("plot") + # "plot" or "view"
  tm_shape(delta_z_all) +
  tm_raster(title = "Elevation change of Sediment, raw", alpha = 1, style = "cont")

# Calculate critical level of detection
lod_crit <- sqrt(t1_cfg$z_mean_accuracy^2 + t0_cfg$z_mean_accuracy^2)

delta_z_cleaned <- discard.uncertain.raster(delta_z_all, lod_crit)

# Plot elevation change with uncertainty assessment
tmap_mode("plot") + # "plot" or "view"
  tm_shape(delta_z_cleaned) +
  tm_raster(title = "Elevation change of Sediment, clean", alpha = 1, style = "cont")

dist <- data.frame(values(delta_z_all)) %>% 
  mutate(discarded = if_else(values.delta_z_all. <= lod_crit, T, F))

hist(dist$values.delta_z_all.)

summary(dist)

# plot(las, size = 1, color = "RGB", bg = "white")


# good values for sediment (rigidness=1)
# class_thres_i <- c(0.9, 0.85, 0.8, 0.75)
# cloth_res_i <- c(1.8, 1.7, 1.6, 1.5)
# 
# class_thres_i <- c(0.4, 0.85)
# cloth_res_i <- c(5.0, 1.7)

# good values for sediment (rigidness=1)
# cloth_res_i <- seq(from = 2.5, to = 7.0, by = 0.1)
# class_thres_i <- seq(from = 0.14, to = 0.8, by = 0.01)
cloth_res_i <- seq(from = 1.9, to = 6.5, by = 0.2)
class_thres_i <- seq(from = 0.3, to = 0.9, by = 0.03)

col <- height.colors(15)
class_id <- targeted_class$Id

df <- data.frame(name=c(""), class=c(""), rigidness=c(""), classthreshold=c(""),
                  clothresolution=c(""), steepslope=c(""), n_obs=c(""), 
                 kappa=c(""), comp_time_sec=c(""), rasterresolution=c(""))

    start_ij <- as_datetime(lubridate::now())
    rigid_n <- 1
    status_sed <- as.character(paste("RGL", n, "_rigid", rigid_n, "_clthres", i, "clothres", j, sep = ""))
    print(status_sed)
    las_ij <- classify.gnd(las, i, j, rigid_n)
    las_ij <- add_attribute(las_ij, FALSE, "ground")
    las_ij$ground <- if_else(las_ij$Classification == LASGROUND, T, F)

    # las$Classification <- LASNONCLASSIFIED

    las_ij <- filter_poi(las_ij, Classification == LASGROUND)
    # plot(las_ij, size = 1, color = "RGB", bg = "white", axis = F)
    
    # Filter points which are not within area of interest
    las_ij <- classify_poi(las_ij, class = LASNOISE, roi = mcdut, inverse_roi = T)
    las_ij <- filter_poi(las_ij, Classification != LASNOISE)
    plot(las_ij, size = 1, color = "RGB", bg = "black", axis = F)
    set.RGLtopview()
    output_png_name <- as.character(paste(status_sed, ".png", sep = ""))
    output_png_path <- file.path(output_path, output_png_name, fsep="/")
    rgl.snapshot(output_png_path)
    rgl.close()

    # Rasterize point cloud    
    DEM_ij <- rasterize_canopy(las_ij, res = raster_res, p2r(), pkg = "raster")

    # Subtract water course
    sed_ij <- DEM_ij * mask_water
    
    # Save plot of raster
    output_sed_rast_name <- as.character(paste(status_sed, ".png", sep = ""))
    output_sed_rast_path <- file.path(output_path, output_sed_rast_name, fsep="/")
    png(output_sed_rast_path, height=nrow(sed_ij), width=ncol(sed_ij)) 
    plot(sed_ij, maxpixels=ncell(sed_ij), legend =F)
    dev.off()
    
    # Not used here, since static water raster is used.
    # raster_ext <- extent(xmin(DEM_ij), xmax(DEM_ij), ymin(DEM_ij), ymax(DEM_ij))
    tar_raw <- raster(nrows=nrow(sed_ij), ncols=ncols(sed_ij), crs=2056,
                      ext=raster_ext, resolution=raster_res, vals=NULL)
    target <- fasterize(targeted_class, tar_raw, field = "Id", fun="sum")
    
    # Normalise raster features
    rater1 <- values(target)
    rater1[rater1 != 0] <- 1
    rater1[is.na(rater1)] <- 0
    rater2 <- values(sed_ij)
    rater2[rater2 != 0] <- 1
    rater2[is.na(rater2)] <- 0
    
    kap <- cohen.kappa(x=cbind(rater1,rater2))

    end_ij <- as_datetime(lubridate::now())
    timespan_ij <- interval(start_ij, end_ij)
    delta_t <- as.numeric(timespan_ij, "seconds")
    obs <- c(status_sed, class_id, rigid_n, i, j, "FALSE", kap$n.obs, kap$kappa, delta_t, raster_res)
    df <- rbind(df, obs)

write_delim(df, file=output_csv_path, delim = ";")

png(output_target_rast_path, height=nrow(target), width=ncol(target)) 
plot(target, maxpixels=ncell(target), legend =F)
dev.off()


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
