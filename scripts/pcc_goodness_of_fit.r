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

classify.gnd <- function(las, class_thres, cloth_res, rigid) {
  mycsf <- csf(F, class_thres, cloth_res, rigid)
  las <- classify_ground(las, mycsf)
  las_gnd <- filter_poi(las, Classification == LASGROUND)
  # plot(las_gnd, size = 1, color = "RGB", bg = "white")
  return(las_gnd)
}

set.RGLtopview <- function(x_scale = 800, y_scale = 800) {
  view3d(theta = 0, phi = 0, zoom = 0.8)
  par3d(windowRect = c(30, 30, x_scale, y_scale))
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
  dir_repo <- "C:/code_wc/pcc_standalone"
  dir_data <- "C:/code_wc/math_gubelyve"
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

mctar_shp_name <- as.character(paste(cfg$mapcurve_target_shp, sep = ""))
mctar_path <- file.path(dir_repo, mctar_shp_name, fsep="/")

mcdut_shp_name <- as.character(paste(cfg$mapcurve_dut_shp, sep = ""))
mcdut_path <- file.path(dir_repo, mcdut_shp_name, fsep="/")

output_json_name <- as.character(paste(output_id, ".json", sep = ""))
output_json_path <- file.path(output_path, output_json_name, fsep="/")

output_target_name <- as.character(paste(output_id, ".png", sep = ""))
output_target_path <- file.path(output_path, output_target_name, fsep="/")

output_csv_name <- as.character(paste(output_id, ".csv", sep = ""))
output_csv_path <- file.path(output_path, output_csv_name, fsep="/")

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
mctar_all <- read_sf(dsn = mctar_path) 
mcdut <- read_sf(dsn = mcdut_path) 
aoi_shp <- read_sf(dsn = aoi_path)

# Intersect target with area of interest
mctar_bb <- st_intersection(mctar_all, bounding_box)

# Separate targets filtered by class
mctar_water <- mctar_bb %>%  filter(Id == 1)
mctar_sed <- mctar_bb %>%  filter(Id == 2)
mctar_veg <- mctar_bb %>%  filter(Id == 3)
mctar_rock <- mctar_bb %>%  filter(Id == 5)

targeted_class <- mctar_water

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
las <- readLAS(data_path, select = "xyzRGBc", filter = cfg$las_filter)

# Filter points which are not within area of interest---------------------------
# las <- classify_poi(las, class = LASNOISE, roi = aoi_shp, inverse_roi = T)
# las <- filter_poi(las, Classification != LASNOISE)
# plot(las, size = 1, color = "RGB", bg = "white", axis = F)

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

# Segment Ground with Cloth Simulation Filter-----------------------------------
# setup csf filter settings
# rigidness: does not seem to have much impact.
# class_threshold and cloth_resolution influence each other. 0.5 x 0.5 is more conservative compared to 0.5 x 1.
# If sharper watercourse desired: increase threshold to 0.7 (leads to more canopy in ground points)
# plot(las, size = 1, color = "RGB", bg = "white")

las_origin <- las

las <- las_origin
par(mfrow=c(1,1))

# good values for sediment (rigidness=1)
# class_thres_i <- c(0.9, 0.85, 0.8, 0.75)
# cloth_res_i <- c(1.8, 1.7, 1.6, 1.5)
# 
# class_thres_i <- c(0.22)
# cloth_res_i <- c(3.5)

# good values for sediment (rigidness=1)
# class_thres_i <- c(0.19, 0.2, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.3, 0.31)
# cloth_res_i <- c(3.6, 3.5, 3.4, 3.3, 3.2, 3.1, 3.0)
cloth_res_i <- seq(from = 3.5, to = 7.0, by = 0.1)
class_thres_i <- seq(from = 0.3, to = 0.8, by = 0.01)

col <- height.colors(15)
raster_res <- 0.5
class_id <- targeted_class$Id
# 
# record <- data.frame(names(c("class", "rigidness", "class threshold",
#                              "cloth resolution", "steep slope", "extent", "n_obs",
#                              "kappa_low","kappa_upper", "kappa")))


df <- data.frame(name=c(""), class=c(""), rigidness=c(""), classthreshold=c(""),
                  clothresolution=c(""), steepslope=c(""), n_obs=c(""), 
                 kappa=c(""), comp_time_sec=c(""), rasterresolution=c(""))

n <- 1L
for (i in class_thres_i) {
  for (j in cloth_res_i) {
    start_ij <- as_datetime(lubridate::now())
    rigid_n <- 3
    status <- as.character(paste("RGL", n, "_rigid", rigid_n, "_clthres", i, "clothres", j, sep = ""))
    print(status)
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
    output_png_name <- as.character(paste(status, ".png", sep = ""))
    output_png_path <- file.path(output_path, output_png_name, fsep="/")
    rgl.snapshot(output_png_path)
    rgl.close()

    # Rasterize point cloud    
    DEM_ij <- rasterize_canopy(las_ij, res = raster_res, p2r(), pkg = "raster")

    raster_ext <- extent(xmin(DEM_ij), xmax(DEM_ij), ymin(DEM_ij), ymax(DEM_ij))
    tar_raw <- raster(nrows=nrow(DEM_ij), ncols=ncols(DEM_ij), crs=2056,
                      ext=raster_ext, resolution=raster_res, vals=NULL)
    target <- fasterize(targeted_class, tar_raw, field = "Id", fun="sum")
    
    # Normalise raster features
    rater1 <- values(target)
    rater1[rater1 != 0] <- 1
    rater1[is.na(rater1)] <- 0
    rater2 <- values(DEM_ij)
    rater2[rater2 != 0] <- 1
    rater2[is.na(rater2)] <- 0
    
    kap <- cohen.kappa(x=cbind(rater1,rater2))

    end_ij <- as_datetime(lubridate::now())
    timespan_ij <- interval(start_ij, end_ij)
    delta_t <- as.numeric(timespan_ij, "seconds")
    obs <- c(status, class_id, rigid_n, i, j, "FALSE", kap$n.obs, kap$kappa, delta_t, raster_res)
    df <- rbind(df, obs)
    n <- n + 1
  }
}
write_delim(df, file=output_csv_path, delim = ";")

plot(target, legend =F)

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
