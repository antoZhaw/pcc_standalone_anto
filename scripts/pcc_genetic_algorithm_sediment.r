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
library(GA) # genetic algorithm for parameter search optimisation


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

cohen.kappa.csf <- function(raw_las, ga_aoi_shp, targets_shp, ga_output_path,
                            global_rigid, ct_sed_n, clr_sed_o, 
                            ct_wat_i, clr_wat_j, raster_res) {
  id_ij <- sample(1:999, 1)
  start_ij <- as_datetime(lubridate::now())
  rig_sed_m <- global_rigid
  rig_wat_h <- global_rigid
  # classify ground
  msg_sed <- as.character(paste("SED", id_ij, "_rig", round(rig_sed_m, 4), "_ct", round(ct_sed_n, 4), "_clr", round(clr_sed_o, 4), sep = ""))
  print(msg_sed)
  las_sed_ij <- classify.gnd(las, ct_sed_n, clr_sed_o, rig_sed_m)
  las_sed_ij <- classify_poi(las_sed_ij, class = LASNOISE, roi = ga_aoi_shp, inverse_roi = T)
  las_sed_ij <- filter_poi(las_sed_ij, Classification != LASNOISE)
  # classify water surface
  msg_wat <- as.character(paste("WAT", id_ij,"_rig", round(rig_wat_h, 4), "_ct", round(ct_wat_i, 4), "_clr", round(clr_wat_j, 4), sep = ""))
  print(msg_wat)
  las_wat_ij <- classify.gnd(las, ct_wat_i, clr_wat_j, rig_wat_h)
  las_wat_ij <- classify_poi(las_wat_ij, class = LASNOISE, roi = ga_aoi_shp, inverse_roi = T)
  las_wat_ij <- filter_poi(las_wat_ij, Classification != LASNOISE)
  # Save plot of classified water surface
  # plot(las_wat_ij, size = 1, color = "RGB", bg = "black", axis = F)
  # set.RGLtopview()
  # output_wat_png_name <- as.character(paste("LAS", msg_wat, ".png", sep = ""))
  # output_wat_png_path <- file.path(ga_output_path, output_wat_png_name, fsep="/")
  # rgl.snapshot(output_wat_png_path)
  # rgl.close()
  # Rasterize both point clouds    
  DEM_sed_ij <- rasterize_canopy(las_sed_ij, res = raster_res, p2r(), pkg = "raster")
  DEM_wat_ij <- rasterize_canopy(las_wat_ij, res = raster_res, p2r(), pkg = "raster")
  
  # Save plot of water raster. Currently disabled since las is already saved.
  # output_wat_rast_name <- as.character(paste("RAS", msg_wat, ".png", sep = ""))
  # output_wat_rast_path <- file.path(ga_output_path, output_wat_rast_name, fsep="/")
  # png(output_wat_rast_path, height=nrow(DEM_wat_ij), width=ncol(DEM_wat_ij))
  # plot(DEM_wat_ij, maxpixels=ncell(DEM_wat_ij), legend =F)
  # dev.off()
  
  # Generate water mask
  msk_wat_ij <- DEM_wat_ij
  msk_wat_ij[msk_wat_ij != 0] <- 0
  msk_wat_ij[is.na(msk_wat_ij)] <- 1            
  # Subtract water course
  sed_ij <- DEM_sed_ij * msk_wat_ij
  
  # Save plot of masked raster
  output_sed_rast_name <- as.character(paste(msg_sed, ".png", sep = ""))
  output_sed_rast_path <- file.path(ga_output_path, output_sed_rast_name, fsep="/")
  png(output_sed_rast_path, height=nrow(sed_ij), width=ncol(sed_ij))
  plot(sed_ij, maxpixels=ncell(sed_ij), legend =F)
  dev.off()
  
  # Generate raster for water comparison
  target_wat <- targets_shp %>%  filter(Id == 1)
  raster_ext_wat <- extent(xmin(DEM_wat_ij), xmax(DEM_wat_ij), ymin(DEM_wat_ij), ymax(DEM_wat_ij))
  tar_raw_wat <- raster(nrows=nrow(DEM_wat_ij), ncols=ncols(DEM_wat_ij), crs=2056,
                        ext=raster_ext_wat, resolution=raster_res, vals=NULL)
  tar_wat <- fasterize(target_wat, tar_raw_wat, field = "Id", fun="sum")

  # Generate raster for sediment comparison
  target_sed <- targets_shp %>%  filter(Id == 2)
  raster_ext_sed <- extent(xmin(sed_ij), xmax(sed_ij), ymin(sed_ij), ymax(sed_ij))
  tar_raw_sed <- raster(nrows=nrow(sed_ij), ncols=ncols(sed_ij), crs=2056,
                        ext=raster_ext_sed, resolution=raster_res, vals=NULL)
  tar_sed <- fasterize(target_sed, tar_raw_sed, field = "Id", fun="sum")
  
  # Normalise raster values for comparison
  rater1 <- values(tar_sed)
  rater1[rater1 != 0] <- 1
  rater1[is.na(rater1)] <- 0
  rater2 <- values(sed_ij)
  rater2[rater2 != 0] <- 1
  rater2[is.na(rater2)] <- 0
  kap_sed <- cohen.kappa(x=cbind(rater1,rater2))
  
  # Normalise raster values for comparison
  rater3 <- values(tar_wat)
  rater3[rater3 != 0] <- 1
  rater3[is.na(rater3)] <- 0
  rater4 <- values(DEM_wat_ij)
  rater4[rater4 != 0] <- 1
  rater4[is.na(rater4)] <- 0
  kap_wat <- cohen.kappa(x=cbind(rater3,rater4))
  
  end_ij <- as_datetime(lubridate::now())
  timespan_ij <- interval(start_ij, end_ij)
  delta_t <- as.numeric(timespan_ij, "seconds")
  total_kappa <- kap_sed$kappa
  iter_wat_msg <- as.character(paste("Total Kappa (wat): ", round(kap_wat$kappa, 4)))
  iter_result_msg <- as.character(paste("Total Kappa (sed): ", round(total_kappa, 4), ", computed in ", round(delta_t, 3), " seconds."))
  print(iter_wat_msg)
  print(iter_result_msg)
  obs <- as.character(paste(msg_sed, rig_sed_m, ct_sed_n, clr_sed_o, "FALSE", kap_sed$kappa, msg_wat, rig_wat_h, ct_wat_i, clr_wat_j, "FALSE", kap_wat$kappa, kap_sed$n.obs, delta_t, raster_res, sep =";"))
  output_csv_name <- as.character(paste("genetic_algo_report.csv", sep = ""))
  output_csv_path <- file.path(ga_output_path, output_csv_name, fsep="/")
  write(obs, file=output_csv_path, append = T)
  popsize_kappa <- total_kappa * 10000
  return(popsize_kappa)
}

# Globals for Configuration-----------------------------------------------------
# Specify dataset
dataset_id <- "1"
wholeset <- T
year <- "2021"
perspective <- "uav"
settype <- if_else(wholeset == T, "wholeset", "subset")
raster_res <- 0.5

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

ga_aoi_shp_name <- as.character(paste(cfg$mapcurve_dut_shp, sep = ""))
ga_aoi_path <- file.path(dir_repo, ga_aoi_shp_name, fsep="/")

output_json_name <- as.character(paste(output_id, ".json", sep = ""))
output_json_path <- file.path(output_path, output_json_name, fsep="/")

output_target_name <- as.character(paste(output_id, ".png", sep = ""))
output_target_path <- file.path(output_path, output_target_name, fsep="/")

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

output_target_wat_name <- as.character(paste(output_id, "-target-wat.png", sep = ""))
output_target_wat_path <- file.path(output_path, output_target_wat_name, fsep="/")

output_target_sed_name <- as.character(paste(output_id, "-target-sed.png", sep = ""))
output_target_sed_path <- file.path(output_path, output_target_sed_name, fsep="/")

data_path <- file.path(dir_data, dir_persp, year, settype, dataset)

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
csf_aoi_shp <- read_sf(dsn = ga_aoi_path) 
aoi_shp <- read_sf(dsn = aoi_path)

# Intersect target with area of interest
targets_aoi_shp <- st_intersection(mctar_all, bounding_box)

# Separate targets_shp filtered by class
mctar_wat <- targets_aoi_shp %>%  filter(Id == 1)
mctar_sed <- targets_aoi_shp %>%  filter(Id == 2)
mctar_veg <- targets_aoi_shp %>%  filter(Id == 3)
mctar_rock <- targets_aoi_shp %>%  filter(Id == 5)

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
#   geom_sf(data = csf_aoi_shp, mapping = aes()) +
#   coord_sf(crs = st_crs(2056))

# Read LAS
# Intensity (ct_wat_i), color information (RGB), number of Returns (r), classification (c)
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
# plot(las, size = 1, color = "RGB", bg = "white")

las_origin <- las

las <- las_origin
par(mfrow=c(1,1))

# good values for sediment (rigidness=1)
# ct_wat_i <- c(0.9, 0.85, 0.8, 0.75)
# clr_wat_j <- c(1.8, 1.7, 1.6, 1.5)
# 
# ct_wat_i <- c(0.4, 0.85)
# clr_wat_j <- c(5.0, 1.7)

# good values for sediment (rigidness=1)
# clr_wat_j <- seq(from = 2.5, to = 7.0, by = 0.1)
# ct_wat_i <- seq(from = 0.14, to = 0.8, by = 0.01)
# clr_wat_j <- seq(from = 3.8, to = 4.0, by = 0.1)
# ct_wat_i <- seq(from = 0.2, to = 0.8, by = 0.2)
# rig_wat_h <- c(1,2,3)
# rig_wat_h <- c(1)

# tbd: tune
# clr_sed_o <- seq(from = 1.9, to = 6.5, by = 0.2)
# ct_sed_n <- seq(from = 0.3, to = 0.9, by = 0.03)
# clr_sed_o <- c(5.5)
# ct_sed_n <- c(1.2)

# rig_sed_m <- c(1,2,3)
# rig_sed_m <- c(1)


df <- as.character(paste("sed_name", "sed_rigidness", "sed_classthreshold",
                 "sed_clothresolution", "sed_steepslope", "sed_kappa",
                 "wat_name", "wat_rigidness", "wat_classthreshold",
                 "wat_clothresolution", "wat_steepslope", "wat_kappa",
                 "n_obs", "comp_time_sec", "rasterresolution", sep =";"))

# Write header before start
output_csv_name <- as.character(paste("genetic_algo_report.csv", sep = ""))
output_csv_path <- file.path(output_path, output_csv_name, fsep="/")
write(df, file=output_csv_path, append = T)

csf_glob_rig <- 3
GA_R3 <- ga(type = "real-valued", 
         fitness =  function(x) -cohen.kappa.csf(las, csf_aoi_shp, targets_aoi_shp, output_path, 
                                                 csf_glob_rig, x[1], x[2], x[3], x[4], x[5]),
         lower = c(0.2, 0.4, 0.2, 2.0, 0.4), 
         upper = c(4, 20, 0.8, 9.0, 0.5), 
         suggestions = c(0.5, 1.9, 0.3, 4.0, 0.5),
         popSize = 1000, maxiter = 180, run = 100,
         maxFitness = 10000,
         optim = TRUE)

csf_glob_rig <- 2
GA_R2 <- ga(type = "real-valued", 
            fitness =  function(x) -cohen.kappa.csf(las, csf_aoi_shp, targets_aoi_shp, output_path, 
                                                    csf_glob_rig, x[1], x[2], x[3], x[4], x[5]),
            lower = c(0.2, 0.4, 0.2, 2.0, 0.4), 
            upper = c(4, 17, 0.8, 9.0, 0.5), 
            suggestions = c(0.5, 1.9, 0.3, 4.0, 0.5),
            popSize = 1000, maxiter = 180, run = 100,
            maxFitness = 10000,
            optim = TRUE)

csf_glob_rig <- 1
GA_R1 <- ga(type = "real-valued", 
            fitness =  function(x) -cohen.kappa.csf(las, csf_aoi_shp, targets_aoi_shp, output_path, 
                                                    csf_glob_rig, x[1], x[2], x[3], x[4], x[5]),
            lower = c(0.2, 0.4, 0.2, 2.0, 0.4), 
            upper = c(4, 17, 0.8, 9.0, 0.5), 
            suggestions = c(0.5, 1.9, 0.3, 4.0, 0.5),
            popSize = 1000, maxiter = 180, run = 100,
            maxFitness = 10000,
            optim = TRUE)

png(output_target_sed_path, height=nrow(tar_sed), width=ncol(tar_sed)) 
plot(tar_sed, maxpixels=ncell(tar_sed), legend =F)
dev.off()

png(output_target_wat_path, height=nrow(tar_wat), width=ncol(tar_wat)) 
plot(tar_wat, maxpixels=ncell(tar_wat), legend =F)
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
