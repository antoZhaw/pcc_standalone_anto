#### Quantitative assessment of terrestrial sediment budgets in a ####
# hydropower-regulated floodplain using point cloud classification
# Documentation:  
# https://github.com/gubely/math_documentation/blob/main/main_math_documentation.pdf

# Author: Yves Gubelmann
# Contact: yves.gubelmann@gmx.ch

# Description: This script searches the best parameter set for a TLS point cloud
# within a given parameter window. For this, the developed terrestrial sediment 
# classification routine based on cloth simulation filter (CSF) uses a genetic
# algorithm to classify the provided raw point cloud and compare the intermediate
# result with a reference map. Intermediate results are saved in a .csv-file.

# libraries ####
library(lidR) # Point cloud classification
library(papeR) # summary tables
library(tidyverse) # tidy essentials (ggplot, purr, tidyr, readr, dplyr)
library(lubridate) # handling dates and time
library(terra) # handling spatial data
library(sf) # handling spatial data
library(rgdal) # export raster objects
library(ncdf4) # export ncdf files
library(janitor) # clean and consistent naming
library(forcats) # handling factor levels
library(raster) # rasterizing vector data
library(RCSF) # for CSF ground classification
library(geometry) # for raserize_canopy function
library(purrr) # for map function
library(rjson) # for JSON generation
library(rgl) # for RGL Viewer control functions
library(fasterize) # for faster rasterization
library(psych) # for cohen's kappa
library(GA) # genetic algorithm for parameter search optimisation


# Functions ####

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

cohen.kappa.csf.sed.tls <- function(raw_las, ga_aoi_shp, targets_shp, 
                            ga_output_path, output_csv_path,
                            las_year, las_persp, las_dataset_id,
                            rig_sed_m, ct_sed_n, clr_sed_o, 
                            raster_res) {
  id_ij <- sample(1:999, 1)
  start_ij <- as_datetime(lubridate::now())
  rig_sed_m <- round(rig_sed_m, 0)
  steep_slope_sed <- if_else(rig_sed_m == 3, T, F)
  # classify ground
  msg_sed <- as.character(paste("SED", id_ij, "_rig", round(rig_sed_m, 4), "_ct", round(ct_sed_n, 4), "_clr", round(clr_sed_o, 4), sep = ""))
  print(msg_sed)
  las_sed_ij <- classify.gnd(las, steep_slope_sed, ct_sed_n, clr_sed_o, rig_sed_m)
  las_sed_ij <- classify_poi(las_sed_ij, class = LASNOISE, roi = ga_aoi_shp, inverse_roi = T)
  las_sed_ij <- filter_poi(las_sed_ij, Classification != LASNOISE)
  # Rasterize both point clouds    
  DEM_sed_ij <- rasterize_canopy(las_sed_ij, res = raster_res, p2r(), pkg = "raster")
  sed_ij <- DEM_sed_ij
  # Save plot of masked raster
  output_sed_rast_name <- as.character(paste(msg_sed, ".png", sep = ""))
  output_sed_rast_path <- file.path(ga_output_path, output_sed_rast_name, fsep="/")
  png(output_sed_rast_path, height=nrow(sed_ij), width=ncol(sed_ij))
  plot(sed_ij, maxpixels=ncell(sed_ij), legend =F)
  dev.off()
  # Generate raster for sediment comparison
  raster_ext_sed <- extent(xmin(sed_ij), xmax(sed_ij), ymin(sed_ij), ymax(sed_ij))
  tar_raw_sed <- raster(nrows=nrow(sed_ij), ncols=ncols(sed_ij), crs=2056,
                        ext=raster_ext_sed, resolution=raster_res, vals=NULL)  
  target_sed <- targets_shp %>%  filter(Id == 2)
  tar_sed <- fasterize(target_sed, tar_raw_sed, field = "Id", fun="sum")
  # Normalise raster values for comparison
  rater1 <- values(tar_sed)
  rater1[rater1 != 0] <- 1
  rater1[is.na(rater1)] <- 0
  rater2 <- values(sed_ij)
  rater2[rater2 != 0] <- 1
  rater2[is.na(rater2)] <- 0
  kap_sed <- cohen.kappa(x=cbind(rater1,rater2))
  end_ij <- as_datetime(lubridate::now())
  timespan_ij <- interval(start_ij, end_ij)
  delta_t <- as.numeric(timespan_ij, "seconds")
  total_kappa <- kap_sed$kappa
  iter_result_msg <- as.character(paste("Total Kappa (sed): ", round(total_kappa, 4), ", computed in ", round(delta_t, 3), " seconds."))
  print(iter_result_msg)
  obs <- as.character(paste(msg_sed, rig_sed_m, ct_sed_n, clr_sed_o, steep_slope_sed, kap_sed$kappa,
                            NA, NA, NA, NA, NA, NA, 
                            kap_sed$n.obs, delta_t, raster_res,
                            NA, NA, NA, NA,
                            xmin(raster_ext_sed), xmax(raster_ext_sed),
                            ymin(raster_ext_sed), ymax(raster_ext_sed),
                            las_year, las_persp, las_dataset_id,
                            ga_output_path,
                            sep =";"))
  write(obs, file=output_csv_path, append = T)
  popsize_kappa <- total_kappa * 10000
  return(popsize_kappa)
}

# Globals for Configuration ####
# Specify dataset
dataset_id <- "4"
wholeset <- T
year <- "2022"
perspective <- "tls"
settype <- if_else(wholeset == T, "wholeset", "subset")

# Internal globals such as paths and IDs ####
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
  dir_data <- "C:/code_wc/math_gubelyve"
}

dir_persp <- if_else(perspective == "tls", "tls_data", "uav_data")
dir_config <-  file.path(dir_repo, "config", fsep="/")
dir_export <-  file.path(dir_repo, "export", fsep="/")

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

output_target_all_name <- as.character(paste(output_id, "-target-all.png", sep = ""))
output_target_all_path <- file.path(output_path, output_target_all_name, fsep="/")

output_ga_sed_report_name <- as.character(paste("ga_sed", perspective, year, dataset_id, "report.csv", sep = "_"))
output_ga_sed_report_path <- file.path(dir_export, output_ga_sed_report_name, fsep="/")

data_path <- file.path(dir_data, dir_persp, year, settype, dataset)

# Read files ####

# empty warnings if existing.
if(length(warnings())!=0){
  assign("last.warning", NULL, envir = baseenv())
}

# Generate rectangular polygon of area of interest (relevant for uav)
# gen_xy <- structure(list(dat = c("BURR-1-1-1", "BURR-1-1-1", 
#                                  "BURR-1-1-1", "BURR-1-1-1"),
#                          Longitude = c(2575011, 2575011, 
#                                       2575517.4, 2575517.4),
#                          Latitude = c(1178366.7, 1178993, 
#                                        1178993, 1178366.7)),
#                     class = "data.frame", row.names = c(NA,-4L))
# 
# bounding_box <- gen_xy %>%
#   st_as_sf(coords = c("Longitude", "Latitude"), crs = 2056) %>%
#   group_by(dat) %>%
#   summarise(geometry = st_combine(geometry)) %>%
#   st_cast("POLYGON")

# Narrowed tls polygon
gen_xy_tls <- structure(list(dat = c("AOI TLS", "AOI TLS", 
                                 "AOI TLS", "AOI TLS"),
                         Longitude = c(2575340, 2575340, 
                                       2575480, 2575480),
                         Latitude = c(1178520, 1178780, 
                                      1178780, 1178520)),
                    class = "data.frame", row.names = c(NA,-4L))

bounding_box_tls <- gen_xy_tls %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 2056) %>%
  group_by(dat) %>%
  summarise(geometry = st_combine(geometry)) %>%
  st_cast("POLYGON")


# Read Shapefiles
mctar_all <- read_sf(dsn = mctar_path)

if(year == "2022"){
  mctar_all <- mctar_all %>% 
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

csf_aoi_shp <- read_sf(dsn = ga_aoi_path) 

# Intersect target with area of interest
targets_aoi_shp <- st_intersection(mctar_all, bounding_box_tls)

ggplot() +
  geom_sf(data = bounding_box_tls, aes(fill = dat), alpha = 0.1) +
  geom_sf(data = targets_aoi_shp, aes(fill = as.factor(Class_name)), alpha = 0.5) +
  geom_sf(data = csf_aoi_shp, aes(fill = as.factor(id)), alpha = 0.5) +
  coord_sf(datum=st_crs(2056)) +
  theme_minimal()

ggsave(output_target_all_path, bg = 'white')

# Information about rasterized target
# target$layer
# st_crs(target)

# Read LAS
# Intensity (ct_wat_i), color information (RGB), number of Returns (r), classification (c)
# of the first point is loaded only to reduce computational time.
las <- readLAS(data_path, select = "xyzRGBc", filter = cfg$las_filter)

# Reset class LASNOISE for further procedure
las$Classification <- LASNONCLASSIFIED


# Check LAS whether it complies with the required
if(has.lasClassification(las)){
  print("LAS file is already classified. Are you sure to continue?")}
# if (length(warnings())>=1) {stop("The read LAS file throws warnings, script stops.")}

par(mfrow=c(1,1))

# Genetic algorithm for parameter optimisation ####
# Write header before start
df <- as.character(paste("sed_name", "sed_rigidness", "sed_classthreshold",
                 "sed_clothresolution", "sed_steepslope", "sed_kappa",
                 "wat_name", "wat_rigidness", "wat_classthreshold",
                 "wat_clothresolution", "wat_steepslope", "wat_kappa",
                 "n_obs", "comp_time_sec", "rasterresolution",
                 "wat_xmin", "wat_xmax", "wat_ymin", "wat_ymax",
                 "sed_xmin", "sed_xmax", "sed_ymin", "sed_ymax", "year", 
                 "perspective", "dataset_id", "output_path", sep =";"))


write(df, file=output_ga_sed_report_path, append = T)

# Start search ####
gar3_start <- as_datetime(lubridate::now())
GA <- ga(type = "real-valued", 
         fitness =  function(x) -cohen.kappa.csf.sed.tls(las, csf_aoi_shp, targets_aoi_shp, 
                                                 output_path, output_ga_sed_report_path,
                                                 year, perspective, dataset_id,
         3, x[1], x[2], 0.2),
         lower = c(0.08, 1.2), 
         upper = c(1, 5), 
         suggestions = c(0.246668090, 3.78168753),
         popSize = 1000, maxiter = 50, run = 10,
         maxFitness = 10000,
         optim = F)

gar3_end <- as_datetime(lubridate::now())
timespan <- interval(gar3_start, gar3_end)
timespan