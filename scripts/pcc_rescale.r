library(lidR) # summary tables

# Globals-----------------------------------------------------------------------

data_path <- "C:/Daten/math_gubelyve/tls_data/2021/subset/saane_20211013_subsample_onlyRGBpts.las"
output_las_path <- "C:/Daten/math_gubelyve/tls_data/2021/subset/saane_20211013_subsample_onlyRGBpts_rescaled.las"

# Read LAS file-----------------------------------------------------------------

las <- readLAS(data_path, select = "xyzRGBci", filter = "-keep_first")
las <- las_rescale(las, xscale = 10e-08, yscale = 10e-08, zscale = 1e-07)

# Save rescaled LAS file
writeLAS(las, file = output_las_path)
