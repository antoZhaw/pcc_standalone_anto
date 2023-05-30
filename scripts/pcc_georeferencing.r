# libraries ####
library(tidyverse)
library(janitor)
library(lubridate)

# globals ####
user <- Sys.getenv("USERNAME")
dir_repo <- if_else(user == "gubelyve", "C:/Daten/math_gubelyve/pcc_standalone", "C:/code_wc/pcc_standalone")
dir_data <- file.path(dir_repo, "data")

# read data ####
dt <- read_delim(file.path(dir_data, "230106_Georeferenzierung_TLS_2021.csv"),",",
                 locale = locale(encoding = "UTF-8",decimal_mark = "."), delim = ";", 
                 col_types = cols()) %>%
  janitor::clean_names() %>% # säubert die Spaltennamen 
  mutate(delta_x = gps_x - ply_x) %>% 
  mutate(delta_y = gps_y - ply_y) %>% 
  mutate(delta_z = gps_z - ply_z) %>% 
  mutate(class = as.factor("saane_20211013_merged_subs_noSP1"))
summary(dt)

p_delta_x <- ggplot(data = dt, mapping = aes(x = class, y = delta_x, color = gps_id)) +
  geom_point(size=1) +
  labs(x = "Datensatz",
       y = "Delta der Koordinaten", 
       title = "Delta der X Koordinate",
       subtitle = "GPS Sample minus .ply Messung ergibt die Translation.")
p_delta_x


p_delta_y <- ggplot(data = dt, mapping = aes(x = class, y = delta_y, color = gps_id)) +
  geom_point(size=1) +
  labs(x = "Datensatz",
       y = "Delta der Koordinaten", 
       title = "Delta der Y Koordinate",
       subtitle = "GPS Sample minus .ply Messung ergibt die Translation.")
p_delta_y

p_delta_z <- ggplot(data = dt, mapping = aes(x = class, y = delta_z, color = gps_id)) +
  geom_point(size=1) +
  labs(x = "Datensatz",
       y = "Delta der Koordinaten", 
       title = "Delta der Z Koordinate",
       subtitle = "GPS Sample minus .ply Messung ergibt die Translation.")
p_delta_z

mean(dt$delta_x)
mean(dt$delta_y)
mean(dt$delta_z)
