# libraries --------------------------------------------------------------------

# install packages

library(tidyverse)
library(janitor)    #saeubert spaltennamen
library(lubridate)

# globals-----------------------------------------------------------------------
user <- Sys.getenv("USERNAME")

dir_repo <- if_else(user == "gubelyve", "C:/Daten/math_gubelyve/pcc_standalone", "C:/code_wc/pcc_standalone")
dir_data <- file.path(dir_repo, "data")

# read data---------------------------------------------------------------------

dt <- read_delim(file.path(dir_data, "221205_Klassifikation.csv"),",",
                 locale = locale(encoding = "UTF-8"), delim = ";", 
                 col_types = cols()) %>%
  janitor::clean_names() %>% # säubert die Spaltennamen 
  mutate(RtoG = r/g) %>% 
  mutate(RtoB = r/b) %>% 
  mutate(GtoB = g/b) %>% 
  mutate(RBtimesGB = RtoB * GtoB) %>% 
  mutate(RGplusB = RtoB + GtoB) %>% 
  mutate(
    main_class = case_when(
      class == "sky"~"sky",
      class == "dark blue sky"~"sky",
      class == "sediment"~"sediment",
      class == "vegetation general"~"vegetation",
      class == "vegetation bright"~"vegetation",
      class == "vegetation dark"~"vegetation",
      class == "cliff bright"~"cliff",
      class == "cliff blue"~"cliff",
      class == "cliff dark"~"cliff",
      TRUE~"default" #Default case
    )
  )

summary(dt)

pRBtimesGB <- ggplot(data = dt, mapping = aes(x = class, y = RBtimesGB, colour = main_class)) +
  geom_boxplot(size=0.5) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_y_continuous(breaks = c(0,1,2,3,4,5), limits = c(0,5)) +
  labs(x = "Klasse",
       y = "R zu G", 
       title = "R/G multipliziert mit B/G nach Klasse",
       subtitle = "Manuell zugewiesene Klassen")

pRBtimesGB

pRGplusB <- ggplot(data = dt, mapping = aes(x = class, y = RGplusB, colour = main_class)) +
  geom_boxplot(size=0.5) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_y_continuous(breaks = c(0,1,2,3,4,5), limits = c(0,5)) +
  labs(x = "Klasse",
       y = "R zu G", 
       title = "Verhältnis R/G nach Klasse",
       subtitle = "Manuell zugewiesene Klassen")

pRGplusB

pscat <- ggplot(data = dt, mapping = aes(x = g, y = b, colour = main_class)) +
  geom_point(size=2) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  # scale_y_continuous(breaks = c(0,50,100,150,200,250), limits = c(0,255)) +
  labs(x = "Klasse",
       y = "R zu G", 
       title = "Verhältnis R/G nach Klasse",
       subtitle = "Manuell zugewiesene Klassen")
pscat

pRtoG <- ggplot(data = dt, mapping = aes(x = class, y = RtoG, colour = main_class)) +
  geom_boxplot(size=0.5) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  # scale_y_continuous(breaks = c(0,50,100,150,200,250), limits = c(0,255)) +
  labs(x = "Klasse",
       y = "R zu G", 
       title = "Verhältnis R/G nach Klasse",
       subtitle = "Manuell zugewiesene Klassen")

pRtoB <- ggplot(data = dt, mapping = aes(x = class, y = RtoB, colour = main_class)) +
  geom_boxplot(size=0.5) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_y_continuous(breaks = c(0,1,2,3,4,5), limits = c(0,5)) +
  labs(x = "Klasse",
       y = "R zu B", 
       title = "Verhältnis R/B nach Klasse",
       subtitle = "Manuell zugewiesene Klassen")

pGtoB <- ggplot(data = dt, mapping = aes(x = class, y = GtoB, colour = main_class)) +
  geom_boxplot(size=0.5) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_y_continuous(breaks = c(0,1,2,3,4,5), limits = c(0,5)) +
  labs(x = "Klasse",
       y = "G zu B", 
       title = "Verhältnis G/B nach Klasse",
       subtitle = "Manuell zugewiesene Klassen")

pR <- ggplot(data = dt, mapping = aes(x = class, y = r, colour = main_class)) +
  geom_boxplot(size=0.5) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_y_continuous(breaks = c(0,50,100,150,200,250), limits = c(0,255)) +
  labs(x = "Klasse",
       y = "R Anteil (0 bis 255)", 
       title = "Farbanteile nach Klasse - Rot",
       subtitle = "Manuell zugewiesene Klassen")

pG <- ggplot(data = dt, mapping = aes(x = class, y = g, colour = main_class)) +
  geom_boxplot(size=0.5) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_y_continuous(breaks = c(0,50,100,150,200,250), limits = c(0,255)) +
  labs(x = "Klasse",
       y = "G Anteil (0 bis 255)", 
       title = "Farbanteile nach Klasse - Grün",
       subtitle = "Manuell zugewiesene Klassen")

pB <- ggplot(data = dt, mapping = aes(x = class, y = b, colour = main_class)) +
  geom_boxplot(size=0.5) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_y_continuous(breaks = c(0,50,100,150,200,250), limits = c(0,255)) +
  labs(x = "Klasse",
       y = "B Anteil (0 bis 255)", 
       title = "Farbanteile nach Klasse - Blau",
       subtitle = "Manuell zugewiesene Klassen")


pRmain <- ggplot(data = dt, mapping = aes(x = main_class, y = r, colour = main_class)) +
  geom_boxplot(size=0.5) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_y_continuous(breaks = c(0,50,100,150,200,250), limits = c(0,255)) +
  labs(x = "Klasse",
       y = "R Anteil (0 bis 255)", 
       title = "Farbanteile nach Klasse - Rot",
       subtitle = "Manuell zugewiesene Klassen")

pGmain <- ggplot(data = dt, mapping = aes(x = main_class, y = g, colour = main_class)) +
  geom_boxplot(size=0.5) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_y_continuous(breaks = c(0,50,100,150,200,250), limits = c(0,255)) +
  labs(x = "Klasse",
       y = "G Anteil (0 bis 255)", 
       title = "Farbanteile nach Klasse - Grün",
       subtitle = "Manuell zugewiesene Klassen")

pBmain <- ggplot(data = dt, mapping = aes(x = main_class, y = b, colour = main_class)) +
  geom_boxplot(size=0.5) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_y_continuous(breaks = c(0,50,100,150,200,250), limits = c(0,255)) +
  labs(x = "Klasse",
       y = "B Anteil (0 bis 255)", 
       title = "Farbanteile nach Klasse - Blau",
       subtitle = "Manuell zugewiesene Klassen")

# pRmain
# pGmain
# pBmain

pR
pG
pB

pRtoG
pRtoB
pGtoB

dt_cliff_blue <- dt %>% filter(class == "cliff blue")
dt_cliff_bright <- dt %>% filter(class == "cliff bright")
dt_cliff_dark <- dt %>% filter(class == "cliff dark")

dt_sky <- dt %>% filter(class == "sky")
dt_darksky <- dt %>% filter(class == "dark blue sky")

dt_veg_dark <- dt %>% filter(class == "vegetation dark")
dt_veg_bright <- dt %>% filter(class == "vegetation bright")
dt_veg_gen <- dt %>% filter(class == "vegetation general")

dt_sediment <- dt %>% filter(class == "sediment")

# Cliff
summary(dt_cliff_blue)
summary(dt_cliff_bright)
summary(dt_cliff_dark)

# Sediment
summary(dt_sediment)

# Sky
summary(dt_sky)
summary(dt_darksky)

# Vegetation
summary(dt_veg_dark)
summary(dt_veg_bright)
summary(dt_veg_gen)

