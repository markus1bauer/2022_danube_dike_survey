# Show map of the Danube dike experiment ####
# Markus Bauer
# Citation: Markus Bauer, Jakob Huber, Katharina Beck, Johannes Kollmann (xxxx)



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Preparation ################################################################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

### Packages ###
library(tidyverse)
library(soiltexture)

### Start ###
rm(list = ls())
setwd("Z:/Documents/0_Donaudeiche/3_Aufnahmen_und_Ergebnisse/2022_Danube_old_dikes/data/raw")


### Load data ###
sites <- read_csv2("data_raw_sites_2017_2018_2019.csv", col_names = T, na = "na", col_types = 
                     cols(
                       .default = col_double(),
                       id = col_factor(),
                       location = col_factor(),
                       ageCategory = col_factor(),
                       HCl = col_factor(),
                       phosphorousClass = col_factor(),
                       potassiumClass = col_factor(),
                       magnesiumClass = col_character()
                     )        
)

sites <- sites %>%
  select(sand, silt, clay, phosphorous) %>%
  rename(SAND = sand, SILT = silt, CLAY = clay)
sites <- as.data.frame(sites)



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Plot ##############################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


TT.plot(
  class.sys = "DE.BK94.TT", 
  tri.data = sites, 
  z.name = "phosphorous")

TT.plot(
  class.sys = "DE.BK94.TT", 
  tri.data = sites
  )

