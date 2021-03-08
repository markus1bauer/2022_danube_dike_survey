# Show figure 2 ####
# Markus Bauer
# Citation: Markus Bauer 



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Preparation ################################################################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


### Packages ###
library(tidyverse)
library(adespatial)

### Start ###
rm(list = ls())
setwd("Z:/Documents/0_Donaudeiche/3_Aufnahmen_und_Ergebnisse/2022_Danube_old_dikes/data/processed")


### Load data ###
speciesPa <- read_csv2("data_processed_tbiPa.csv", col_names = T, na = "na", col_types = 
                         cols(
                           .default = col_double(),
                           id = col_factor()
                           )
)
speciesPa <- column_to_rownames(speciesPa, var = "id")
speciesAbu <- read_csv2("data_processed_tbiAbu.csv", col_names = T, na = "na", col_types = 
                          cols(
                            .default = col_double(),
                            id = col_factor()
                            )
)
speciesAbu <- column_to_rownames(speciesAbu, var = "id")

### Make sub data frames of 2017 and 2019 ###
speciesPa17 <- speciesPa[1:62,]
speciesPa19 <- speciesPa[63:124,]

speciesAbu17 <- speciesAbu[1:62,]
speciesAbu19 <- speciesAbu[63:124,]
rm(speciesAbu, speciesPa)

### Models ###
resPa <- TBI(speciesPa17, speciesPa19, method = "sorensen", 
              nperm = 9999, test.t.perm = T, clock = T)
resAbu <- TBI(speciesAbu17, speciesAbu19, method = "%diff", 
              nperm = 9999, test.t.perm = T, clock = T)


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# C Plots ################################################################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


par(mfrow = c(1,2))
plot(resPa, type = "BC", xlim = c(0,0.6), ylim = c(0,0.5), 
     col.bg = "grey50", main = "Species richness 2017-2019")
plot(resAbu, type = "BC", xlim = c(0,0.6), ylim = c(0,0.5), 
     col.bg = "grey50", main = "Species abundance 2017-2019")

