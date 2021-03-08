# Calculate TBI ####
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
                     ) %>%
  column_to_rownames(var = "id")
speciesAbu <- read_csv2("data_processed_tbiAbu.csv", col_names = T, na = "na", col_types = 
                          cols(
                            .default = col_double(),
                            id = col_factor()
                          )
                        ) %>%
  column_to_rownames(var = "id")

### Make sub data frames of 2017 and 2019 ###
half <- as.numeric(count(speciesPa)/2)
speciesPa17 <- slice_head(speciesPa, n = half)
speciesPa19 <- slice_tail(speciesPa, n = half)
speciesAbu17 <- slice_head(speciesAbu, n = half)
speciesAbu19 <- slice_tail(speciesAbu, n = half)
rm(speciesAbu, speciesPa, half)



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Statistics ################################################################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


## 1 Dataset 1984-2018 #####################################################################################

### 1.1 Presence-absence data 2017-2019 --------------------------------------------------------------------------------------------
####Total
resPa <- TBI(speciesPa17, speciesPa19, method = "sorensen", 
              nperm = 9999, test.t.perm = T, clock = T)
#### FFH6510 vs. FFH 6210
#resPa65 <- TBI(speciesPa17, speciesPa19, method = "sorensen",
#               nperm = 9999, test.t.perm = T, clock = T)
#resPa62 <- TBI(speciesPa17, speciesPa19, method = "sorensen", 
#               nperm = 9999, test.t.perm = T, clock = T)
####Test plots
par(mfrow = c(1,1))
plot(resPa, type = "BC")
resPa #xxx
#par(mfrow = c(1,3))
#plot(resPa65, type = "BC", xlim = c(0,0.8), ylim = c(0,0.8))
#plot(resPa62, type = "BC", xlim = c(0,0.8), ylim = c(0,0.8))
#resPa65 #xxx
#resPa62 #xxx

### 1.2 Abundance data 2017-2019 --------------------------------------------------------------------------------------------
####Total
resAbu <- TBI(speciesAbu17, speciesAbu19, method = "%diff", 
               nperm = 9999, test.t.perm = T, clock = T)

resAbu$BCD.summary #xxx
resAbu$t.test_B.C #xxx
bcdAbu <- as_tibble(resAbu$BCD.mat)
bcdAbu$id <- speciesAbu17 %>%
  rownames_to_column(var = "id") %>%
  select(id); bcdAbu

#### Test plots
par(mfrow = c(1,1))
plot(resAbu, type = "BC")
