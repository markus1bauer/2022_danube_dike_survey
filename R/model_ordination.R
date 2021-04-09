# Model for ordination ####
# Markus Bauer



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Preparation ################################################################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


### Packages ###
library(here)
library(tidyverse)
library(vegan)

### Start ###
rm(list = ls())
setwd(here("data/processed"))

### Load data ###
sites <- read_csv2("data_processed_sites.csv", col_names = T, na = "na", col_types = 
                    cols(
                      .default = col_guess(),
                      id = col_factor(),
                      location = col_factor(),
                      block = col_factor(),
                      plot = col_factor(),
                      exposition = col_factor(),
                      ffh = col_factor(),
                      vegetationCov = col_double()
                    )) %>%
  select(id, plot, block, location, surveyYear, constructionYear, plotAge, exposition, PC1, PC2, PC3, ffh, changeType, vegetationCov, targetRichness, targetCov, ffh6510Richness, ffh6210Richness) %>%
  mutate(surveyYearF = as_factor(surveyYear))

species <- read_csv2("data_processed_species.csv", col_names = T, na = "na", col_types = 
                      cols(
                        .default = col_double(),
                        name = col_factor()
                      )) %>%  
  pivot_longer(-name, "id", "value") %>%
  pivot_wider(id, name) %>%
  arrange(id) %>%
  column_to_rownames("id")

sitesFFH <- sites %>%
  filter(ffh != "6210")
speciesFFH <- species %>%
  rownames_to_column("id") %>%
  semi_join(sitesFFH, by = "id") %>%
  column_to_rownames("id")
  

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Statistics ################################################################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


### 1 NMDY #####################################################################################

set.seed(1)
(ordi <- metaMDS(species, try = 99, previous.best = T, na.rm = T))
stressplot(ordi)
# stress: 0.286 (wisconsin(sqrt()))
#R²linear = .593
#R²non-metric = .918

### 2 Environmental factors #####################################################################################

#### a Vectors ----------------------------------------------------------------------------------------
(ef_vector1 <- envfit(ordi ~  surveyYear + plotAge + PC1 + PC2 + PC3 + vegetationCov + targetRichness + ffh6510Richness + ffh6210Richness + targetCov, 
              data = sites, permu = 999, na.rm = T))
plot(ordi, type = "n"); plot(ef_vector1, add = T, p. = .99)
(ef_vector2 <- envfit(ordi ~  surveyYear + plotAge + PC1 + PC2 + PC3, 
                      data = sites, permu = 999, na.rm = T))
plot(ordi, type = "n"); plot(ef_vector2, add = T, p. = .99)

#### b factors ----------------------------------------------------------------------------------------
(ef_factor1 <- envfit(ordi ~  plot + block + surveyYear + exposition + location + ffh + changeType, 
              data = sites, permu = 999, na.rm = T))
plot(ordi, type = "n"); ordiellipse(ordi, sites$location, kind = "sd", draw = "lines", label = T)
plot(ordi, type = "n"); ordiellipse(ordi, sites$plot, kind = "sd", draw = "lines", label = T)
plot(ordi, type = "n"); ordiellipse(ordi, sites$surveyYear, kind = "sd", draw = "lines", label = T)
plot(ordi, type = "n"); ordiellipse(ordi, sites$constructionYear, kind = "sd", draw = "lines", label = T)
plot(ordi, type = "n"); ordiellipse(ordi, sites$exposition, kind = "sd", draw = "lines", label = T)
plot(ordi, type = "n"); ordiellipse(ordi, sites$ffh, kind = "sd", draw = "lines", label = T)
plot(ordi, type = "n"); ordiellipse(ordi, sites$changeType, kind = "sd", draw = "lines", label = T)

