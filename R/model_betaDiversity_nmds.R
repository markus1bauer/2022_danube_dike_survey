# Model for NMDS ####
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
sites <- read_csv("data_processed_sites.csv", col_names = T, na = c("na", "NA", ""), col_types = 
                    cols(
                      .default = "d",
                      id = "f",
                      locationAbb = "f",
                      block = "f",
                      plot = "f",
                      exposition = "f"
                    )) %>%
  select(id, plot, block, locationAbb, surveyYear, exposition, vegetationCov, targetCov, graminoidCovratio, speciesRichness, targetRichness, shannon, eveness, accumulatedCov, syn_total, syn_trend, syn_detrend) %>%
  mutate(surveyYearF = as_factor(surveyYear)) %>%
  filter(accumulatedCov > 0)

species <- read_csv("data_processed_species.csv", col_names = T, na = c("na", "NA", ""), col_types = 
                      cols(
                        .default = "d",
                        name = "f"
                      )) %>%  
  mutate(across(where(is.numeric), ~replace(., is.na(.), 0))) %>%
  pivot_longer(-name, "id", "value") %>%
  pivot_wider(id, name) %>%
  arrange(id) %>%
  semi_join(sites, by = "id") %>%
  column_to_rownames("id")



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Statistics ################################################################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


### 1 NMDY #####################################################################################

### Claculate ###
set.seed(1)
(ordi <- metaMDS(species, dist = "bray", binary = F,
                 try = 99, previous.best = T, na.rm = T))
### Stress ###
stressplot(ordi)
#wisconsin(sqrt(species))
#stress: .281
#R²linear = .612
#R²non-metric = .920
### Goodness of fit ###
gof <- goodness(ordi)
plot(ordi, type = "t", main = "Goodness of fit")
points(ordi, display = "sites", cex = gof * 300)



### 2 Environmental factors #####################################################################################

#### a Vectors ----------------------------------------------------------------------------------------
(ef_vector1 <- envfit(ordi ~  vegetationCov + graminoidCovratio + targetRichness + targetCov + speciesRichness + eveness + shannon + syn_total + syn_trend + syn_detrend, 
              data = sites, permu = 999, na.rm = T))
plot(ordi, type = "n"); plot(ef_vector1, add = T, p. = .99)
(ef_vector2 <- envfit(ordi ~  graminoidCovratio + speciesRichness + eveness, 
                      data = sites, permu = 999, na.rm = T))
plot(ordi, type = "n"); plot(ef_vector2, add = T, p. = .99)

#### b Factors ----------------------------------------------------------------------------------------
(ef_factor1 <- envfit(ordi ~  surveyYearF + exposition, 
              data = sites, permu = 999, na.rm = T))
plot(ordi, type = "n"); ordiellipse(ordi, sites$surveyYearF, kind = "sd", draw = "lines", label = T)
plot(ordi, type = "n"); ordiellipse(ordi, sites$exposition, kind = "sd", draw = "lines", label = T)

