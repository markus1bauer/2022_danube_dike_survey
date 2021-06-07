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
species <- read_csv("data_processed_species.csv", col_names = T, na = "na", col_types = 
                      cols(
                        .default = "d",
                        name = "f"
                      )) %>%  
  pivot_longer(-name, "id", "value") %>%
  pivot_wider(id, name) %>%
  arrange(id) 

sites <- read_csv("data_processed_sites.csv", col_names = T, na = "na", col_types = 
                    cols(
                      .default = "d",
                      id = "f",
                      locationAbb = "f",
                      block = "f",
                      plot = "f",
                      exposition = "f",
                      ffh = "f",
                      changeType = "f"
                    )) %>%
  select(id, plot, block, locationAbb, surveyYear, constructionYear, plotAge, exposition, PC1, PC2, PC3, ffh, changeType, vegetationCov, targetCov, graminoidCovratio, speciesRichness, targetRichness, ffh6510Richness, ffh6210Richness, shannon, eveness) %>%
  mutate(surveyYearF = as_factor(surveyYear)) %>%
  semi_join(species, by = "id")

species <- species %>%
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
# stress: 0.286
#R²linear = .593
#R²non-metric = .918
### Goodness of fit ###
gof <- goodness(ordi)
plot(ordi, type = "t", main = "Goodness of fit")
points(ordi, display = "sites", cex = gof * 300)
### Plots ###
plot(ordi, display = "species")
plot(ordi, display = "sites")



### 2 Environmental factors #####################################################################################

#### a Vectors ----------------------------------------------------------------------------------------
(ef_vector1 <- envfit(ordi ~  surveyYear + constructionYear + plotAge + PC1 + PC2 + PC3 + vegetationCov + graminoidCovratio + targetRichness + ffh6510Richness + ffh6210Richness + targetCov + speciesRichness + eveness + shannon, 
              data = sites, permu = 999, na.rm = T))
plot(ordi, type = "n"); plot(ef_vector1, add = T, p. = .99)
(ef_vector2 <- envfit(ordi ~  surveyYear + plotAge + PC1 + PC2 + PC3, 
                      data = sites, permu = 999, na.rm = T))
plot(ordi, type = "n"); plot(ef_vector2, add = T, p. = .99)

#### b Factors ----------------------------------------------------------------------------------------
(ef_factor1 <- envfit(ordi ~  plot + block + surveyYearF + locationAbb + exposition + ffh + changeType, 
              data = sites, permu = 999, na.rm = T))
plot(ordi, type = "n"); ordiellipse(ordi, sites$locationAbb, kind = "sd", draw = "lines", label = T)
plot(ordi, type = "n"); ordiellipse(ordi, sites$plot, kind = "sd", draw = "lines", label = T)
plot(ordi, type = "n"); ordiellipse(ordi, sites$surveyYear, kind = "sd", draw = "lines", label = T)
plot(ordi, type = "n"); ordiellipse(ordi, sites$constructionYear, kind = "sd", draw = "lines", label = T)
plot(ordi, type = "n"); ordiellipse(ordi, sites$exposition, kind = "sd", draw = "lines", label = T)
plot(ordi, type = "n"); ordiellipse(ordi, sites$ffh, kind = "sd", draw = "lines", label = T)
plot(ordi, type = "n"); ordiellipse(ordi, sites$changeType, kind = "sd", draw = "lines", label = T)

