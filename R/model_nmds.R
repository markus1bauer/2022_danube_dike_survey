# Model for NMDS ####
# Markus Bauer



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Preparation ################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


### Packages ###
library(here)
library(tidyverse)
library(vegan)

### Start ###
rm(list = ls())
setwd(here("data", "processed"))

### Load data ###
sites_dikes <- read_csv("data_processed_sites_spatial.csv", col_names = TRUE,
                  na = c("na", "NA", ""), col_types =
                    cols(
                      .default = "?",
                      id = "f",
                      locationAbb = "f",
                      block = "f",
                      plot = "f",
                      exposition = "f"
                    )) %>%
  select(id, surveyYear, locationYear, latitude, longitude, riverkm,
         distanceRiver, constructionYear, plotAge, exposition, side, PC1soil,
         PC2soil, speciesRichness, accumulatedCov) %>%
  mutate(surveyYearF = as_factor(surveyYear)) %>%
  filter(accumulatedCov > 0)

sites_splot <- read_csv("data_processed_sites_splot.csv", col_names = TRUE,
                  na = c("na", "NA", ""), col_types =
                    cols(
                      .default = "?"
                    ))

sites <- sites_dikes %>%
  bind_rows(sites_splot) %>%
  mutate(reference = if_else(is.na(reference), "no", reference))

species_dikes <- read_csv("data_processed_species.csv", col_names = TRUE,
                    na = c("na", "NA", ""), col_types =
                      cols(
                        .default = "d",
                        name = "f"
                      ))

species_splot <- read_csv("data_processed_species_splot.csv", col_names = TRUE,
                        na = c("na", "NA", ""), col_types =
                          cols(
                            .default = "?"
                          ))

species <- species_dikes %>%
  full_join(species_splot, by = "name") %>%
  mutate(across(where(is.numeric), ~replace(., is.na(.), 0))) %>%
  pivot_longer(cols = -name, names_to = "id", values_to = "value") %>%
  pivot_wider(names_from = "name", values_from = "value") %>%
  arrange(id) %>%
  semi_join(sites, by = "id") %>%
  column_to_rownames("id")


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Statistics ################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


### 1 NMDY ####################################################################

### Calculate ###
set.seed(1)
(ordi <- metaMDS(species, dist = "bray", binary = FALSE,
                 try = 99, previous.best = TRUE, na.rm = TRUE))
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



### 2 Environmental factors ###################################################

#### a Vectors ----------------------------------------------------------------
(ef_vector1 <- envfit(ordi ~
                        speciesRichness + eveness + shannon +
                        accumulatedCov + graminoidCovratio +
                        targetRichness + targetCovratio,
              data = sites, 
              permu = 999, 
              na.rm = TRUE))
plot(ordi, type = "n")
plot(ef_vector1, add = TRUE, p. = .99)
(ef_vector2 <- envfit(ordi ~
                        speciesRichness +
                        accumulatedCov, 
                      data = sites, 
                      permu = 999, 
                      na.rm = T))
plot(ordi, type = "n")
plot(ef_vector2, add = TRUE, p. = .99)

#### b Factors ----------------------------------------------------------------
(ef_factor1 <- envfit(ordi ~  surveyYearF + exposition + reference, 
              data = sites, permu = 999, na.rm = TRUE))
plot(ordi, type = "n")
ordiellipse(ordi, sites$surveyYearF, kind = "sd", draw = "lines", label = TRUE)
plot(ordi, type = "n")
ordiellipse(ordi, sites$exposition, kind = "sd", draw = "lines", label = TRUE)
plot(ordi, type = "p")
ordiellipse(ordi, sites$reference, kind = "sd", draw = "lines", label = TRUE)
