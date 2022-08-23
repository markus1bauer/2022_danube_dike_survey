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
sites_dikes <- read_csv("data_processed_sites_spatial_nmds.csv",
                        col_names = TRUE, na = c("na", "NA", ""), col_types =
                    cols(
                      .default = "?",
                      id = "f"
                    )) %>%
  select(id, survey_year, species_richness, eveness, shannon, target_richness,
         accumulated_cover, target_cover_ratio, graminoid_cover_ratio) %>%
  mutate(survey_year_factor = as_factor(survey_year))

sites_splot <- read_csv("data_processed_sites_splot.csv", col_names = TRUE,
                  na = c("na", "NA", ""), col_types =
                    cols(
                      .default = "?"
                    ))

sites <- sites_dikes %>%
  bind_rows(sites_splot) %>%
  mutate(reference = if_else(is.na(reference), "no", reference),
         esy = if_else(is.na(esy), "dike", esy))

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

rm(list = setdiff(ls(), c("sites", "species")))

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
                        species_richness + #eveness + shannon +
                        accumulated_cover + graminoid_cover_ratio +
                        target_richness + target_cover_ratio,
              data = sites, 
              permu = 999, 
              na.rm = TRUE))
plot(ordi, type = "n")
plot(ef_vector1, add = TRUE, p. = .99)
(ef_vector2 <- envfit(ordi ~
                        species_richness +
                        accumulated_cover, 
                      data = sites, 
                      permu = 999, 
                      na.rm = T))
plot(ordi, type = "n")
plot(ef_vector2, add = TRUE, p. = .99)

#### b Factors ----------------------------------------------------------------
(ef_factor1 <- envfit(ordi ~  survey_year_factor + exposition + esy, 
              data = sites, permu = 999, na.rm = TRUE))
plot(ordi, type = "n")
ordiellipse(ordi, sites$survey_year_factor, kind = "sd", draw = "lines", label = TRUE)
plot(ordi, type = "n")
ordiellipse(ordi, sites$exposition, kind = "sd", draw = "lines", label = TRUE)
plot(ordi, type = "p")
ordiellipse(ordi, sites$esy, kind = "sd", draw = "lines", label = TRUE)
