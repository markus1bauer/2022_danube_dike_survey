# Beta diversity on dike grasslands
# Non-metric multidimensional scaling (NMDS) ordination ####

# Markus Bauer
# 2024-02-26



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Preparation ################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



### Packages ###
library(here)
library(tidyverse)
library(vegan)

### Start ###
rm(list = ls())

#### * Load data sites ####

sites <- read_csv(
  here("data", "processed", "data_processed_sites_spatial.csv"),
  col_names = TRUE, na = c("na", "NA", ""),
  col_types =
    cols(
      .default = "?",
      id = "f"
    )
) %>%
  select(
    id, survey_year, orientation, exposition, esy,
    location_construction_year,
    species_richness, eveness, shannon,
    ellenberg_richness, ellenberg_cover_ratio,
    accumulated_cover, graminoid_cover_ratio, ruderal_cover
  ) %>%
  filter(esy == "R1A" | esy == "R22" | esy == "R" | esy == "V38" | esy == "?")

#### * Load data species ####

species <- read_csv(
  here("data", "processed", "data_processed_species.csv"),
  col_names = TRUE, na = c("na", "NA", ""),
  col_types =
    cols(
      .default = "d",
      name = "f"
    )
) %>%
  mutate(across(where(is.numeric), ~replace(., is.na(.), 0))) %>%
  pivot_longer(cols = -name, names_to = "id", values_to = "value") %>%
  pivot_wider(names_from = "name", values_from = "value") %>%
  arrange(id) %>%
  semi_join(sites, by = "id") %>%
  column_to_rownames("id")

rm(list = setdiff(ls(), c("sites", "species", "theme_mb")))



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Statistics ################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



### 1 NMDS ####################################################################


### Calculate ###
# set.seed(11)
# ordi <- metaMDS(
#   species, dist = "bray", binary = TRUE,
#   try = 99, previous.best = TRUE, na.rm = TRUE
#   )
base::load(here("outputs", "models", "model_nmds.Rdata"))
ordi

### Stress ###
stressplot(ordi)
goodness_of_fit <- goodness(ordi)
plot(ordi, type = "t", main = "Goodness of fit")
points(ordi, display = "sites", cex = goodness_of_fit * 300)



### 2 Environmental factors ###################################################


#### a Vectors ----------------------------------------------------------------

(ef_vector1 <- envfit(
  ordi ~ species_richness + eveness + shannon +
    accumulated_cover + graminoid_cover_ratio + ruderal_cover +
    ellenberg_richness + ellenberg_cover_ratio + survey_year,
  data = sites,
  permu = 999,
  na.rm = TRUE
  ))
plot(ordi, type = "n")
plot(ef_vector1, add = TRUE, p. = .99)
(ef_vector2 <- envfit(
  ordi ~ ellenberg_richness + graminoid_cover_ratio + ruderal_cover +
    survey_year,
  data = sites,
  permu = 999,
  na.rm = TRUE
  ))
plot(ordi, type = "n")
plot(ef_vector2, add = TRUE, p. = .99)


#### b Factors ----------------------------------------------------------------

(ef_factor1 <- envfit(
  ordi ~  orientation + exposition + esy + as_factor(survey_year),
  data = sites, permu = 999, na.rm = TRUE
  ))
plot(ordi, type = "n")
ordiellipse(ordi, sites$orientation, kind = "sd", draw = "lines", label = TRUE)
plot(ordi, type = "n")
ordiellipse(ordi, sites$exposition, kind = "sd", draw = "lines", label = TRUE)
plot(ordi, type = "n")
ordiellipse(ordi, sites$esy, kind = "sd", draw = "lines", label = TRUE)
plot(ordi, type = "n")
ordiellipse(
  ordi, as_factor(sites$survey_year), kind = "sd", draw = "lines", label = TRUE
)



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# C Save ######################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



save(ordi, file = here("outputs", "models", "model_nmds.Rdata"))
save(
  ef_vector2, file = here("outputs", "models", "model_nmds_envfit_vector.Rdata")
  )
