# Beta diversity on dike grasslands
# Plot Fig. 2 ####
# Markus Bauer
# 2022-09-13



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


### Functions ###
theme_mb <- function() {
  theme(
    panel.background = element_rect(fill = "white"),
    text = element_text(size = 9, color = "black"),
    strip.text = element_text(size = 10),
    axis.text = element_text(angle = 0, hjust = 0.5, size = 9,
                             color = "black"),
    axis.title = element_text(angle = 0, hjust = 0.5, size = 9,
                              color = "black"),
    axis.line = element_line(),
    legend.key = element_rect(fill = "white"),
    legend.position = "right",
    legend.margin = margin(0, 0, 0, 0, "cm"),
    plot.margin = margin(0, 0, 0, 0, "cm")
  )
}

#### * Load data sites ####

sites_dikes <- read_csv("data_processed_sites_spatial.csv",
                        col_names = TRUE, na = c("na", "NA", ""), col_types =
                          cols(
                            .default = "?",
                            id = "f"
                          )) %>%
  select(id, survey_year, orientation, exposition, esy,
         location_construction_year,
         species_richness, eveness, shannon,
         ellenberg_richness, ellenberg_cover_ratio,
         accumulated_cover, graminoid_cover_ratio, ruderal_cover)

sites_splot <- read_csv("data_processed_sites_splot.csv", col_names = TRUE,
                        na = c("na", "NA", ""), col_types =
                          cols(
                            .default = "?"
                          ))

sites <- sites_dikes %>%
  bind_rows(sites_splot) %>%
  mutate(
    reference = as.character(survey_year),
    reference = if_else(
      esy == "E12a", "Dry grassland", if_else(
        esy == "E22", "Hay meadow", reference
      )
    ),
    esy = if_else(
      esy == "E22", "R22-ref", if_else(
        esy == "E12a", "R1A-ref", if_else(
          esy == "?", NA_character_, if_else(
            esy == "R21", "R", esy
          )
        )
      )
    )
  ) %>%
  select(-givd_id, -longitude, -latitude)

#### * Load data species ####

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



### 1 Overview ################################################################


sites %>%
  mutate(esy = replace_na(esy, "no class")) %>%
  filter(!str_detect(esy, "ref")) %>%
  count(esy, survey_year) %>%
  pivot_wider(names_from = "survey_year", values_from = "n") %>%
  write_csv(here("outputs", "statistics",
                 "vegetation_classification.csv"))

sites %>%
  mutate(esy = replace_na(esy, "no class")) %>%
  filter(!str_detect(esy, "ref")) %>%
  count(esy, location_construction_year) %>%
  pivot_wider(names_from = "location_construction_year", values_from = "n")



### 2 NMDS ####################################################################


### Calculate ###
table(sites$esy)
set.seed(1)
(ordi <- metaMDS(species, dist = "bray", binary = TRUE,
                 try = 99, previous.best = TRUE, na.rm = TRUE))
### Stress ###
stressplot(ordi)
goodness_of_fit <- goodness(ordi)
plot(ordi, type = "t", main = "Goodness of fit")
points(ordi, display = "sites", cex = goodness_of_fit * 300)



### 3 Environmental factors ###################################################


#### a Vectors ----------------------------------------------------------------

(ef_vector1 <- envfit(
  ordi ~ species_richness + eveness + shannon +
    accumulated_cover + graminoid_cover_ratio + ruderal_cover +
    ellenberg_richness + ellenberg_cover_ratio,
  data = sites, 
  permu = 999, 
  na.rm = TRUE
  ))
plot(ordi, type = "n")
plot(ef_vector1, add = TRUE, p. = .99)
(ef_vector2 <- envfit(
  ordi ~ ellenberg_richness + ellenberg_cover_ratio +
    graminoid_cover_ratio + ruderal_cover, 
  data = sites, 
  permu = 999, 
  na.rm = TRUE
  ))
plot(ordi, type = "n")
plot(ef_vector2, add = TRUE, p. = .99)


#### b Factors ----------------------------------------------------------------

(ef_factor1 <- envfit(
  ordi ~  orientation + exposition + esy + reference, 
  data = sites, permu = 999, na.rm = TRUE
  ))
plot(ordi, type = "n")
ordiellipse(ordi, sites$orientation, kind = "sd", draw = "lines", label = TRUE)
plot(ordi, type = "n")
ordiellipse(ordi, sites$exposition, kind = "sd", draw = "lines", label = TRUE)
plot(ordi, type = "n")
ordiellipse(ordi, sites$esy, kind = "sd", draw = "lines", label = TRUE)
plot(ordi, type = "n")
ordiellipse(ordi, sites$reference, kind = "sd", draw = "lines", label = TRUE)
