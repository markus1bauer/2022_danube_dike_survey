# Beta diversity on dike grasslands
# Plot Fig. 2 ####
# Markus Bauer
# 2022-09-05



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Preparation ################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



### Start ###
rm(list = ls())
setwd(here("data", "processed"))

### Packages ###
library(here)
library(tidyverse)
library(vegan)

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

veganCovEllipse <- function(cov, center = c(0, 0), scale = 1, npoints = 100) {
  theta <- (0:npoints) * 2 * pi / npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

#### * Load data sites ####

sites_dikes <- read_csv("data_processed_sites_spatial.csv",
                        col_names = TRUE, na = c("na", "NA", ""), col_types =
                          cols(
                            .default = "?",
                            id = "f"
                          )) %>%
  select(id, survey_year, orientation, exposition,
         species_richness, eveness, shannon,
         target_richness, target_cover_ratio,
         accumulated_cover, graminoid_cover_ratio, ruderal_cover) %>%
  mutate(survey_year_factor = as_factor(survey_year),
         target_richness_group = if_else(
           target_richness < 10, "<10", if_else(
             target_richness >= 10 & target_richness < 20, "10-19", if_else(
               target_richness >= 20 & target_richness < 30, "20-29", if_else(
                 target_richness >= 30, ">30", "warning"
               )
             )
           )
         ),
         target_richness_group = fct_relevel(target_richness_group,
                                             "<10", "10-19", "20-29", ">30"))

sites_splot <- read_csv("data_processed_sites_splot.csv", col_names = TRUE,
                        na = c("na", "NA", ""), col_types =
                          cols(
                            .default = "?"
                          ))

sites <- sites_dikes %>%
  bind_rows(sites_splot) %>%
  mutate(
    survey_year = if_else(is.na(survey_year), 0, survey_year),
    esy = if_else(
      survey_year == 2017, "2017", if_else(
        survey_year == 2018, "2018", if_else(
          survey_year == 2019, "2019", if_else(
            survey_year == 2021, "2021", if_else(
              esy == "E12a", "Dry grassland", if_else(
                esy == "E22", "Hay meadow", "warning"
              )
            )
          )
        )
      )
    )
  )

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



### 1 NMDS ####################################################################


### Calculate ###
set.seed(1)
(ordi <- metaMDS(species, dist = "bray", binary = TRUE,
                 try = 99, previous.best = TRUE, na.rm = TRUE))
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
    target_richness + target_cover_ratio,
  data = sites, 
  permu = 999, 
  na.rm = TRUE
  ))
plot(ordi, type = "n")
plot(ef_vector1, add = TRUE, p. = .99)
(ef_vector2 <- envfit(
  ordi ~ target_richness + graminoid_cover_ratio + ruderal_cover, 
  data = sites, 
  permu = 999, 
  na.rm = TRUE
  ))
plot(ordi, type = "n")
plot(ef_vector2, add = TRUE, p. = .99)


#### b Factors ----------------------------------------------------------------

(ef_factor1 <- envfit(
  ordi ~  survey_year_factor + orientation + exposition + esy, 
  data = sites, permu = 999, na.rm = TRUE
))
plot(ordi, type = "n")
ordiellipse(ordi, sites$survey_year_factor, kind = "sd", draw = "lines",
            label = TRUE)
plot(ordi, type = "n")
ordiellipse(ordi, sites$orientation, kind = "sd", draw = "lines", label = TRUE)
plot(ordi, type = "n")
ordiellipse(ordi, sites$exposition, kind = "sd", draw = "lines", label = TRUE)
plot(ordi, type = "n")
ordiellipse(ordi, sites$esy, kind = "sd", draw = "lines", label = TRUE)
