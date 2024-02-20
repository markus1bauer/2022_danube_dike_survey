# Beta diversity on dike grasslands
# Habitat types of restoration outcomes ####

# Markus Bauer
# 2024-02-20



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Preparation ################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



### Packages ###
library(here)
library(tidyverse)
library(vegan)

### Start ###
rm(list = ls())

### Load data ###

sites_dikes <- read_csv(
  here("data", "processed", "data_processed_sites_spatial.csv"),
  col_names = TRUE, na = c("na", "NA", ""),
  col_types =
    cols(
      .default = "?",
      id = "f"
    )
) %>%
  select(id, survey_year, orientation, exposition, esy)

sites_splot <- read_csv(
  here("data", "processed", "data_processed_sites_splot.csv"),
  col_names = TRUE, na = c("na", "NA", ""),
  col_types =
    cols(
      .default = "?"
    )
)

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



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Statistics ################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



### 1 Data exploration #########################################################


### a Table of raw data -------------------------------------------------------

table(sites$esy)
data <- sites %>%
  mutate(esy = replace_na(esy, "no class")) %>%
  filter(!str_detect(esy, "ref")) %>%
  count(esy, survey_year, exposition) %>%
  mutate(ratio = round(n/41, digits = 2))
data_table <- data %>%
  select(-n) %>%
  pivot_wider(names_from = "survey_year", values_from = "ratio") %>%
  arrange(esy, exposition)
data_table

### b Graphs of raw data -------------------------------------------------------

data %>%
  ggplot(aes(y = ratio, x = survey_year)) +
  facet_grid(vars(esy), vars(exposition)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(limits = c(0, 0.4))


### 2 Save ####################################################################

data_table %>%
  write_csv(here("outputs", "statistics", "vegetation_classification.csv"))


