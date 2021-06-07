# Model for synchrony ####
# Markus Bauer



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Preparation ################################################################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


### Packages ###
library(here)
library(tidyverse)
library(adespatial)
remotes::install_github("larsito/tempo")
library(tempo)

### Start ###
rm(list = ls())
setwd(here("data/processed"))

### Load data ###
sites <- read_csv("data_processed_sites.csv", col_names = T, na = c("na", "NA"), col_types = 
                    cols(
                      .default = "?",
                      id = "f",
                      locationAbb = "f",
                      block = "f",
                      plot = "f",
                      exposition = "f",
                      side = "f",
                      ffh = "f",
                      vegetationCov = "d",
                      locationYear = "f"
                    )) %>%
  select(id, plot, block, locationAbb, surveyYear, constructionYear, plotAge, locationYear, riverkm,
         exposition, side, PC1, PC2, PC3) %>%
  mutate(surveyYearF = as_factor(surveyYear)) %>%
  mutate(constructionYearF = as_factor(constructionYear)) %>%
  column_to_rownames("id")

species <- read_csv("data_processed_species.csv", col_names = T, na = "na", col_types = 
                      cols(
                        .default = "d",
                        name = "f"
                      )) %>%  
  pivot_longer(-name, "id", "value") %>%
  pivot_wider(id, name) %>%
  arrange(id) %>%
  column_to_rownames("id")

data <- species %>%
  rownames_to_column(var = "id") %>%
  mutate(plot = factor(str_sub(id, 1, 3))) %>%
  column_to_rownames(var = "id") %>%
  add_count(plot) %>%
  filter(n > 2) %>%
  select(-n) #%>%
  #group_by(plot) %>%
  #group_split(.keep = F)
data <- split(data, data$plot, drop = T) %>%
  map(~ (.x %>% select(-plot)))
sync_indices <- lapply(data, calc_sync)
data <- as.data.frame(do.call("rbind", sync_indices))

