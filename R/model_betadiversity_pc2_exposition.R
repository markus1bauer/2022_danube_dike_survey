# Model for Variation partitioning ####
# Markus Bauer



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Preparation ################################################################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


### Packages ###
library(here)
library(tidyverse)
library(naniar) #are_na()
library(vegan)
library(adespatial)
library(lme4)

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
  select(id, plot, block, locationYear, constructionYear, longitude, latitude,
         exposition, side, PC1soil, PC2soil, PC3soil,
         locationAbb, riverkm, distanceRiver,
         MEM1, MEM2, MEM3, MEM4, MEM5, MEM6,
         surveyYear, plotAge, PC1constructionYear, PC2constructionYear, PC3constructionYear,
         accumulatedCov) %>%
  mutate(surveyYearF = as_factor(surveyYear),
         expositionN = as.double(exposition),
         sideN = as.double(side),
         locationAbbN = as.double(locationAbb)) %>%
  filter(accumulatedCov > 0)

species <- read_csv("data_processed_species.csv", col_names = T, na = c("na", "NA", ""), col_types = 
                      cols(
                        .default = "d",
                        name = "f"
                      )) %>%  
  mutate(across(where(is.numeric), ~replace(., is.na(.), 0))) %>%
  pivot_longer(-name, "id", "value") %>%
  pivot_wider(id, name) %>%
  semi_join(sites, by = "id") %>%
  arrange(id) %>%
  column_to_rownames("id")



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Statistics ################################################################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


### 1 Beta diversity #####################################################################################

### * check collinearity ####
data <- sites %>%
  select(where(is.numeric), -ends_with("N"), -accumulatedCov, -constructionYear, -surveyYear)
#GGally::ggpairs(data, lower = list(continuous = "smooth_loess"))
#--> MEM1 ~ riverkm has r > 0.7 (Dormann et al. 2013 Ecography) --> MEM1 has to be excluded

### * Calculate: Baselga presence-absence ####
beta <- beta.div.comp(species, coef = "BS", quant = F)
beta$Note
beta$part #total = 0.33, substitution = 0.28, subsets = 0.05
beta_total <- beta$D %>% # sörensen dissimilarity
  as.matrix() %>%
  as.data.frame()
beta_substitution <- beta$repl %>% # replacement / simpson dissimilarity
  as.matrix() %>%
  as.data.frame()
beta_subsets <- beta$rich %>% # nestedness
  as.matrix() %>%
  as.data.frame()

### * Combine with sites ####
beta_total <- beta_total %>%
  select(X39_m_2017) %>%
  rownames_to_column(var = "id") %>%
  rename(beta_total = "X39_m_2017")
beta_substitution <- beta_substitution %>%
  select(X39_m_2017) %>%
  rownames_to_column(var = "id") %>%
  rename(beta_substitution = "X39_m_2017")
beta_subsets <- beta_subsets %>%
  select(X39_m_2017) %>%
  rownames_to_column(var = "id") %>%
  rename(beta_subsets = "X39_m_2017")
data <- sites %>%
  left_join(beta_total, by = "id") %>%
  left_join(beta_substitution, by = "id") %>%
  left_join(beta_subsets, by = "id") %>%
  filter(exposition == "south" | exposition == "north")

### Calculate ####
ggplot(data, aes(y = beta_total, x = PC1soil, color = exposition)) +
  #geom_point() +
  geom_smooth(method = "loess") +
  facet_wrap(~surveyYearF)
ggplot(data, aes(y = beta_substitution, x = PC1soil, color = exposition)) +
  #geom_point() +
  geom_smooth(method = "loess") +
  facet_wrap(~surveyYearF)
ggplot(data, aes(y = beta_subsets, x = PC1soil, color = exposition)) +
  #geom_point() +
  geom_smooth(method = "loess") +
  facet_wrap(~surveyYearF)

### riverkm ####
beta_total <- beta_total %>%
  select(X55_m_2017) %>%
  rownames_to_column(var = "id") %>%
  rename(beta_total = "X55_m_2017")
beta_substitution <- beta_substitution %>%
  select(X55_m_2017) %>%
  rownames_to_column(var = "id") %>%
  rename(beta_substitution = "X55_m_2017")
beta_subsets <- beta_subsets %>%
  select(X55_m_2017) %>%
  rownames_to_column(var = "id") %>%
  rename(beta_subsets = "X55_m_2017")
data <- sites %>%
  left_join(beta_total, by = "id") %>%
  left_join(beta_substitution, by = "id") %>%
  left_join(beta_subsets, by = "id") %>%
  filter(exposition == "south" | exposition == "north")

ggplot(data, aes(y = beta_total, x = riverkm, color = exposition)) +
  #geom_point() +
  geom_smooth(method = "loess") +
  facet_wrap(~surveyYearF)
ggplot(data, aes(y = beta_substitution, x = riverkm, color = exposition)) +
  #geom_point() +
  geom_smooth(method = "loess") +
  facet_wrap(~surveyYearF)
ggplot(data, aes(y = beta_subsets, x = riverkm, color = exposition)) +
  #geom_point() +
  geom_smooth(method = "loess") +
  facet_wrap(~surveyYearF)
