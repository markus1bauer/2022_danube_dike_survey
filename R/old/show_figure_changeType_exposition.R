# Model for species richness ####
# Markus Bauer



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Preparation ################################################################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

### Packages ###
library(here)
library(tidyverse)
library(ggbeeswarm)
library(lmerTest)
library(DHARMa)
library(emmeans)

### Start ###
rm(list = ls())
setwd(here("data/processed"))

### Load data ###
sites <- read_csv2("data_processed_sites.csv", col_names = T, na = "na", col_types = 
                     cols(
                       .default = col_guess(),
                       id = col_factor(),
                       location = col_factor(),
                       block = col_factor(),
                       plot = col_factor(),
                       exposition = col_factor(c("north", "south")),
                       ffh = col_factor(),
                       changeType = col_factor(c("FFH6510", "better", "change", "worse", "any-FFH", "non-FFH"))
                     )) %>%
  select(surveyYear, constructionYear, id, block, location, plot, exposition, plotAge, PC1, PC2, PC3, changeType, ffh) %>%
  mutate(plotAge = scale(plotAge, scale = T, center = T)) %>%
  mutate(surveyYear = scale(surveyYear, scale = T, center = T)) %>%
  mutate(surveyYearF = as_factor(surveyYear)) %>%
  mutate(constructionYearF = as_factor(constructionYear)) %>%
  filter(ffh != "6210", 
         changeType != "any-FFH", 
         !is.na(changeType),
         exposition != "east",
         exposition != "west") %>%
  mutate(changeType = fct_collapse(changeType, "change/worse" = c("change", "worse")))




#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Statistics ################################################################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


### 1 Data exploration #####################################################################################

#### a Graphs ---------------------------------------------------------------------------------------------
#2way
ggplot(sites, aes(x = changeType, fill = factor(exposition))) + 
  geom_bar(position = "fill") +
  scale_fill_manual(values = c("royalblue4", "darkorange2"))

### Save ###
setwd(here("outputs/figures"))
ggsave("figure_changeType_exposition_(800dpi_8x7cm).tiff",
       dpi = 800, width = 8, height = 7, units = "cm")

