# Model for PERMANOVA ####
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
sites <- read_csv("data_processed_sites.csv", col_names = T, na = c("na", "NA", ""), col_types = 
                    cols(
                      .default = "d",
                      id = "f",
                      locationAbb = "f",
                      block = "f",
                      plot = "f",
                      exposition = "f"
                    )) %>%
  select(id, plot, block, locationAbb, surveyYear, exposition, vegetationCov, targetCov, graminoidCovratio, speciesRichness, targetRichness, shannon, eveness, accumulatedCov) %>%
  mutate(surveyYearF = as_factor(surveyYear)) %>%
  filter(accumulatedCov > 0)

species <- read_csv("data_processed_species.csv", col_names = T, na = c("na", "NA", ""), col_types = 
                      cols(
                        .default = "d",
                        name = "f"
                      )) %>%  
  mutate(across(where(is.numeric), ~replace(., is.na(.), 0))) %>%
  pivot_longer(-name, "id", "value") %>%
  pivot_wider(id, name) %>%
  arrange(id) %>%
  semi_join(sites, by = "id") %>%
  column_to_rownames("id")



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Statistics ################################################################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


### 1 PERMDISP #####################################################################################

permdisp <- betadisper(d = vegdist(species), group = sites$exposition)
permutest(permdisp, pairwise = T) # south is different from other directions
permdisp <- betadisper(d = vegdist(species), group = sites$surveyYearF)
permutest(permdisp, pairwise = T) # different


### 2 PERMANOVA #####################################################################################

(permanova <- adonis(species ~ exposition + surveyYearF,
                     data = sites, 
                     strata = sites$plot, 
                     permutations = 999, 
                     method = "bray"))
densityplot(permustats(permanova))
pairwiseAdonis::pairwise.adonis2(species ~ exposition + surveyYearF,
                 data = sites)
table <- fortify(permanova$aov.tab) %>%
  rownames_to_column(var = "variable")

### Save ###
write_csv(table, here("outputs/tables/table_permanova.csv"))
