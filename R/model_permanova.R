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
sites <- read_csv2("data_processed_sites.csv", col_names = T, na = "na", col_types = 
                     cols(
                       .default = col_guess(),
                       id = col_factor(),
                       location = col_factor(),
                       block = col_factor(),
                       plot = col_factor(),
                       exposition = col_factor(),
                       ffh = col_factor(),
                       changeType = col_factor(),
                       vegetationCov = col_double()
                     )) %>%
  select(id, plot, block, location, surveyYear, constructionYear, plotAge, exposition, PC1, PC2, PC3, ffh, changeType, vegetationCov, targetRichness, targetCov, ffh6510Richness, ffh6210Richness) %>%
  mutate(surveyYearF = as_factor(surveyYear))

species <- read_csv2("data_processed_species.csv", col_names = T, na = "na", col_types = 
                       cols(
                         .default = col_double(),
                         name = col_factor()
                       )) %>%  
  pivot_longer(-name, "id", "value") %>%
  pivot_wider(id, name) %>%
  arrange(id) %>%
  column_to_rownames("id")

sitesFFH <- sites %>%
  filter(ffh != "6210")
speciesFFH <- species %>%
  rownames_to_column("id") %>%
  semi_join(sitesFFH, by = "id") %>%
  column_to_rownames("id")



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Statistics ################################################################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


### 1 PERMDISP #####################################################################################

(permdisp <- betadisper(d = vegdist(species), group = sites$exposition))
permutest(permdisp, pairwise = T) # similar
(permdisp <- betadisper(d = vegdist(speciesFFH), group = sitesFFH$ffh))
permutest(permdisp, pairwise = T) # not similar
(permdisp <- betadisper(d = vegdist(speciesFFH), group = sitesFFH$surveyYearF))
permutest(permdisp, pairwise = T) # similar
(permdisp <- betadisper(d = vegdist(speciesFFH), group = sitesFFH$location))
permutest(permdisp, pairwise = T) # not similar


### 2 PERMANOVA #####################################################################################

(permanova <- adonis(speciesFFH ~ plotAge + constructionYear + surveyYear * location + 
                       exposition + PC1 + PC2 + PC3 + 
                       changeType, 
                     data = sitesFFH, 
                     strata = sitesFFH$plot, 
                     permutations = 999, 
                     method = "bray"))
densityplot(permustats(permanova)) 
table <- fortify(permanova$aov.tab) %>%
  rownames_to_column(var = "variable")

### Save ###
setwd(here("data/tables"))
write_csv2(table, "table_permanova.csv")
