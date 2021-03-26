# Model for ordination ####
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
                      .default = col_double(),
                      id = col_factor(),
                      location = col_factor(),
                      block = col_factor(),
                      plot = col_factor(),
                      position = col_factor(),
                      dataset = col_factor(),
                      side = col_factor(),
                      exposition = col_factor(),
                      phosphorousClass = col_factor(c("A", "B", "C", "D", "E")),
                      biotopeType = col_factor(),
                      baykompv = col_factor(),
                      ffh = col_factor(),
                      min8 = col_factor(),
                      min9 = col_factor()
                    )        
)

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


### 1 NMDS #####################################################################################

#### a ordination ----------------------------------------------------------------------------------------
set.seed(1)
(ordi <- metaMDS(species, try = 99, previous.best = T, na.rm = T))
stressplot(ordi)
# stress: 0.286 (wisconsin(sqrt()))
#R²linear = .593
#R²non-metric = .918

#### b environmental factors ----------------------------------------------------------------------------------------
### vectors ###
(ef_vector1 <- envfit(ordi ~  vegetationCov + plotAge + topsoilDepth + phosphorous + pH + calciumcarbonatPerc + humusPerc + NtotalPerc + cnRatio + sandPerc + siltPerc + clayPerc + NtotalConc + ufc + targetRichness + ffh6510Richness + ffh6210Richness + targetCov, 
              data = sites, permu = 999, na.rm = T))
plot(ordi, type = "n"); plot(ef_vector1, add = T, p. = .99)
(ef_vector2 <- envfit(ordi ~  vegetationCov + plotAge + topsoilDepth + phosphorous + NtotalConc + NtotalPerc + cnRatio + sandPerc + ffh6510Richness + ffh6210Richness, 
                      data = sites, permu = 999, na.rm = T))
plot(ordi, type = "n"); plot(ef_vector2, add = T, p. = .99)
### factors ###
(ef_factor1 <- envfit(ordi ~  plot + block + side + surveyYear + exposition + location + ffh + biotopeType + biotopePoints + baykompv + phosphorousClass, 
              data = sites, permu = 999, na.rm = T))
plot(ordi, type = "n"); ordiellipse(ordi, sites$location, kind = "sd", draw = "lines", label = T)
plot(ordi, type = "n"); ordiellipse(ordi, sites$surveyYear, kind = "sd", draw = "lines", label = T)
plot(ordi, type = "n"); ordiellipse(ordi, sites$exposition, kind = "sd", draw = "lines", label = T)
plot(ordi, type = "n"); ordiellipse(ordi, sites$ffh, kind = "sd", draw = "lines", label = T)
plot(ordi, type = "n"); ordiellipse(ordi, sites$plot, kind = "sd", draw = "lines", label = T)


### 2 PERMANOVA #####################################################################################

#### a exposition ----------------------------------------------------------------------------------------
(permdisp <- betadisper(d = vegdist(species), group = sites$exposition))
permutest(permdisp, pairwise = T) # similar
(permanova <- adonis(species ~ exposition, data = sites, 
                          strata = sites$plot, permutations = 999, method = "bray"))
densityplot(permustats(permanova)) # R² = .04, n.s.

#### b ffh ----------------------------------------------------------------------------------------
(permdisp <- betadisper(d = vegdist(speciesFFH), group = sitesFFH$ffh))
permutest(permdisp, pairwise = T) # not similar
(permanova <- adonis(speciesFFH ~ ffh, data = sitesFFH, 
                     strata = sitesFFH$plot, permutations = 999, method = "bray"))
densityplot(permustats(permanova)) # p= .004, R² = .03
