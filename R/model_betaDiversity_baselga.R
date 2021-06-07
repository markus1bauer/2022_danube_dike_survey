# Model for Variation partitioning ####
# Markus Bauer



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Preparation ################################################################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


### Packages ###
library(here)
library(tidyverse)
library(adespatial)

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

sitessoil <- sites %>%
  select(PC1, PC2, PC3, exposition, side)
sitesspacetime <- sites %>%
  select(locationYear, riverkm, surveyYearF)
sitesspace <- sites %>%
  select(locationAbb, riverkm)
sitestime <- sites %>%
  select(surveyYearF, constructionYear)



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Statistics ################################################################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


### 1 Beta diversity #####################################################################################

### Calculate beta diversity: Baselga presence-absence ###
(beta <- beta.div.comp(species, coef = "BS", quant = F))
### Make data frames ###
totalbeta <- beta$D %>% # sörensen dissimilarity
  as.matrix() %>%
  as.data.frame()
substitution <- beta$repl %>% # replacement
  as.matrix()
subsets <- beta$rich %>% # nestedness
  as.matrix()

### Calculate beta diversity: Baselga abundance ###
(beta <- beta.div.comp(species, coef = "BS", quant = T))
### Make data frames ###
totalbeta <- beta$D %>% # Bray-Curtis dissimilarity
  as.matrix() #%>%
  #as.data.frame() %>%
  #rownames_to_column("id") %>%
  #mutate(id = factor(id))
substitution <- beta$repl %>% # balanced variation
  as.matrix() 
subsets <- beta$rich %>% # abundance gradients
  as.matrix()


### 2 RDA #####################################################################################

### a Overall beta diversity +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
m1totalbeta <- dbrda(totalbeta ~ PC1 + PC2 + PC3 + exposition + side +
                 Condition(locationAbb + surveyYearF + plotAge + locationYear + riverkm),
               data = sites)
anova(m1totalbeta)
RsquareAdj(m1totalbeta)

m2totalbeta <- dbrda(totalbeta ~ locationAbb + surveyYearF + plotAge + locationYear + riverkm +
                 Condition(PC1 + PC2 + PC3 + exposition + side),
               data = sites)
anova(m2totalbeta)
RsquareAdj(m2totalbeta)
### Variation partitioning ###
(m3totalbeta <- varpart(totalbeta, sitessoil, sitesspacetime))
plot(m3totalbeta, 
     Xnames = c("Soil", "Space x time"),
     cutoff = 0, digits = 2, bg = NA)
(m4totalbeta <- varpart(totalbeta, sitessoil, sitesspace, sitestime))
plot(m4totalbeta, 
     Xnames = c("Soil", "Space", "Time"),
     cutoff = 0, digits = 1, bg = NA)

### b Substitution +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
m1substitution <- dbrda(substitution ~ PC1 + # .01 p = 0.001
                 Condition(PC2 + PC3 + exposition + side + locationAbb + riverkm + surveyYearF + constructionYearF),
               data = sites)
anova(m1substitution)
RsquareAdj(m1substitution)
m1substitution <- dbrda(substitution ~ PC2 + # .002 p = 0.077
                          Condition(PC1 + PC3 + exposition + side + locationAbb + riverkm + surveyYearF + constructionYearF),
                        data = sites)
anova(m1substitution)
RsquareAdj(m1substitution)
m1substitution <- dbrda(substitution ~ exposition + # .04 p = .001
                          Condition(PC1 + PC2 + PC3 + side + locationAbb + riverkm + surveyYearF + constructionYearF),
                        data = sites)
anova(m1substitution)
RsquareAdj(m1substitution)
m1substitution <- dbrda(substitution ~ side + # .01 p = 0 .001
                          Condition(PC1 + PC2 + PC3 + exposition + locationAbb + riverkm + surveyYearF + constructionYearF),
                        data = sites)
anova(m1substitution)
RsquareAdj(m1substitution)
m1substitution <- dbrda(substitution ~ locationAbb + # .1 p= .001
                          Condition(PC1 + PC2 + PC3 + exposition + side + riverkm + surveyYearF + constructionYearF),
                        data = sites)
anova(m1substitution)
RsquareAdj(m1substitution)
m1substitution <- dbrda(substitution ~ riverkm + # .006 p = .016
                          Condition(PC1 + PC2 + PC3 + exposition + side + locationAbb + surveyYearF + constructionYearF),
                        data = sites)
anova(m1substitution)
RsquareAdj(m1substitution)
m1substitution <- dbrda(substitution ~ surveyYearF + # .003 p = .001
                          Condition(PC1 + PC2 + PC3 + exposition + side + locationAbb + riverkm + constructionYearF),
                        data = sites)
anova(m1substitution)
RsquareAdj(m1substitution)
m1substitution <- dbrda(substitution ~ constructionYearF + # .04 p = .001
                          Condition(PC1 + PC2 + PC3 + exposition + side + locationAbb + riverkm + surveyYear),
                        data = sites)
anova(m1substitution)
RsquareAdj(m1substitution)
m1substitution <- dbrda(substitution ~  locationYear + # .29 p = .001
                          Condition(PC1 + PC2 + PC3 + exposition + side + riverkm + surveyYearF),
                        data = sites)
anova(m1substitution)
RsquareAdj(m1substitution)
### Variation partitioning ###
(m3substitution <- varpart(substitution, sitessoil, sitesspacetime))
plot(m3substitution, 
     Xnames = c("Site", "Space x time"),
     cutoff = 0, digits = 2, bg = NA)
(m4substitution <- varpart(substitution, sitessoil, sitesspace, sitestime))
plot(m4substitution, 
     Xnames = c("Site", "Space", "Time"),
     cutoff = 0, digits = 2, bg = NA)

### c Subsets +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
m1subsets <- dbrda(subsets ~ PC1 + # . p = .
                          Condition(PC2 + PC3 + exposition + side + locationAbb + riverkm + surveyYearF + constructionYearF),
                        data = sites)
anova(m1subsets)
RsquareAdj(m1subsets)
m1subsets <- dbrda(subsets ~ PC2 + # .03 p = 0.061
                          Condition(PC1 + PC3 + exposition + side + locationAbb + riverkm + surveyYearF + constructionYearF),
                        data = sites)
anova(m1subsets)
RsquareAdj(m1subsets)
m1subsets <- dbrda(subsets ~ exposition + # . p = .
                          Condition(PC1 + PC2 + PC3 + side + locationAbb + riverkm + surveyYearF + constructionYearF),
                        data = sites)
anova(m1subsets)
RsquareAdj(m1subsets)
m1subsets <- dbrda(subsets ~ side + # .04 p = 0 .05
                          Condition(PC1 + PC2 + PC3 + exposition + locationAbb + riverkm + surveyYearF + constructionYearF),
                        data = sites)
anova(m1subsets)
RsquareAdj(m1subsets)
m1subsets <- dbrda(subsets ~ locationAbb + # . p= .
                          Condition(PC1 + PC2 + PC3 + exposition + side + riverkm + surveyYearF + constructionYearF),
                        data = sites)
anova(m1subsets)
RsquareAdj(m1subsets)
m1subsets <- dbrda(subsets ~ riverkm + # . p = .
                          Condition(PC1 + PC2 + PC3 + exposition + side + locationAbb + surveyYearF + constructionYearF),
                        data = sites)
anova(m1subsets)
RsquareAdj(m1subsets)
m1subsets <- dbrda(subsets ~ surveyYearF + # .003 p = .001
                          Condition(PC1 + PC2 + PC3 + exposition + side + locationAbb + riverkm + constructionYearF),
                        data = sites)
anova(m1subsets)
RsquareAdj(m1subsets)
m1subsets <- dbrda(subsets ~ constructionYearF + # .12 p = .001
                          Condition(PC1 + PC2 + PC3 + exposition + side + locationAbb + riverkm + surveyYear),
                        data = sites)
anova(m1subsets)
RsquareAdj(m1subsets)
m1subsets <- dbrda(subsets ~  locationYear + # . p = .
                          Condition(PC1 + PC2 + PC3 + exposition + side + riverkm + surveyYearF),
                        data = sites)
anova(m1subsets)
RsquareAdj(m1subsets)
### Variation partitioning ###
(m3subsets <- varpart(subsets, sitessoil, sitesspacetime))
plot(m3subsets, 
     Xnames = c("Site", "Space x time"),
     cutoff = 0, digits = 2, bg = NA)
(m4subsets <- varpart(subsets, sitessoil, sitesspace, sitestime))
plot(m4subsets, 
     Xnames = c("Site", "Space", "Time"),
     cutoff = 0, digits = 1, bg = NA)


