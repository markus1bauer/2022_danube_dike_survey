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
  select(id, plot, block, locationYear, constructionYear,
         exposition, side, PC1soil, PC2soil, PC3soil,
         locationAbb, riverkm,
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

sites <- sites %>%
  column_to_rownames("id")

sites_soil <- sites %>%
  select(PC1soil, PC2soil, PC3soil, expositionN, sideN)
sites_space <- sites %>%
  select(locationAbb, riverkm)
sites_history <- sites %>%
  select(surveyYearF, plotAge, PC1constructionYear, PC2constructionYear, PC3constructionYear)



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Statistics ################################################################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


### 1 Beta diversity #####################################################################################

### * check collinearity ####
data <- sites %>%
  select(where(is.numeric), -ends_with("N"), -accumulatedCov, -constructionYear, -surveyYear)
GGally::ggpairs(data, lower = list(continuous = "smooth_loess"))
#--> no correaltion r > 0.7 (Dormann et al. 2013 Ecography)

### * Calculate: Baselga presence-absence ####
(beta <- beta.div.comp(species, coef = "BS", quant = F)) #total = 0.33, substitution = 0.28, subsets = 0.05
beta_total <- beta$D %>% # sörensen dissimilarity
  as.matrix() %>%
  as.data.frame()
beta_substitution <- beta$repl %>% # replacement / simpson dissimilarity
  as.matrix()
beta_subsets <- beta$rich %>% # nestedness
  as.matrix()


### 2 db-RDA #####################################################################################

### a Overall beta diversity -----------------------------------------------------------------------------------

m1_total_varpart <- varpart(beta_total, sites_soil, sites_space, sites_history)
plot(m1_total_varpart, 
     Xnames = c("Site", "Space", "History"),
     cutoff = 0.01, digits = 1, bg = NA, id.size = 1)


### b Substitution --------------------------------------------------------------------------------------------

### * full model ####
m1 <- dbrda(beta_substitution ~ PC1soil + PC2soil + PC3soil + exposition + side + 
                                       locationAbb + riverkm + 
                                       surveyYearF + plotAge + PC1constructionYear + PC2constructionYear + PC3constructionYear,
            data = sites)
anova(m1, permutations = how(nperm = 999)) #P = .001
(r2adj <- RsquareAdj(m1)$adj.r.squared) #R2adj = .51


### * forward selection ####
m1 <- dbrda(beta_substitution ~ PC1soil + PC2soil + PC3soil + expositionN + sideN, 
            data = sites)
(r2adj <- RsquareAdj(m1)$adj.r.squared)
sel <- forward.sel(beta_substitution, 
                   sites_soil,
                   adjR2thresh = r2adj,
                   nperm = 9999)
sel$p_adj <- p.adjust(sel$pval, method = 'holm', n = ncol(data));sel #https://www.davidzeleny.net/anadat-r/doku.php/en:forward_sel_examples
env.sign <- sort(sel$order)
env.red <- data[ ,c(env.sign)]


m <- dbrda(beta_substitution ~ exposition, data = sites); anova(m); RsquareAdj(m) # p = .001, R2 = .06
m <- dbrda(beta_substitution ~ side, data = sites); anova(m); RsquareAdj(m) # p = .001, R2 = .045
m <- dbrda(beta_substitution ~ locationAbb, data = sites); anova(m); RsquareAdj(m) # p = .001, R2 = .34
m <- dbrda(beta_substitution ~ surveyYearF, data = sites); anova(m); RsquareAdj(m) # p = .001, R2 = .046
# --> locationAbb, exposition, riverkm, surveyYearF, side, PC3cY, plotAge, PC1cY, PC1soil, PC2soil, PC3soil, PC2cY
### step 1 ###
m <- dbrda(beta_substitution ~ locationAbb, data = sites)
anova(m, permutations = how(nperm = 999)) # .001
RsquareAdj(m) #R2adj = .29
### step 2 ###
m <- dbrda(beta_substitution ~ locationAbb + exposition, data = sites)
anova(m, permutations = how(nperm = 999)) # .001
RsquareAdj(m) #R2adj = .34
### step 3 ###
m <- dbrda(beta_substitution ~ locationAbb + exposition + riverkm, data = sites)
anova(m, permutations = how(nperm = 999)) # .001
RsquareAdj(m) #R2adj = .40
### step 4 ###
m <- dbrda(beta_substitution ~ locationAbb + exposition + riverkm + surveyYearF, data = sites)
anova(m, permutations = how(nperm = 999)) # .001
RsquareAdj(m) #R2adj = .36


### * partial db-RDA ####
m1_substitution <- dbrda(beta_substitution ~ PC2soil + # R2adj = .01 p = .001
                          Condition(PC1soil + PC3soil + exposition + side + 
                                      locationAbb + riverkm + 
                                      surveyYearF + plotAge + PC1constructionYear + PC2constructionYear + PC3constructionYear),
                        data = sites)
anova(m1_substitution, permutations = how(nperm = 999))
RsquareAdj(m1_substitution)
m1_substitution <- dbrda(beta_substitution ~ PC3soil + # R2adj = .01 p = .001
                           Condition(PC1soil + PC2soil + exposition + side + 
                                       locationAbb + riverkm + 
                                       surveyYearF + plotAge + PC1constructionYear + PC2constructionYear + PC3constructionYear),
                         data = sites)
anova(m1_substitution, permutations = how(nperm = 999))
RsquareAdj(m1_substitution)
m1_substitution <- dbrda(beta_substitution ~ exposition + #R2adj = .06 p = .001
                          Condition(PC1soil + PC2soil + PC3soil + side + 
                                      locationAbb + riverkm + 
                                      surveyYearF + plotAge + PC1constructionYear + PC2constructionYear + PC3constructionYear),
                        data = sites)
anova(m1_substitution, permutations = how(nperm = 999))
RsquareAdj(m1_substitution)
m1_substitution <- dbrda(beta_substitution ~ side + #R2adj = . p = 0 .001
                          Condition(PC1soil + PC2soil + PC3soil + exposition + 
                                      locationAbb + riverkm + 
                                      surveyYearF + plotAge + PC1constructionYear + PC2constructionYear + PC3constructionYear),
                        data = sites)
anova(m1_substitution, permutations = how(nperm = 999))
RsquareAdj(m1_substitution)
m1_substitution <- dbrda(beta_substitution ~ locationAbb + #R2adj = . p= .001
                          Condition(PC1soil + PC2soil + PC3soil + exposition + side + 
                                      riverkm + 
                                      surveyYearF + plotAge + PC1constructionYear + PC2constructionYear + PC3constructionYear),
                        data = sites)
anova(m1_substitution, permutations = how(nperm = 999))
RsquareAdj(m1_substitution)
m1_substitution <- dbrda(beta_substitution ~ riverkm + #R2adj = . p = .016
                          Condition(PC1soil + PC2soil + PC3soil + exposition + side + 
                                      locationAbb + 
                                      surveyYearF + plotAge + PC1constructionYear + PC2constructionYear + PC3constructionYear),
                        data = sites)
anova(m1_substitution, permutations = how(nperm = 999))
RsquareAdj(m1_substitution)
m1_substitution <- dbrda(beta_substitution ~ surveyYearF + #R2adj = . p = .001
                          Condition(PC1soil + PC2soil + PC3soil + exposition + side + 
                                      locationAbb + riverkm + 
                                      plotAge + PC1constructionYear + PC2constructionYear + PC3constructionYear),
                        data = sites)
anova(m1_substitution, permutations = how(nperm = 999))
RsquareAdj(m1_substitution)
m1_substitution <- dbrda(beta_substitution ~ plotAge + #R2adj = . p = .001
                          Condition(PC1soil + PC2soil + PC3soil + exposition + side + 
                                      locationAbb + riverkm + 
                                      surveyYearF + PC1constructionYear + PC2constructionYear + PC3constructionYear),
                        data = sites)
anova(m1_substitution, permutations = how(nperm = 999))
RsquareAdj(m1_substitution)
m1_substitution <- dbrda(beta_substitution ~  locationYear + #R2adj = . p = .001
                          Condition(PC1soil + PC2soil + PC3soil + exposition + side + 
                                      riverkm + 
                                      surveyYearF),
                        data = sites)
anova(m1_substitution, permutations = how(nperm = 999))
RsquareAdj(m1_substitution)

### * Variation partitioning ####
m1_substitution_varpart <- varpart(beta_substitution, sites_soil, sites_space, sites_history)
tiff(here("outputs/figures/figure_betadiversity_substitution_presence_(800dpi_8x8cm).tiff"),
     res = 72, width = 12, height = 12, units = "cm", compression = "none")
plot(m1_substitution_varpart, 
     Xnames = c("Site", "Space", "History"),
     cutoff = 0.01, digits = 2, bg = NA)
dev.off()


### c Subsets --------------------------------------------------------------------------------------------

### * partial db-RDA ####
m1_subsets <- dbrda(beta_subsets ~ PC1soil + #R2adj = . p = .
                          Condition(PC2soil + PC3soil + exposition + side + locationAbb + riverkm + surveyYearF + constructionYearF),
                        data = sites)
anova(m1_subsets, permutations = how(nperm = 999))
RsquareAdj(m1_subsets)
m1_subsets <- dbrda(beta_subsets ~ PC2soil + #R2adj = . p = 0.061
                          Condition(PC1soil + PC3soil + exposition + side + locationAbb + riverkm + surveyYearF + constructionYearF),
                        data = sites)
anova(m1_subsets, permutations = how(nperm = 999))
RsquareAdj(m1_subsets)
m1_subsets <- dbrda(beta_subsets ~ exposition + #R2adj = . p = .
                          Condition(PC1soil + PC2soil + PC3soil + side + locationAbb + riverkm + surveyYearF + constructionYearF),
                        data = sites)
anova(m1_subsets, permutations = how(nperm = 999))
RsquareAdj(m1_subsets)
m1_subsets <- dbrda(beta_subsets ~ side + #R2adj = . p = 0 .05
                          Condition(PC1soil + PC2soil + PC3soil + exposition + locationAbb + riverkm + surveyYearF + constructionYearF),
                        data = sites)
anova(m1_subsets, permutations = how(nperm = 999))
RsquareAdj(m1_subsets)
m1_subsets <- dbrda(beta_subsets ~ locationAbb + #R2adj = . p= .
                          Condition(PC1soil + PC2soil + PC3soil + exposition + side + riverkm + surveyYearF + constructionYearF),
                        data = sites)
anova(m1_subsets, permutations = how(nperm = 999))
RsquareAdj(m1_subsets)
m1_subsets <- dbrda(beta_subsets ~ riverkm + #R2adj = . p = .
                          Condition(PC1soil + PC2soil + PC3soil + exposition + side + locationAbb + surveyYearF + constructionYearF),
                        data = sites)
anova(m1_subsets, permutations = how(nperm = 999))
RsquareAdj(m1_subsets)
m1_subsets <- dbrda(beta_subsets ~ surveyYearF + #R2adj = . p = .001
                          Condition(PC1soil + PC2soil + PC3soil + exposition + side + locationAbb + riverkm + constructionYearF),
                        data = sites)
anova(m1_subsets, permutations = how(nperm = 999))
RsquareAdj(m1_subsets)
m1_subsets <- dbrda(beta_subsets ~ constructionYearF + #R2adj = .12 p = .001
                          Condition(PC1soil + PC2soil + PC3soil + exposition + side + locationAbb + riverkm + surveyYear),
                        data = sites)
anova(m1_subsets, permutations = how(nperm = 999))
RsquareAdj(m1_subsets)
m1_subsets <- dbrda(beta_subsets ~  locationYear + #R2adj = . p = .
                          Condition(PC1soil + PC2soil + PC3soil + exposition + side + riverkm + surveyYearF),
                        data = sites)
anova(m1_subsets, permutations = how(nperm = 999))
RsquareAdj(m1_subsets)

### * Variation partitioning ####
m1_subsets_varpart <- varpart(beta_subsets, sites_soil, sites_space, sites_history)
tiff(here("outputs/figures/figure_betadiversity_subsets_presence_(800dpi_8x8cm).tiff"),
     res = 72, width = 12, height = 12, units = "cm", compression = "none")
plot(m1_subsets_varpart, 
     Xnames = c("Site", "Space", "History"),
     cutoff = 0.01, digits = 2, bg = NA, id.size = 1)
dev.off()
