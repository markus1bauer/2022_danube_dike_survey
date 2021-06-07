# Model for species richness ####
# Markus Bauer



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Preparation ################################################################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

### Packages ###
library(here)
library(tidyverse)
library(ggbeeswarm)
library(ade4)
library(DHARMa)
library(emmeans)

### Start ###
rm(list = ls())
setwd(here("data/processed"))

### Load data ###
sites <- read_csv2("data_processed_sites.csv", col_names = T, na = c("na", "NA"), col_types = 
                     cols(
                       .default = col_guess(),
                       id = "f",
                       locationAbb = "f",
                       block = "f",
                       plot = "f",
                       exposition = "f",
                       side = "f",
                       ffh = "f",
                       changeType = col_factor(c("FFH6510", "any-FFH", "better", "change", "worse", "non-FFH"))
                     )) %>%
  select(targetRichness, id, surveyYear, constructionYear, plotAge, locationAbb, block, plot, side, exposition, PC1, PC2, PC3, changeType, ffh) %>%
  mutate(n = targetRichness) %>%
  mutate(location = factor(locationAbb, levels = unique(locationAbb[order(constructionYear)]))) %>%
  mutate(plotAge = scale(plotAge, scale = T, center = F)) %>%
  mutate(surveyYear = scale(surveyYear, scale = T, center = F)) %>%
  mutate(surveyYearF = as_factor(surveyYear)) %>%
  mutate(constructionYear = scale(constructionYear, scale = T, center = F)) %>%
  mutate(constructionYearF = as_factor(constructionYear)) %>%
  filter(exposition != "east",
         exposition != "west") %>%
  group_by(plot) %>%
  mutate(count = n()) %>%
  ungroup() %>%
  slice_max(count, n = 1) %>%
  mutate(exposition = factor(exposition)) %>%
  mutate(locationAbb = factor(locationAbb)) %>%
  mutate(plot = factor(plot)) %>%
  mutate(block = factor(block)) %>%
  select(-count) %>%
  select(id, exposition, side, PC1, PC2, PC3)

traits <- read_csv2("data_processed_traits.csv", col_names = T, na = c("na", "NA"), col_types = 
                       cols(
                         .default = "?",
                         name = "f",
                         sla = "d",
                         seedmass = "d",
                         height = "d"
                       )) %>%  
  select(name, sla, seedmass, height) %>%
  drop_na()

species <- read_csv2("data_processed_species.csv", col_names = T, na = c("na", "NA"), col_types = 
                       cols(
                         .default = "d",
                         name = "f"
                       )) %>%
  semi_join(traits, by = "name") %>%
  pivot_longer(-name, "id", "value") %>%
  pivot_wider(id, name) %>%
  arrange(id) %>%
  semi_join(sites, by = "id") %>%
  column_to_rownames("id")

sites <- sites %>% column_to_rownames("id")
traits <- traits %>% column_to_rownames("name")



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Statistics ################################################################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


## 1 RLQ #############################################################################################

### Preliminary analyses; CA, Hill-Smith and PCA) ###
L <- dudi.coa(species, scannf = F)
R <- dudi.hillsmith(sites, scannf = F, row.w = L$lw)
Q <- dudi.pca(traits, scannf = F, row.w = L$cw)
### RLQ analyses ###
rlq <- rlq(R, L, Q, scannf = F)
plot(rlq)
### Site (L) scores ###
s.label(rlq$lR, boxes = T)
### Species (Q) scores ###
s.label(rlq$lQ, boxes = F)
# Environmental variables:
s.arrow(rlq$l1, xlim = c(-1.4,2))
# Species traits:
s.arrow(rlq$c1, xlim = c(-1.4,0.6))
# Global test:
randtest(rlq, nrepet = 999, 
         modeltype = 6)


## 2 Fourth-corner analysis #############################################################################################

m1 <- fourthcorner(sites, species, traits,
                   modeltype = 6, nrep = 49999, p.adjust.method.D = "none", p.adjust.method.G = "none");
m1.adj <-  p.adjust.4thcorner(m1, 
                              p.adjust.method.G = "fdr", 
                              p.adjust.method.D = "fdr",
                              p.adjust.D = "global")
summary(m1.adj)
plot(m1.adj, alpha = 0.05, stat = "D2") #nothing significant
plot(m1.adj, alpha = 0.05, stat = "D2", x.rlq = rlq, type = "biplot")
