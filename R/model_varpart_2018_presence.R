# Beta diversity on dike grasslands
# Variation partitioning of 2018 (presence-absence data) ####
# Markus Bauer
# 2022-01-11



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Preparation #########################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


### Packages ###
library(here)
library(tidyverse)
library(vegan)
library(adespatial)

### Start ###
rm(list = ls())
setwd(here("data", "processed"))

### Load data ###
sites <- read_csv("data_processed_sites_spatial.csv",
  col_names = TRUE,
  na = c("na", "NA"), col_types =
    cols(
      .default = "?",
      id = "f",
      locationAbb = "f",
      block = "f",
      plot = "f",
      exposition = "f",
      side = "f",
      locationYear = "f"
    )
) %>%
  select(
    id, plot, block, locationYear, constructionYear, longitude, latitude,
    exposition, side, PC1soil, PC2soil, PC3soil,
    locationAbb, riverkm, distanceRiver, MEM2_2018,
    surveyYear, plotAge, PC1constructionYear, PC2constructionYear,
    PC3constructionYear,
    accumulatedCov
  ) %>%
  filter(surveyYear == 2018) %>%
  mutate(
    surveyYearF = as_factor(surveyYear),
    exposition_numeric = as.double(exposition),
    side_numeric = as.double(side),
    locationAbb_numeric = as.double(locationAbb)
  )

species <- read_csv("data_processed_species.csv",
  col_names = TRUE,
  na = c("na", "NA", ""), col_types =
    cols(
      .default = "d",
      name = "f"
    )
) %>%
  mutate(across(where(is.numeric), ~ replace(., is.na(.), 0))) %>%
  pivot_longer(-name, "id", "value") %>%
  pivot_wider(id, name) %>%
  semi_join(sites, by = "id") %>%
  arrange(id) %>%
  column_to_rownames("id")

sites <- sites %>%
  column_to_rownames("id")

sites_soil <- sites %>%
  select(PC1soil, PC2soil, PC3soil, exposition_numeric, side_numeric)
sites_space <- sites %>%
  select(locationAbb_numeric, distanceRiver, riverkm, MEM2_2018)
sites_history <- sites %>%
  select(plotAge, PC1constructionYear, PC2constructionYear)



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Statistics ##########################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


## 1 Beta diversity #####################################################

### * check collinearity ####
data <- sites %>% select(
  where(is.numeric), -ends_with("numeric"),
  -accumulatedCov, -constructionYear, -surveyYear
)
GGally::ggpairs(data, lower = list(continuous = "smooth_loess"))
#--> distanceRiver ~ PC3constructionYear r > 0.7 (Dormann et al. 2013 Ecography)
#--> PC3constructionYear removed

### * Calculate: Baselga presence-absence ####
beta <- beta.div.comp(species, coef = "BS", quant = FALSE)
beta$Note
beta$part # total = 0.331, substitution = 0.281, subsets = 0.050
beta_total <- beta$D %>% # sörensen dissimilarity
  as.matrix() %>%
  as.data.frame()
beta_substitution <- beta$repl %>% # replacement / simpson dissimilarity
  as.matrix()
beta_subsets <- beta$rich %>% # nestedness
  as.matrix()


## 2 db-RDA #############################################################

### a Overall variation partitioning ------------------------------------

(m1_total_varpart <- varpart(beta_total, sites_soil, sites_space, sites_history))
plot(
  m1_total_varpart,
  Xnames = c("Site", "Space", "History"),
  cutoff = 0.01, digits = 1, bg = NA, id.size = 1
  )

### b Substitution ------------------------------------------------------

### * linear trend in data ####
# m1 <- dbrda(beta_substitution ~ longitude + latitude, data = sites)
# anova(m1) #sig
# beta_substitution_detrended <- resid(lm(beta_substitution ~ longitude + latitude, data = sites))
#--> this trend is captured by riverkm

### * full model ####
m1 <- dbrda(beta_substitution ~ PC1soil + PC2soil + PC3soil + exposition + side +
  locationAbb + riverkm + distanceRiver + MEM2_2018 +
  plotAge + PC1constructionYear + PC2constructionYear,
data = sites
)
anova(m1, permutations = how(nperm = 9999)) # P = 1e-04***
(r2adj <- RsquareAdj(m1)$adj.r.squared) # R2adj = .307

### * forward selection ####
### Soil ###
m1 <- dbrda(beta_substitution ~ PC1soil + PC2soil + PC3soil + exposition + side,
  data = sites
)
r2adj <- RsquareAdj(m1)$adj.r.squared
sel <- forward.sel(beta_substitution,
  sites_soil,
  adjR2thresh = r2adj,
  nperm = 9999
)
sel$p_adj <- p.adjust(sel$pvalue, method = "holm", n = ncol(sites_soil))
sel # https://www.davidzeleny.net/anadat-r/doku.php/en:forward_sel_examples
sites_soil_selected <- sites %>%
  select(PC3soil)
### Space ###
m1 <- dbrda(beta_substitution ~ MEM2_2018 + locationAbb + riverkm + distanceRiver,
  data = sites
)
r2adj <- RsquareAdj(m1)$adj.r.squared
sel <- forward.sel(beta_substitution,
  sites_space,
  adjR2thresh = r2adj,
  nperm = 9999
)
sel$p_adj <- p.adjust(sel$pvalue, method = "holm", n = ncol(sites_space))
sel # https://www.davidzeleny.net/anadat-r/doku.php/en:forward_sel_examples
sites_space_selected <- sites %>%
  select(locationAbb_numeric)
### History ###
m1 <- dbrda(beta_substitution ~ plotAge + PC1constructionYear + PC2constructionYear,
  data = sites
)
(r2adj <- RsquareAdj(m1)$adj.r.squared)
sel <- forward.sel(beta_substitution,
  sites_history,
  adjR2thresh = r2adj,
  nperm = 9999
)
sel$p_adj <- p.adjust(sel$pvalue, method = "holm", n = ncol(sites_history))
sel # https://www.davidzeleny.net/anadat-r/doku.php/en:forward_sel_examples
sites_history_selected <- sites %>%
  select()

### * Variation partitioning ####
(m1_substitution_varpart <- varpart(beta_substitution, sites_soil_selected, sites_space_selected))
tiff(
  here("outputs", "figures", "figure_4b_2018_800dpi_12x12cm.tiff"),
  res = 72, width = 12, height = 12, units = "cm", compression = "none"
  )
plot(
  m1_substitution_varpart,
  Xnames = c("Site", "Space"),
  cutoff = 0.01, digits = 2, bg = NA
  )
dev.off()

### * partial db-RDA ####
### Soil / PC3soil ###
m1_substitution <- dbrda(beta_substitution ~ PC3soil +
  Condition(locationAbb),
data = sites
)
anova(m1_substitution, permutations = how(nperm = 9999)) # p = .579
RsquareAdj(m1_substitution) # R2adj = .017
### Space / locationAbb ###
m1_substitution <- dbrda(beta_substitution ~ locationAbb +
  Condition(PC3soil),
data = sites
)
anova(m1_substitution, permutations = how(nperm = 9999)) # p = 1e-04
RsquareAdj(m1_substitution) # R2adj = .179

### c Subsets ----------------------------------------------------------

### * full model ####
m1 <- dbrda(beta_subsets ~ PC1soil + PC2soil + PC3soil + exposition + side +
  locationAbb + riverkm + distanceRiver +
  plotAge + PC1constructionYear + PC2constructionYear,
data = sites
)
anova(m1, permutations = how(nperm = 9999)) # P = .523
(r2adj <- RsquareAdj(m1)$adj.r.squared) # R2adj = -.015

### * forward selection ####
# --> now forward selection because full model is not significant
