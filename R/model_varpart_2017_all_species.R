# Beta diversity on dike grasslands
# Variation partitioning of 2017 - All species ####
# Markus Bauer
# 2022-09-15



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Preparation ###############################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



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
      location_abb = "f",
      block = "f",
      plot = "f",
      exposition = "f",
      orientation = "f",
      location_construction_year = "f"
    )) %>%
  filter(survey_year == 2017) %>%
  select(
    id, plot, block, longitude, latitude,
    botanist, location_construction_year, construction_year,
    exposition, orientation, pc1_soil, pc2_soil, pc3_soil,
    location_abb, river_km, river_distance,
    survey_year, plot_age, pc1_construction_year, pc2_construction_year,
    pc3_construction_year,
    accumulated_cover
    ) %>%
  mutate(
    survey_year_factor = as_factor(survey_year),
    exposition_numeric = as.double(exposition),
    orientation_numeric = as.double(orientation),
    location_abb_numeric = as.double(location_abb),
    botanist_numeric = as.double(as_factor(botanist))
  )

species <- read_csv("data_processed_species.csv",
  col_names = TRUE,
  na = c("na", "NA", ""), col_types =
    cols(
      .default = "d",
      name = "f"
    )) %>%
  mutate(across(where(is.numeric), ~ replace(., is.na(.), 0))) %>%
  pivot_longer(-name, names_to = "id", values_to = "value") %>%
  pivot_wider(id, names_from = "name", values_from = "value") %>%
  semi_join(sites, by = "id") %>%
  arrange(id) %>%
  column_to_rownames("id")

sites <- sites %>%
  column_to_rownames("id")

sites_soil <- sites %>%
  select(pc1_soil, pc2_soil, pc3_soil, exposition_numeric, orientation_numeric)
sites_space <- sites %>%
  select(location_abb_numeric, river_distance, river_km)
sites_history <- sites %>%
  select(
    plot_age, pc1_construction_year, pc2_construction_year,
    pc3_construction_year
    )



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Statistics ################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



## 1 Calculate beta diversity ##################################################


### * Check collinearity ####
data <- sites %>% select(
  where(is.numeric), -ends_with("numeric"),
  -accumulated_cover, -construction_year, -survey_year, -longitude, -latitude
  )
GGally::ggpairs(data, lower = list(continuous = "smooth_loess"))
#--> no relationship has r > 0.7 (Dormann et al. 2013 Ecography)

### * Calculate: Baselga presence-absence ####
beta <- beta.div.comp(species, coef = "BS", quant = FALSE)
beta$Note
beta$part # total = 0.317, substitution = 0.275, subsets = 0.041
beta_total <- beta$D %>% # sörensen dissimilarity
  as.matrix() %>%
  as.data.frame()
beta_substitution <- beta$repl %>% # replacement / simpson dissimilarity
  as.matrix()
beta_subsets <- beta$rich %>% # nestedness
  as.matrix()



## 2 db-RDA ###################################################################


### a Overall variation partitioning ------------------------------------------

m1_total_varpart <- varpart(beta_total, sites_soil, sites_space, sites_history)
plot(
  m1_total_varpart,
  Xnames = c("Site", "Space", "History"),
  cutoff = 0.01, digits = 1, bg = NA, id.size = 1
  )


### b Substitution ------------------------------------------------------------

### * Check linear trend in data ####
# m1 <- dbrda(beta_substitution ~ longitude + latitude, data = sites)
# anova(m1) #sig
# beta_substitution_detrended <- resid(lm(beta_substitution ~ longitude + latitude, data = sites))
#--> this trend is captured by river_km
### * full model ####
m1 <- dbrda(
  beta_substitution ~
    pc1_soil + pc2_soil + pc3_soil + exposition + orientation +
    location_abb + river_km + river_distance +
    plot_age + pc1_construction_year + pc2_construction_year +
    pc3_construction_year,
  data = sites
  )
anova(m1, permutations = how(nperm = 9999)) # P: .001
(r2adj <- RsquareAdj(m1)$adj.r.squared) # R2adj: .501

### * Forward selection soil ####
m1 <- dbrda(
  beta_substitution ~ pc1_soil + pc2_soil + pc3_soil + exposition + orientation,
  data = sites
  )
r2adj <- RsquareAdj(m1)$adj.r.squared
sel <- forward.sel(
  beta_substitution,
  sites_soil,
  adjR2thresh = r2adj,
  nperm = 9999
  )
sel$p_adj <- p.adjust(sel$pvalue, method = "holm", n = ncol(sites_soil))
sel # https://www.davidzeleny.net/anadat-r/doku.php/en:forward_sel_examples
sites_soil_selected <- sites %>%
  select(pc3_soil, orientation_numeric, exposition_numeric, pc1_soil)

### * Forward selection space ####
m1 <- dbrda(
  beta_substitution ~ location_abb + river_km + river_distance,
  data = sites
  )
r2adj <- RsquareAdj(m1)$adj.r.squared
sel <- forward.sel(
  beta_substitution,
  sites_space,
  adjR2thresh = r2adj,
  nperm = 9999
  )
sel$p_adj <- p.adjust(sel$pvalue, method = "holm", n = ncol(sites_space))
sel # https://www.davidzeleny.net/anadat-r/doku.php/en:forward_sel_examples
sites_space_selected <- sites %>%
  select(location_abb_numeric)

### * Forward selection history ####
m1 <- dbrda(
  beta_substitution ~ plot_age + pc1_construction_year + pc2_construction_year +
  pc3_construction_year,
  data = sites
  )
r2adj <- RsquareAdj(m1)$adj.r.squared
sel <- forward.sel(
  beta_substitution,
  sites_history,
  adjR2thresh = r2adj,
  nperm = 9999
  )
sel$p_adj <- p.adjust(sel$pvalue, method = "holm", n = ncol(sites_history))
sel # https://www.davidzeleny.net/anadat-r/doku.php/en:forward_sel_examples
sites_history_selected <- sites %>%
  select(pc1_construction_year)

### * Variation partitioning ####
m1_substitution_varpart <- varpart(
  beta_substitution, sites_soil_selected,
  sites_space_selected, sites_history_selected
)
tiff(
  here("outputs", "figures", "figure_3a_2017_800dpi_12x12cm.tiff"),
  res = 72, width = 12, height = 12, units = "cm", compression = "none"
  )
plot(
  m1_substitution_varpart,
  Xnames = c("Site", "Space", "History"),
  cutoff = 0.01, digits = 2, bg = NA
  )
dev.off()

### * Partial db-RDA soil ####
m1_substitution <- dbrda(
  beta_substitution ~ pc3_soil + orientation + exposition + pc1_soil +
    Condition(location_abb + pc1_construction_year),
data = sites
)
anova(m1_substitution, permutations = how(nperm = 9999)) # p: 2e-04 ***
RsquareAdj(m1_substitution) # R2adj: .141

### * Partial db-RDA single variables ####
### Space = location_abb ###
m1_substitution <- dbrda(
  beta_substitution ~ location_abb +
    Condition(pc3_soil + orientation + exposition + pc1_soil +
                pc1_construction_year),
  data = sites
  )
anova(m1_substitution, permutations = how(nperm = 9999)) # p: 1e-04 ***
RsquareAdj(m1_substitution) # R2adj: .275
### History = pc1_construction_year###
m1_substitution <- dbrda(
  beta_substitution ~ pc1_construction_year +
  Condition(pc3_soil + orientation + exposition + pc1_soil +
    location_abb),
  data = sites
  )
anova(m1_substitution, permutations = how(nperm = 9999)) # p: 1.5e-03**
RsquareAdj(m1_substitution) # R2adj: .049
### pc3_soil ###
m1_substitution <- dbrda(
  beta_substitution ~ pc3_soil +
  Condition(orientation + exposition + pc1_soil +
    location_abb +
    pc1_construction_year),
  data = sites
  )
anova(m1_substitution, permutations = how(nperm = 9999)) # p: 6.9e-03**
RsquareAdj(m1_substitution) # R2adj: .041
### Side ###
m1_substitution <- dbrda(
  beta_substitution ~ orientation +
  Condition(pc3_soil + exposition + pc1_soil +
    location_abb +
    pc1_construction_year),
  data = sites
  )
anova(m1_substitution, permutations = how(nperm = 9999)) # p: 7.0e-02.
RsquareAdj(m1_substitution) # R2adj: .019
### exposition ###
m1_substitution <- dbrda(
  beta_substitution ~ exposition +
  Condition(pc3_soil + orientation + pc1_soil +
    location_abb +
    pc1_construction_year),
  data = sites
  )
anova(m1_substitution, permutations = how(nperm = 9999)) # p: 3e-04***
RsquareAdj(m1_substitution) # R2adj: .069
### pc1_soil ###
m1_substitution <- dbrda(
  beta_substitution ~ pc1_soil +
  Condition(pc3_soil + orientation + exposition +
    location_abb +
    pc1_construction_year),
  data = sites
  )
anova(m1_substitution, permutations = how(nperm = 9999)) # p: 4.9e-01
RsquareAdj(m1_substitution) # R2adj: .000


### c Subsets ------------------------------------------------------------------

### * Full model ####
m1 <- dbrda(
  beta_subsets ~ pc1_soil + pc2_soil + pc3_soil + exposition + orientation +
  location_abb + river_km + river_distance +
  plot_age + pc1_construction_year + pc2_construction_year +
    pc3_construction_year,
  data = sites
  )
anova(m1, permutations = how(nperm = 999)) # P: .987
(r2adj <- RsquareAdj(m1)$adj.r.squared) # R2adj: -.580

### * Forward selection ####
# --> no forward selection because full model is not significant
