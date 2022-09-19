# Beta diversity on dike grasslands
# Variation partitioning of 2019 - All species ####
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
  na = c("", "na", "NA"), col_types =
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
  select(
    id, plot, block, location_construction_year, construction_year,
    longitude, latitude,
    exposition, orientation, pc1_soil, pc2_soil, pc3_soil,
    location_abb, river_km, river_distance,
    survey_year, plot_age, pc1_construction_year, pc2_construction_year,
    pc3_construction_year,
    accumulated_cover
  ) %>%
  filter(survey_year == 2019) %>%
  mutate(
    survey_year_factor = as_factor(survey_year),
    exposition_numeric = as.double(exposition),
    orientation_numeric = as.double(orientation),
    location_abb_numeric = as.double(location_abb)
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
  select(plot_age, pc1_construction_year, pc2_construction_year)



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
#--> MEM1_2019 ~ river_km and river_distance ~ pc3_construction_year have r > 0.7 (Dormann et al. 2013 Ecography)
#--> MEM1 and pc3_construction_year removed

### * Calculate: Baselga presence-absence ####
beta <- beta.div.comp(species, coef = "BS", quant = FALSE)
beta$Note
beta$part # total = 0.319, substitution = 0.275, subsets = 0.040
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


### b Substitution -------------------------------------------------------------

### * Check for linear trend in data ####
# m1 <- dbrda(beta_substitution ~ longitude + latitude, data = sites)
# anova(m1) #sig
# beta_substitution_detrended <- resid(lm(beta_substitution ~ longitude + latitude, data = sites))
#--> this trend is captured by river_km

### * full model ####
m1 <- dbrda(
  beta_substitution ~ pc1_soil + pc2_soil + pc3_soil + exposition +
    orientation +
  location_abb + river_km + river_distance +
  plot_age + pc1_construction_year + pc2_construction_year +
    pc3_construction_year,
  data = sites
  )
anova(m1, permutations = how(nperm = 9999)) # P: 1e-04
(r2adj <- RsquareAdj(m1)$adj.r.squared) # R2adj: .309

### * Forward selection soil ####
m1 <- dbrda(
  beta_substitution ~ pc1_soil + pc2_soil + pc3_soil + exposition + orientation,
  data = sites
  )
(r2adj <- RsquareAdj(m1)$adj.r.squared)
sel <- forward.sel(
  beta_substitution,
  sites_soil,
  adjR2thresh = r2adj,
  nperm = 9999
  )
sel$p_adj <- p.adjust(sel$pvalue, method = "holm", n = ncol(sites_soil))
sel # https://www.davidzeleny.net/anadat-r/doku.php/en:forward_sel_examples
sites_soil_selected <- sites %>%
  select(exposition_numeric, pc3_soil)

### * Forward selection space ####
m1 <- dbrda(
  beta_substitution ~ location_abb + river_km + river_distance,
  data = sites
  )
(r2adj <- RsquareAdj(m1)$adj.r.squared)
sel <- forward.sel(
  beta_substitution,
  sites_space,
  adjR2thresh = r2adj,
  nperm = 9999
  )
sel$p_adj <- p.adjust(sel$pvalue, method = "holm", n = ncol(sites_space))
sel # https://www.davidzeleny.net/anadat-r/doku.php/en:forward_sel_examples
sites_space_selected <- sites %>%
  select()

### * Forward selection history ####
m1 <- dbrda(
  beta_substitution ~ plot_age + pc1_construction_year + pc2_construction_year,
  data = sites
  )
(r2adj <- RsquareAdj(m1)$adj.r.squared)
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
  beta_substitution, sites_soil_selected, sites_history_selected
  )
tiff(
  here("outputs", "figures", "figure_4c_2019_800dpi_12x12cm.tiff"),
  res = 72, width = 12, height = 12, units = "cm", compression = "none"
  )
plot(
  m1_substitution_varpart,
  Xnames = c("Soil", "History"),
  cutoff = 0.01, digits = 2, bg = NA
  )
dev.off()

### * Partial db-RDA soil ####
m1_substitution <- dbrda(
  beta_substitution ~ exposition + pc3_soil +
  Condition(pc1_construction_year),
  data = sites
  )
anova(m1_substitution, permutations = how(nperm = 9999)) # p: 1.1e-02
RsquareAdj(m1_substitution) # R2adj: .062

### * Partial db-RDA single variables ####
### History / pc1_construction_year ###
m1_substitution <- dbrda(
  beta_substitution ~ pc1_construction_year +
  Condition(exposition + pc3_soil),
  data = sites
  )
anova(m1_substitution, permutations = how(nperm = 9999)) # p: 9.0e-03
RsquareAdj(m1_substitution) # R2adj: .034
### Exposition ###
m1_substitution <- dbrda(
  beta_substitution ~ exposition +
  Condition(pc3_soil +
    pc1_construction_year),
  data = sites
  )
anova(m1_substitution, permutations = how(nperm = 9999)) # p: 3.2e-03
RsquareAdj(m1_substitution) # R2adj: .046
### pc3_soil ###
m1_substitution <- dbrda(
  beta_substitution ~ pc3_soil +
  Condition(exposition +
    pc1_construction_year),
  data = sites
  )
anova(m1_substitution, permutations = how(nperm = 9999)) # p: 2.7e-01
RsquareAdj(m1_substitution) # R2adj: .006


### c Subsets ------------------------------------------------------------------

### * Full model ####
m1 <- dbrda(
  beta_subsets ~ pc1_soil + pc2_soil + pc3_soil + exposition + orientation +
  location_abb + river_km + river_distance +
  plot_age + pc1_construction_year + pc2_construction_year,
  data = sites
  )
anova(m1, permutations = how(nperm = 9999)) # P: 2.5e-01
(r2adj <- RsquareAdj(m1)$adj.r.squared) # R2adj: .185

### * Forward selection ####
# --> now forward selection because full model is not significant
