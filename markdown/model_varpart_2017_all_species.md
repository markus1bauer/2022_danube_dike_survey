Analysis of Bauer et al. (unpublished) Beta diversity on dike
grasslands: <br> Spatial variation 2017
================
<b>Markus Bauer</b> <br>
<b>2023-01-17</b>

- <a href="#preparation" id="toc-preparation">Preparation</a>
- <a href="#statistics" id="toc-statistics">Statistics</a>
  - <a href="#calculate-beta-diversity"
    id="toc-calculate-beta-diversity">Calculate beta diversity</a>
    - <a href="#check-collinearity" id="toc-check-collinearity">Check
      collinearity</a>
    - <a href="#calculate-baselga-presence-absence"
      id="toc-calculate-baselga-presence-absence">Calculate: Baselga
      presence-absence</a>
  - <a href="#db-rda-replacement-component"
    id="toc-db-rda-replacement-component">db-RDA: Replacement component</a>
    - <a href="#check-linear-trend-in-data"
      id="toc-check-linear-trend-in-data">Check linear trend in data</a>
    - <a href="#full-model" id="toc-full-model">Full model</a>
    - <a href="#forward-selection-soil"
      id="toc-forward-selection-soil">Forward selection: Soil</a>
    - <a href="#forward-selection-space"
      id="toc-forward-selection-space">Forward selection: Space</a>
    - <a href="#forward-selection-history"
      id="toc-forward-selection-history">Forward selection: History</a>
    - <a href="#variation-partitioning"
      id="toc-variation-partitioning">Variation partitioning</a>
    - <a href="#partial-db-rda" id="toc-partial-db-rda">Partial db-RDA</a>
  - <a href="#db-rda-nestedness-component"
    id="toc-db-rda-nestedness-component">db-RDA: Nestedness component</a>
    - <a href="#full-model-1" id="toc-full-model-1">Full model</a>
    - <a href="#forward-selection" id="toc-forward-selection">Forward
      selection</a>

<br/> <br/> <b>Markus Bauer</b>

Technichal University of Munich, TUM School of Life Sciences, Chair of
Restoration Ecology, Emil-Ramann-Straße 6, 85354 Freising, Germany

<markus1.bauer@tum.de>

ORCiD ID: [0000-0001-5372-4174](https://orcid.org/0000-0001-5372-4174)
<br> [Google
Scholar](https://scholar.google.de/citations?user=oHhmOkkAAAAJ&hl=de&oi=ao)
<br> GitHub: [markus1bauer](https://github.com/markus1bauer)

To compare different models, you only have to change the models in
section ‘Load models’

# Preparation

Borcard, Gillet & Legendre (2018) Numerical Ecology with R. 2nd edition.
Springer, Cham. [DOI:
10.1007/978-3-319-71404-2](https://doi.org/10.1007/978-3-319-71404-2)
Chapter 6.3

#### Packages

``` r
library(here)
library(tidyverse)
library(vegan)
library(adespatial)
```

#### Load data

``` r
sites <- read_csv(here("data", "processed", "data_processed_sites_spatial.csv"),
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
    location_abb, river_km, river_distance, biotope_distance, biotope_area,
    mem1_2017,
    survey_year, plot_age, pc1_construction_year, pc2_construction_year,
    pc3_construction_year,
    accumulated_cover
    ) %>%
  mutate(
    survey_year_factor = as_factor(survey_year),
    exposition_numeric = as.double(exposition),
    orientation_numeric = as.double(orientation),
    location_abb_numeric = as.double(location_abb),
    botanist_numeric = as.double(as_factor(botanist)),
    biotope_area = if_else(is.na(biotope_area), 0, biotope_area)
  )

species <- read_csv(here("data", "processed", "data_processed_species.csv"),
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
```

# Statistics

## Calculate beta diversity

### Check collinearity

Exclude r \> 0.7 <br> Dormann et al. 2013 Ecography [DOI:
10.1111/j.1600-0587.2012.07348.x](https://doi.org/10.1111/j.1600-0587.2012.07348.x)

``` r
sites %>%
  select(
    where(is.numeric), -ends_with("numeric"),
    -accumulated_cover, -construction_year, -survey_year,
    ) %>%
  GGally::ggpairs(lower = list(continuous = "smooth_loess"))
```

![](model_varpart_2017_all_species_files/figure-gfm/collinearity-1.png)<!-- -->

→ Remove longitude, latitude, biotope_area, mem1_2017,
pc3_construction_year

``` r
sites_soil <- sites %>%
  select(pc1_soil, pc2_soil, pc3_soil, exposition_numeric, orientation_numeric)
sites_space <- sites %>%
  select(location_abb_numeric, river_km, biotope_area)
sites_history <- sites %>%
  select(plot_age, pc1_construction_year, pc2_construction_year)
```

### Calculate: Baselga presence-absence

``` r
beta <- beta.div.comp(species, coef = "BS", quant = FALSE)
beta$Note
```

    ## [1] "Baselga family, Sorensen"

``` r
beta$part
```

    ##      BDtotal         Repl          Nes Repl/BDtotal  Nes/BDtotal 
    ##   0.31742082   0.27585014   0.04157069   0.86903605   0.13096395

``` r
beta_total <- beta$D %>% # = Soerensen dissimilarity
  as.matrix() %>%
  as.data.frame()
beta_substitution <- beta$repl %>% # = Replacement / Simpson dissimilarity
  as.matrix()
beta_subsets <- beta$rich %>% # = Nestedness
  as.matrix()
```

## db-RDA: Replacement component

### Check linear trend in data

``` r
m1 <- dbrda(beta_substitution ~ longitude + latitude, data = sites)
anova(m1)
```

    ## Permutation test for dbrda under reduced model
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Model: dbrda(formula = beta_substitution ~ longitude + latitude, data = sites)
    ##          Df SumOfSqs      F Pr(>F)   
    ## Model     2   0.5874 1.9218  0.008 **
    ## Residual 38   5.8073                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
beta_substitution_detrended <- resid(lm(beta_substitution ~ longitude + latitude, data = sites))
```

→ this trend is captured by river_km

### Full model

``` r
m1 <- dbrda(
  beta_substitution ~
    pc1_soil + pc2_soil + pc3_soil + exposition + orientation +
    location_abb + river_km + biotope_area +
    plot_age + pc1_construction_year + pc2_construction_year,
  data = sites
  )
anova(m1, permutations = how(nperm = 9999))
```

    ## Permutation test for dbrda under reduced model
    ## Permutation: free
    ## Number of permutations: 9999
    ## 
    ## Model: dbrda(formula = beta_substitution ~ pc1_soil + pc2_soil + pc3_soil + exposition + orientation + location_abb + river_km + biotope_area + plot_age + pc1_construction_year + pc2_construction_year, data = sites)
    ##          Df SumOfSqs      F Pr(>F)    
    ## Model    18   4.5851 3.0968  1e-04 ***
    ## Residual 22   1.8096                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
(r2adj <- RsquareAdj(m1)$adj.r.squared)
```

    ## [1] 0.485484

### Forward selection: Soil

``` r
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
```

    ## Testing variable 1
    ## Testing variable 2
    ## Testing variable 3
    ## Testing variable 4
    ## Testing variable 5
    ## Procedure stopped (adjR2thresh criteria) adjR2cum = 0.134818 with 5 variables (> 0.132081)

``` r
sel$p_adj <- p.adjust(sel$pvalue, method = "holm", n = ncol(sites_soil))
sel # https://www.davidzeleny.net/anadat-r/doku.php/en:forward_sel_examples
```

    ##             variables order         R2      R2Cum   AdjR2Cum        F pvalue
    ## 1            pc3_soil     3 0.07271581 0.07271581 0.04893929 3.058304 0.0017
    ## 2 orientation_numeric     5 0.05047211 0.12318792 0.07703991 2.187402 0.0140
    ## 3  exposition_numeric     4 0.04639787 0.16958578 0.10225490 2.067307 0.0221
    ## 4            pc1_soil     1 0.04526762 0.21485340 0.12761489 2.075579 0.0173
    ##    p_adj
    ## 1 0.0085
    ## 2 0.0560
    ## 3 0.0560
    ## 4 0.0560

``` r
sites_soil_selected <- sites %>%
  select(pc3_soil, orientation_numeric, exposition_numeric, pc1_soil)
```

### Forward selection: Space

``` r
m1 <- dbrda(
  beta_substitution ~ location_abb + river_km + biotope_area,
  data = sites
  )
r2adj <- RsquareAdj(m1)$adj.r.squared
sel <- forward.sel(
  beta_substitution,
  sites_space,
  adjR2thresh = r2adj,
  nperm = 9999
  )
```

    ## Testing variable 1
    ## Testing variable 2
    ## Procedure stopped (alpha criteria): pvalue for variable 2 is 0.076200 (> 0.050000)

``` r
sel$p_adj <- p.adjust(sel$pvalue, method = "holm", n = ncol(sites_space))
sel # https://www.davidzeleny.net/anadat-r/doku.php/en:forward_sel_examples
```

    ##              variables order         R2      R2Cum   AdjR2Cum        F pvalue
    ## 1 location_abb_numeric     1 0.07510911 0.07510911 0.05139396 3.167136 0.0012
    ##    p_adj
    ## 1 0.0036

``` r
sites_space_selected <- sites %>%
  select(location_abb_numeric)
```

### Forward selection: History

``` r
m1 <- dbrda(
  beta_substitution ~ plot_age + pc1_construction_year + pc2_construction_year,
  data = sites
  )
r2adj <- RsquareAdj(m1)$adj.r.squared
sel <- forward.sel(
  beta_substitution,
  sites_history,
  adjR2thresh = r2adj,
  nperm = 9999
  )
```

    ## Testing variable 1
    ## Testing variable 2
    ## Testing variable 3
    ## Procedure stopped (alpha criteria): pvalue for variable 3 is 0.564200 (> 0.050000)

``` r
sel$p_adj <- p.adjust(sel$pvalue, method = "holm", n = ncol(sites_history))
sel # https://www.davidzeleny.net/anadat-r/doku.php/en:forward_sel_examples
```

    ##               variables order         R2      R2Cum   AdjR2Cum        F pvalue
    ## 1 pc1_construction_year     2 0.05751076 0.05751076 0.03334437 2.379783 0.0120
    ## 2              plot_age     1 0.04308209 0.10059284 0.05325562 1.820220 0.0484
    ##    p_adj
    ## 1 0.0360
    ## 2 0.0968

``` r
sites_history_selected <- sites %>%
  select(pc1_construction_year)
```

### Variation partitioning

``` r
m1_substitution_varpart <- varpart(
  beta_substitution, sites_soil_selected,
  sites_space_selected, sites_history_selected
)
plot(
  m1_substitution_varpart,
  Xnames = c("Site", "Space", "History"),
  cutoff = 0.01, digits = 2, bg = NA
  )
```

![](model_varpart_2017_all_species_files/figure-gfm/varpart-1.png)<!-- -->

### Partial db-RDA

#### Soil

``` r
m1_substitution <- dbrda(
  beta_substitution ~ pc3_soil + orientation + exposition + pc1_soil +
    Condition(location_abb + pc1_construction_year),
data = sites
)
anova(m1_substitution, permutations = how(nperm = 9999))
```

    ## Permutation test for dbrda under reduced model
    ## Permutation: free
    ## Number of permutations: 9999
    ## 
    ## Model: dbrda(formula = beta_substitution ~ pc3_soil + orientation + exposition + pc1_soil + Condition(location_abb + pc1_construction_year), data = sites)
    ##          Df SumOfSqs      F Pr(>F)    
    ## Model     4   1.0806 2.8627  1e-04 ***
    ## Residual 27   2.5479                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
RsquareAdj(m1_substitution)
```

    ## $r.squared
    ## [1] 0.1689783
    ## 
    ## $adj.r.squared
    ## [1] 0.1418725

#### Space = location_abb

``` r
m1_substitution <- dbrda(
  beta_substitution ~ location_abb +
    Condition(pc3_soil + orientation + exposition + pc1_soil +
                pc1_construction_year),
  data = sites
  )
anova(m1_substitution, permutations = how(nperm = 9999))
```

    ## Permutation test for dbrda under reduced model
    ## Permutation: free
    ## Number of permutations: 9999
    ## 
    ## Model: dbrda(formula = beta_substitution ~ location_abb + Condition(pc3_soil + orientation + exposition + pc1_soil + pc1_construction_year), data = sites)
    ##          Df SumOfSqs      F Pr(>F)    
    ## Model     8   2.2978 3.0437  1e-04 ***
    ## Residual 27   2.5479                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
RsquareAdj(m1_substitution)
```

    ## $r.squared
    ## [1] 0.3593184
    ## 
    ## $adj.r.squared
    ## [1] 0.2757304

#### History = pc1_construction_year

``` r
m1_substitution <- dbrda(
  beta_substitution ~ pc1_construction_year +
  Condition(pc3_soil + orientation + exposition + pc1_soil +
    location_abb),
  data = sites
  )
anova(m1_substitution, permutations = how(nperm = 9999))
```

    ## Permutation test for dbrda under reduced model
    ## Permutation: free
    ## Number of permutations: 9999
    ## 
    ## Model: dbrda(formula = beta_substitution ~ pc1_construction_year + Condition(pc3_soil + orientation + exposition + pc1_soil + location_abb), data = sites)
    ##          Df SumOfSqs      F Pr(>F)   
    ## Model     1  0.31815 3.3715 0.0024 **
    ## Residual 27  2.54788                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
RsquareAdj(m1_substitution)
```

    ## $r.squared
    ## [1] 0.04975221
    ## 
    ## $adj.r.squared
    ## [1] 0.04999346

#### PC3_soil

``` r
m1_substitution <- dbrda(
  beta_substitution ~ pc3_soil +
  Condition(orientation + exposition + pc1_soil +
    location_abb +
    pc1_construction_year),
  data = sites
  )
anova(m1_substitution, permutations = how(nperm = 9999))
```

    ## Permutation test for dbrda under reduced model
    ## Permutation: free
    ## Number of permutations: 9999
    ## 
    ## Model: dbrda(formula = beta_substitution ~ pc3_soil + Condition(orientation + exposition + pc1_soil + location_abb + pc1_construction_year), data = sites)
    ##          Df SumOfSqs     F Pr(>F)   
    ## Model     1  0.27885 2.955 0.0048 **
    ## Residual 27  2.54788                
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
RsquareAdj(m1_substitution)
```

    ## $r.squared
    ## [1] 0.04360581
    ## 
    ## $adj.r.squared
    ## [1] 0.04121288

#### Side

``` r
m1_substitution <- dbrda(
  beta_substitution ~ orientation +
  Condition(pc3_soil + exposition + pc1_soil +
    location_abb +
    pc1_construction_year),
  data = sites
  )
anova(m1_substitution, permutations = how(nperm = 9999))
```

    ## Permutation test for dbrda under reduced model
    ## Permutation: free
    ## Number of permutations: 9999
    ## 
    ## Model: dbrda(formula = beta_substitution ~ orientation + Condition(pc3_soil + exposition + pc1_soil + location_abb + pc1_construction_year), data = sites)
    ##          Df SumOfSqs     F Pr(>F)  
    ## Model     1  0.18288 1.938 0.0646 .
    ## Residual 27  2.54788               
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
RsquareAdj(m1_substitution)
```

    ## $r.squared
    ## [1] 0.02859859
    ## 
    ## $adj.r.squared
    ## [1] 0.01977401

#### Exposition

``` r
m1_substitution <- dbrda(
  beta_substitution ~ exposition +
  Condition(pc3_soil + orientation + pc1_soil +
    location_abb +
    pc1_construction_year),
  data = sites
  )
anova(m1_substitution, permutations = how(nperm = 9999))
```

    ## Permutation test for dbrda under reduced model
    ## Permutation: free
    ## Number of permutations: 9999
    ## 
    ## Model: dbrda(formula = beta_substitution ~ exposition + Condition(pc3_soil + orientation + pc1_soil + location_abb + pc1_construction_year), data = sites)
    ##          Df SumOfSqs     F Pr(>F)    
    ## Model     1  0.40643 4.307  2e-04 ***
    ## Residual 27  2.54788                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
RsquareAdj(m1_substitution)
```

    ## $r.squared
    ## [1] 0.06355758
    ## 
    ## $adj.r.squared
    ## [1] 0.06971542

#### PC1_soil

``` r
m1_substitution <- dbrda(
  beta_substitution ~ pc1_soil +
  Condition(pc3_soil + orientation + exposition +
    location_abb +
    pc1_construction_year),
  data = sites
  )
anova(m1_substitution, permutations = how(nperm = 9999))
```

    ## Permutation test for dbrda under reduced model
    ## Permutation: free
    ## Number of permutations: 9999
    ## 
    ## Model: dbrda(formula = beta_substitution ~ pc1_soil + Condition(pc3_soil + orientation + exposition + location_abb + pc1_construction_year), data = sites)
    ##          Df SumOfSqs      F Pr(>F)
    ## Model     1   0.0916 0.9707 0.4949
    ## Residual 27   2.5479

``` r
RsquareAdj(m1_substitution)
```

    ## $r.squared
    ## [1] 0.01432502
    ## 
    ## $adj.r.squared
    ## [1] -0.0006168162

## db-RDA: Nestedness component

### Full model

``` r
m1 <- dbrda(
  beta_subsets ~ pc1_soil + pc2_soil + pc3_soil + exposition + orientation +
  location_abb + river_km + river_distance + biotope_distance +
  plot_age + pc1_construction_year + pc2_construction_year,
  data = sites
  )
anova(m1, permutations = how(nperm = 999))
```

    ## Permutation test for dbrda under reduced model
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Model: dbrda(formula = beta_subsets ~ pc1_soil + pc2_soil + pc3_soil + exposition + orientation + location_abb + river_km + river_distance + biotope_distance + plot_age + pc1_construction_year + pc2_construction_year, data = sites)
    ##          Df SumOfSqs      F Pr(>F)
    ## Model    19 0.041094 0.2101  0.982
    ## Residual 21 0.216182

``` r
(r2adj <- RsquareAdj(m1)$adj.r.squared)
```

    ## [1] -0.6005206

### Forward selection

→ no forward selection because full model is not significant
