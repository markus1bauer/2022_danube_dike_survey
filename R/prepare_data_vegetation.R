# Prepare vegetation data ####
# Markus Bauer

## Content
# 1 Create simple variables
# 2 Coverages
# 3 Species richness
# 4 Biotope types
# 5 CWM off Ellenberg
# 6 Load functional plant traits
# 7 Calculate CWM and FDis
# 8 Beta diversity
# 9 Environmental variables


### Packages ###
library(here)
library(tidyverse)
library(naniar) #are_na()
library(lubridate) #modify dates
library(vegan)
library(FD) #dbFD()
library(adespatial)
#remotes::install_github("larsito/tempo")
library(tempo)


### Start ###3-trs5mF2122
#installr::updateR(browse_news = F, install_R = T, copy_packages = T, copy_Rprofile.site = T, keep_old_packages = T, update_packages = T, start_new_R = F, quit_R = T, print_R_versions = T, GUI = F)
#sessionInfo()
rm(list = ls())
setwd(here("data/raw"))



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Load data ##############################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


### 1 Sites #####################################################################################

sites <- read_csv("data_raw_sites.csv", col_names = T, na = c("", "NA", "na"), col_types = 
                    cols(
                      .default = "?",
                      id = "f",
                      block = "f",
                      location = "f",
                      side = "f",
                      exposition = "f",
                      ageCategory = "f",
                      surveyDate_2017 = col_date(format = "%Y-%m-%d"),
                      surveyDate_2018 = col_date(format = "%Y-%m-%d"),
                      surveyDate_2019 = col_date(format = "%Y-%m-%d"),
                      HCl = "f",
                      humusLevel = "f"
                    )) %>%
  mutate(across(starts_with("vegetationCov") | 
                  starts_with("vegetationHeight_") | 
                  starts_with("opensoilCov_") | 
                  starts_with("mossCov_") | 
                  starts_with("litterCov_") | 
                  starts_with("botanist_"),
                ~ as.character(.x))) %>%
  pivot_longer(starts_with("vegetationCov") | 
                 starts_with("vegetationHeight_") | 
                 starts_with("opensoilCov_") | 
                 starts_with("mossCov_") | 
                 starts_with("litterCov_") | 
                 starts_with("botanist_"), 
               names_to = c("x", "surveyYear"),
               names_sep = "_",
               values_to ="n") %>%
  pivot_wider(names_from = x, values_from = n) %>%
  mutate(across(c(surveyYear, vegetationCov, vegetationHeight, opensoilCov, mossCov, litterCov),
                ~ as.numeric(.x))) %>%
  mutate(surveyYearF = factor(surveyYear),
         surveyYearFminus = factor(surveyYear - 1),
         constructionYearF = factor(constructionYear),
         constructionYearFplus = factor(constructionYear + 1)) %>%
  mutate(id = str_c(id, surveyYear, sep = "_"), 
         .keep = "all",
         id = paste0("X", id),
         plot = str_sub(id, start = 2, end = 3),
         position = str_sub(id, start = 5, end = 5),
         locationAbb = str_sub(location, 1, 3),
         locationAbb = str_to_upper(locationAbb),
         locationAbb = factor(locationAbb, levels = unique(locationAbb[order(constructionYear)])),
         locationYear = str_c(locationAbb, constructionYear, sep = "-")) %>%
  select(-starts_with("surveyDate_"), -starts_with("topsoilDepth_"), -cnLevel, -ends_with("Class"), -starts_with("sceleton"), -position)


### 2 Species #####################################################################################

species <- data.table::fread("20210915_data_raw_species.csv", 
                             sep = ",",
                             dec = ".",
                             skip = 0,
                             header = T,
                             na.strings = c("", "NA", "na"),
                             colClasses = list(
                               character = "name"
                             )) %>%
  ### Check that each species occurs at least one time ###
  group_by(name) %>%
  arrange(name) %>%
  mutate(total = sum(c_across(starts_with("X")), na.rm = T),
         presence = if_else(total > 0, 1, 0),
         name = factor(name)) %>%
  filter(presence == 1) %>%
  ungroup() %>%
  select(name, sort(tidyselect::peek_vars()), -total, -presence) %>%
  select(name, all_of(sites$id)) %>%
  na_if(0)

### Create list with species names and their frequency ###
specieslist <- species %>%
  mutate_if(is.numeric, ~1 * (. != 0)) %>%
  mutate(sum = rowSums(across(where(is.numeric)), na.rm = T), .keep = "unused") %>%
  group_by(name) %>%
  summarise(sum = sum(sum))
#write_csv(specieslist, "specieslist.csv")


### 3 Traits #####################################################################################

traits <- read_csv("data_raw_traits.csv", col_names = T, na = c("", "NA", "na"), col_types = 
                      cols(
                        .default = "f",
                        name = "c",
                        sociology = "d",
                        l = "d",
                        t = "d",
                        k = "d",
                        f = "d",
                        r = "d",
                        n = "d"
                      )) %>%
  separate(name, c("genus", "species", "ssp", "subspecies"), "_", remove = F, extra = "drop", fill = "right") %>%
  mutate(genus = str_sub(genus, 1, 4),
         species = str_sub(species, 1, 4),
         subspecies = str_sub(subspecies, 1, 4),
         name = factor(name)) %>%
  unite(abb, genus, species, subspecies, sep = "") %>%
  mutate(abb = str_replace(abb, "NA", ""),
         abb = as_factor(abb)) %>%
  arrange(name)
### Check congruency of traits and species table ###
traits[duplicated(traits$abb),]
#traits$name[which(!(traits$name %in% species$name))]
species$name[which(!(species$name %in% traits$name))]
traits <- semi_join(traits, species, by = "name")


### 4 Check data frames #####################################################################################

### Check typos ###
sites %>%
  filter(!str_detect(id, "_seeded$")) %>%
  janitor::tabyl(vegetationCov)
#sites %>% filter(vegetationCov == 17)
species %>%
  select(-name) %>%
  unlist() %>%
  janitor::tabyl()
species %>% # Check special typos
  pivot_longer(-name, names_to = "id", values_to = "value") %>%
  filter(value == 8)

### Compare vegetationCov and accumulatedCov ###
species %>%
  summarise(across(where(is.double),~sum(.x, na.rm = T))) %>%
  pivot_longer(cols = everything(), names_to = "id", values_to = "value") %>%
  mutate(id = factor(id)) %>%
  full_join(sites, by = "id") %>% 
  mutate(diff = (value - vegetationCov)) %>%
  select(id, surveyYear, value, vegetationCov, diff) %>%
  filter(diff > 50 | diff < -30) %>%
  arrange(diff) %>%
  print(n = 100)

### Check plots over time ###
species %>%
  select(name, starts_with("X66")) %>%
  filter(if_any(starts_with("X"), ~ . > 0)) %>%
  print(n = 70)

### Check missing data ###
miss_var_summary(sites, order = T)
vis_miss(sites, cluster = F)
vis_miss(traits, cluster = F, sort_miss = T)

rm(list = setdiff(ls() ,c("species", "traits", "sites")))



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Create variables ##############################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


### 1 Create simple variables #####################################################################################

traits <- traits %>%
  mutate(leanIndicator = if_else(
    !(is.na(table30)) | !(is.na(table33)) | !(is.na(table34)), "yes", "no"
    )) %>%
  mutate(target = if_else(
    targetHerb == "yes" | targetGrass == "yes", "yes", "no"
    )) %>%
  mutate(ruderal = if_else(
    sociology >= 3300 & sociology < 3700, "yes", "no"
    )) %>%
  mutate(targetEllenberg = if_else(
    sociology >= 5300 & sociology < 5400, "dry_grassland", if_else(
      sociology >= 5400 & sociology < 6000, "hay_meadow", if_else(
        sociology >= 5100 & sociology < 5200, "nardus_grassland", if_else(
          sociology >= 5200 & sociology < 5300, "sand_grasland", "no"
          )))))

sites <- sites %>%
  mutate(conf.low = c(1:length(id)),
         conf.high = c(1:length(id)),
         fmMass = round(fmDepth * fmDbd * 10, 3),
         NtotalConc = round(fmMass * NtotalPerc / 100, 3),
         plotAge = surveyYear - constructionYear) %>%
  select(-fmDepth, -fmMass)


### 2 Coverages #####################################################################################

cover <- left_join(species, traits, by = "name") %>%
  select(name, family, target, targetHerb, targetArrhenatherion, leanIndicator, nitrogenIndicator, ruderalIndicator, table33, starts_with("X")) %>%
  pivot_longer(names_to = "id", values_to = "n", cols = starts_with("X")) %>%
  group_by(id)

### * graminoid, herb, and total coverage) ####
cover_total_and_graminoid <- cover %>%
  group_by(id, family) %>%
  summarise(total = sum(n, na.rm = T)) %>%
  mutate(type = if_else(family == "Poaceae" | family == "Cyperaceae" | family == "Juncaceae", "graminoidCov", "herbCov")) %>%
  group_by(id, type) %>%
  summarise(total = sum(total, na.rm = T)) %>%
  spread(type, total) %>%
  mutate(accumulatedCov = graminoidCov + herbCov,
         accumulatedCov = round(accumulatedCov, 1)) %>%
  ungroup()

### * Target species' coverage ####
cover_target <- cover %>%
  filter(target == "yes") %>%
  summarise(targetCov = sum(n, na.rm = T)) %>%
  mutate(targetCov = round(targetCov, 1)) 
  ungroup() %>%

### * Target herb species' coverage ####
cover_targetHerb <- cover %>%
  filter(targetHerb == "yes") %>%
  summarise(targetHerbCov = sum(n, na.rm = T)) %>%
  mutate(targetHerbCov = round(targetHerbCov, 1))
  ungroup() %>%

### * Arrhenatherum species' cover ratio ####
cover_targetArrhenatherion <- cover %>%
  filter(targetArrhenatherion == "yes") %>%
  summarise(arrhCov = sum(n, na.rm = T)) %>%
  mutate(arrhCov = round(arrhCov, 1)) 
  ungroup() %>%

### * Lean indicator's coverage ####
cover_leanIndicator <- cover %>%
  filter(leanIndicator == "yes") %>%
  summarise(leanCov = sum(n, na.rm = T)) %>%
  mutate(leanCov = round(leanCov, 1))
  ungroup() %>%

### * Nitrogen indicator's coverage ####
cover_nitrogenIndicator <- cover %>%
  filter(nitrogenIndicator == "yes") %>%
  summarise(nitrogenCov = sum(n, na.rm = T)) %>%
  mutate(nitrogenCov = round(nitrogenCov, 1))
  ungroup() %>%

### * Ruderal indicator's coverage ####
cover_ruderalIndicator <- cover %>%
  filter(ruderalIndicator == "yes") %>%
  summarise(ruderalCov = sum(n, na.rm = T)) %>%
  mutate(ruderalCov = round(ruderalCov, 1))
  ungroup() %>%
  
### * Table 33 species' coverage ####
cover_table33 <- cover %>%
  mutate(table33 = if_else(table33 == "4" | table33 == "3" | table33 == "2", "table33Cov", "other")) %>%
  filter(table33 == "table33Cov") %>%
  summarise(table33Cov = sum(n, na.rm = T)) %>%
  mutate(table33Cov = round(table33Cov, 1))
  ungroup() %>%

### * implement in sites data set ####
sites <- sites %>%
  right_join(cover_total_and_graminoid, by = "id") %>%
  right_join(cover_target, by = "id") %>%
  right_join(cover_targetHerb, by = "id") %>%
  right_join(cover_targetArrhenatherion, by = "id") %>%
  right_join(cover_leanIndicator, by = "id") %>%
  right_join(cover_nitrogenIndicator, by = "id") %>%
  right_join(cover_ruderalIndicator, by = "id") %>%
  right_join(cover_table33, by = "id") %>%
  ### Calcute the ratio of target species richness of total species richness
  mutate(targetCovratio = targetCov / accumulatedCov,
         graminoidCovratio = graminoidCov / accumulatedCov,
         targetCovratio = round(targetCovratio, 3),
         graminoidCovratio = round(graminoidCovratio, 3))

rm(list = setdiff(ls(), c("sites", "species", "traits")))

### Pruefe auf Unregelmäßigkeiten bei Gesamtschätzungen ###
sites %>%
  mutate(difference = accumulatedCov - vegetationCov) %>%
  select(id, difference) %>%
  filter(difference > 20 | difference < 5)


### 3 Alpha diversity #####################################################################################

### a Species richness ----------------------------------------------------------------------------------------------------
speciesRichness <- left_join(species, traits, by = "name") %>%
  select(name, rlg, rlb, target, targetHerb, targetArrhenatherion, ffh6510, ffh6210, nitrogenIndicator, leanIndicator, table33, table34, starts_with("X")) %>%
  pivot_longer(names_to = "id", values_to = "n", cols = starts_with("X")) %>%
  mutate(n = if_else(n > 0, 1, 0)) %>%
  group_by(id)

### * total species richness ####
speciesRichness_all <- speciesRichness %>%
  summarise(speciesRichness = sum(n, na.rm = T)) %>%
  ungroup()

### * red list Germany (species richness) ####
speciesRichness_rlg <- speciesRichness %>%
  filter(rlg == "1" | rlg == "2" | rlg == "3" | rlg == "V") %>%
  summarise(rlgRichness = sum(n, na.rm = T)) %>%
  ungroup()

### * red list Bavaria (species richness) ####
speciesRichness_rlb <- speciesRichness %>%
  filter(rlb == "1" | rlb == "2" | rlb == "3" | rlb == "V") %>%
  summarise(rlbRichness = sum(n, na.rm = T)) %>%
  ungroup()

### * target species (species richness) ####
speciesRichness_target <- speciesRichness %>%
  filter(target == "yes") %>%
  summarise(targetRichness = sum(n, na.rm = T)) %>%
  ungroup()

### * target herb species (species richness) ####
speciesRichness_targetHerb <- speciesRichness %>%
  filter(targetHerb == "yes") %>%
  summarise(targetHerbRichness = sum(n, na.rm = T)) %>%
  ungroup()

### * Arrhenatherion species (species richness) ####
speciesRichness_arrh <- speciesRichness %>%
  filter(targetArrhenatherion == "yes") %>%
  summarise(arrhRichness = sum(n, na.rm = T)) %>%
  ungroup()

### * ffh6510 species (species richness) ####
speciesRichness_ffh6510 <- speciesRichness %>%
  filter(ffh6510 == "yes") %>%
  summarise(ffh6510Richness = sum(n, na.rm = T)) %>%
  ungroup()

### * ffh6210 species (species richness) ####
speciesRichness_ffh6210 <- speciesRichness %>%
  filter(ffh6210 == "yes") %>%
  summarise(ffh6210Richness = sum(n, na.rm = T)) %>%
  ungroup()

### * leanIndicator species (species richness) ####
speciesRichness_leanIndicator <- speciesRichness %>%
  filter(leanIndicator == "yes") %>%
  summarise(leanIndicatorRichness = sum(n, na.rm = T)) %>%
  ungroup()

### * table33 (species richness) ####
speciesRichness_table33_2 <- speciesRichness %>%
  filter(table33 == "2") %>%
  summarise(table33_2Richness = sum(n, na.rm = T)) %>%
  ungroup()
speciesRichness_table33_3 <- speciesRichness %>%
  filter(table33 == "3") %>%
  summarise(table33_3Richness = sum(n, na.rm = T)) %>%
  ungroup()
speciesRichness_table33_4 <- speciesRichness %>%
  filter(table33 == "4") %>%
  summarise(table33_4Richness = sum(n, na.rm = T)) %>%
  ungroup()

### * table34 (species richness) ####
speciesRichness_table34_2 <- speciesRichness %>%
  filter(table34 == "2") %>%
  summarise(table34_2Richness = sum(n, na.rm = T)) %>%
  ungroup()
speciesRichness_table34_3 <- speciesRichness %>%
  filter(table34 == "3") %>%
  summarise(table34_3Richness = sum(n, na.rm = T)) %>%
  ungroup()

### * implement in sites data set (species richness) ####
sites <- sites %>%
  right_join(speciesRichness_all, by = "id") %>%
  right_join(speciesRichness_rlg, by = "id") %>%
  right_join(speciesRichness_rlb, by = "id") %>%
  right_join(speciesRichness_target, by = "id") %>%
  right_join(speciesRichness_targetHerb, by = "id") %>%
  right_join(speciesRichness_arrh, by = "id") %>%
  right_join(speciesRichness_ffh6510, by = "id") %>%
  right_join(speciesRichness_ffh6210, by = "id") %>%
  right_join(speciesRichness_leanIndicator, by = "id") %>%
  right_join(speciesRichness_table33_2, by = "id") %>%
  right_join(speciesRichness_table33_3, by = "id") %>%
  right_join(speciesRichness_table33_4, by = "id") %>%
  right_join(speciesRichness_table34_2, by = "id") %>%
  right_join(speciesRichness_table34_3, by = "id") %>%
### Calcute the ratio of target species richness of total species richness
  mutate(targetRichratio = targetRichness / speciesRichness,
         targetRichratio = round(targetRichratio, 3))

rm(list = setdiff(ls(), c("sites", "species", "traits")))


### b Species eveness and shannon ----------------------------------------------------------------------
sites <- species %>%
  pivot_longer(-name, names_to = "id", values_to = "value") %>%
  pivot_wider(names_from = "name", values_from = "value") %>%
  column_to_rownames("id") %>%
  diversity(index = "shannon") %>%
  as_tibble(rownames = NA) %>%
  rownames_to_column(var = "id") %>%
  mutate(id = factor(id)) %>%
  rename(shannon = value) %>%
  right_join(sites, by = "id") %>%
  mutate(eveness = shannon / log(speciesRichness))
  


### 4 Biotope types #####################################################################################

### a Calculate types -------------------------------------------------------------------------------------------
biotopetypes <- sites %>%
  select(id, table33_2Richness, table33_3Richness, table33_4Richness, table33Cov, table34_2Richness, table34_3Richness, targetRichness, targetHerbRichness, arrhRichness, targetCov, leanCov, arrhCov, targetHerbCov, nitrogenCov) %>%
  mutate(table33Rich_proof = if_else(
    table33_2Richness >= 2 | table33_2Richness + table33_3Richness >= 3 | table33_2Richness + table33_3Richness + table33_4Richness >= 4, "yes", "no"
    )) %>%
  mutate(table33Cov_proof = if_else(
    table33Cov >= 25, "yes", "no"
    )) %>%
  mutate(table33_proof = if_else(
    table33Rich_proof == "yes" & table33Cov_proof == "yes", "yes", "no"
    )) %>%
  mutate(table34Rich_proof = if_else(
    table34_2Richness >= 2 | table34_2Richness + table34_3Richness >= 3, "yes", "no"
    )) %>%
  mutate(G312_GT6210_type = if_else(
    table33_proof == "yes" & table34Rich_proof == "yes", "yes", "no"
    )) %>%
  mutate(GE_proof = if_else(
    targetHerbCov >= 12.5 & targetRichness >= 20 & nitrogenCov < 25, "yes", "no"
    )) %>%
  mutate(G214_GE6510_type = if_else(
    GE_proof == "yes" & arrhRichness >= 1 & arrhCov > 0.5 & leanCov >= 25 & targetHerbCov >= 12.5, "yes", "no"
    )) %>%
  mutate(G212_LR6510_type = if_else(
    GE_proof == "yes" & arrhRichness >= 1 & arrhCov > 0.5 & leanCov < 25, "yes", "no"
    )) %>%
  mutate(G214_GE00BK_type = if_else(
    GE_proof == "yes" & arrhRichness == 0 & leanCov >= 25 & targetHerbCov >= 12.5, "yes", if_else(
      GE_proof == "yes" & arrhCov >= 0.5 & leanCov >= 25 & targetHerbCov >= 12.5, "yes", "no"
    )
      )) %>%
  mutate(G213_GE00BK_type = if_else(
    GE_proof == "yes" & arrhRichness == 0 & leanCov >= 25 & targetHerbCov < 12.5, "yes", if_else(
      GE_proof == "yes" & arrhCov >= 0.5 & leanCov >= 25 & targetHerbCov < 12.5, "yes", "no"
    )
  )) %>%
  mutate(G213_type = if_else(
    leanCov >= 25, "yes", "no"
    )) %>%
  mutate(G212_type = if_else(
    leanCov >= 1 & leanCov < 25 & targetHerbCov >= 12.5 & targetHerbRichness >= 10, "yes", "no"
    )) %>%
  mutate(G211_type = if_else(
    leanCov >= 1 & leanCov < 25 & targetHerbCov >= 1 & targetHerbRichness >= 5, "yes", "no"
    )) %>%
  mutate(biotopeType = if_else(
    G312_GT6210_type == "yes", "G312-GT6210", if_else(
      G214_GE6510_type == "yes", "G214-GE6510", if_else(
        G212_LR6510_type == "yes", "G212-LR6510", if_else(
          G214_GE00BK_type == "yes", "G214-GE00BK", if_else(
            G213_GE00BK_type == "yes", "G213-GE00BK", if_else(
              G213_type == "yes", "G213", if_else(
                G212_type == "yes", "G212", if_else(
                  G211_type == "yes", "G211", "other"
                  ))))))))) %>%
  mutate(biotopeType = as_factor(biotopeType)) %>%
  select(id, biotopeType, -ends_with("proof"), -starts_with("table")) %>%
  mutate(ffh6510 = str_match(biotopeType, "6510")) %>%
  mutate(ffh6210 = str_match(biotopeType, "6210")) %>%
  mutate(baykompv = as_factor(str_sub(biotopeType, start = 1, end = 4))) %>%
  unite(ffh, ffh6510, ffh6210, sep = "") %>%
  mutate(ffh = str_replace(ffh, "NA", "")) %>%
  mutate(ffh = as_factor(str_replace(ffh, "NA", "non-FFH"))) %>%
  mutate(biotopeType = as_factor(biotopeType)) %>%
  mutate(biotopePoints = if_else(
    biotopeType == "G312-GT6210", 13, if_else(
      biotopeType == "G214-GE6510", 12, if_else(
        biotopeType == "G212-LR6510", 9, if_else(
          biotopeType == "G214-GE00BK", 12, if_else(
            biotopeType == "G213-GE00BK", 9, if_else(
              biotopeType == "G213", 8, if_else(
                biotopeType == "G212", 8, if_else(
                  biotopeType == "G211", 6, 0
                ))))))))) %>%
  mutate(min8 = as_factor(if_else(biotopePoints >= 8, "yes", "no"))) %>%
  mutate(min9 = as_factor(if_else(biotopePoints >= 9, "yes", "no")))
sites <- left_join(sites, biotopetypes, by = "id") %>%  
#test <- left_join(sites, biotopetypes, by = "id") %>%
 # mutate_if(is.factor, as.character) %>%
  #select(id, biotopeType.x, biotopeType.y, arrhCov, arrhRichness, leanCov) %>%
  #mutate(similar = if_else(biotopeType.x == biotopeType.y, "yes", "no")) %>%
  #filter(similar == "no") %>%
  #select(-similar)
  #rm(test)
  select(-targetHerbCov, -arrhCov, -leanCov, -nitrogenCov, -table33Cov, -targetHerbRichness, -arrhRichness, -leanIndicatorRichness, -table33_2Richness, -table33_3Richness, table33_4Richness, table34_2Richness, table34_3Richness)
traits <- traits %>%
  select(-targetArrhenatherion, -table30, -table33, -table34, -nitrogenIndicator, -nitrogenIndicator2, leanIndicator)
rm(biotopetypes)  

### b Calculate constance -------------------------------------------------------------------------------------------
data <- sites %>%
  select(id, plot, surveyYear, ffh) %>%
  group_by(plot) %>%
  mutate(count = n()) %>%
  filter(count == max(count)) %>%
  pivot_wider(id_cols = -id, names_from = "surveyYear", values_from = "ffh") %>%  #group_by(plot) %>%
  rename(x17 = "2017", x18 = "2018", x19 = "2019") %>%
  mutate(changeType = ifelse((x17 == "non-FFH" & x18 != "non-FFH" & x19 != "non-FFH"), "better", 
                       ifelse((x17 != "non-FFH" & x18 == "non-FFH" & x19 != "non-FFH") | (x17 == "non-FFH" & x18 != "non-FFH" & x19 == "non-FFH") | (x17 == "non-FFH" & x18 == "non-FFH" & x19 != "non-FFH") | (x17 != "non-FFH" & x18 != "non-FFH" & x19 == "non-FFH"), "change", 
                              ifelse((x17 == "non-FFH" & x18 == "non-FFH" & x19 == "non-FFH"), "non-FFH", 
                                     ifelse(x17 == "6510" & x18 == "6510" & x19 == "6510", "FFH6510", 
                                            ifelse(x17 == "6210" & x18 == "6210" & x19 == "6210", "FFH6210",
                                                   ifelse(x17 != "non-FFH" & x18 != "non-FFH" & x19 != "non-FFH", "any-FFH", "worse"))))))) %>%
  select(plot, changeType)
sites <- data %>%
  right_join(sites, by = "plot")


### 5 Beta diversity #####################################################################################

### a NMDS -------------------------------------------------------------------------------------------
### Prepare data ###
data <- species %>%
  pivot_longer(-name, "id", "value") %>%
  pivot_wider(id, name) %>%
  arrange(id) %>%
  column_to_rownames("id")
### Calculate NMDS ###
set.seed(1)
(nmds <- metaMDS(data, dist = "bray", binary = F,
                 try = 50, previous.best = T, na.rm = T))
### Add to sites ###
sites <- nmds %>%
  scores() %>%
  as.data.frame() %>%
  rownames_to_column(var = "id") %>%
  as_tibble() %>%
  select(id, NMDS1, NMDS2) %>%
  right_join(sites, by = "id")

### b PERMDISP -------------------------------------------------------------------------------------------
(permdisp <- betadisper(d = vegdist(nmds), 
                        group = sites$plot))
sites <- mutate(sites, 
                permdisp = as.numeric(permdisp$distances))

### c TBI -------------------------------------------------------------------------------------------
### remove sites which are not surveyed every year ###
tbisites <- sites %>%
  select(id, plot, block, locationAbb, surveyYear, constructionYear) %>%
  add_count(plot) %>%
  filter(n == max(n)) %>%
  select(-n) 
which(table(tbisites$plot) == 0) # these plots are not complete
tbispecies <- species %>%
  pivot_longer(-name, "id", "value") %>%
  pivot_wider(id, name) %>%
  arrange(id) %>%
  semi_join(tbisites, by = "id") %>%
  column_to_rownames("id") %>%
  select(which(!colSums(., na.rm = T) %in% 0)) %>% # excludes species with zero abundances
  rownames_to_column(var = "id") %>%
  mutate(year = factor(str_match(id, "\\d\\d\\d\\d"))) %>%
  mutate(plot = factor(str_match(id, "\\d\\d"))) %>%
  select(-id)
### Separate each year in several tibbles ###
species17 <- tbispecies %>%
  filter(year == 2017) %>%
  column_to_rownames("plot") %>%
  select(-year)
species18 <- tbispecies %>%
  filter(year == 2018) %>%
  column_to_rownames("plot") %>%
  select(-year)
species19 <- tbispecies %>%
  filter(year == 2019) %>%
  column_to_rownames("plot") %>%
  select(-year)
rm(tbisites, tbispecies)

### c Synchrony -------------------------------------------------------------------------------------------
data <- species %>%
  pivot_longer(-name, "id", "value") %>%
  pivot_wider(id, name) %>%
  arrange(id) %>%
  mutate(plot = factor(str_sub(id, 1, 3))) %>%
  column_to_rownames(var = "id") %>%
  add_count(plot) %>%
  filter(n > 2) %>%
  select(-n)
data <- data %>%
  split(data$plot, drop = T) %>%
  map(~ (.x %>% select(-plot)))
data <- lapply(data, calc_sync)
data <- do.call("rbind", sync_indices) %>%
  as.data.frame() %>%
  rownames_to_column("plot") %>%
  as_tibble() %>%
  mutate(plot = str_extract(plot, "\\d\\d")) %>%
  select(plot, log_varrat, log_varrat_t3, syn_total, syn_trend, syn_detrend)
sites <- full_join(sites, data, by = "plot")


### 6 Environmental variables #####################################################################################

### a Soil PCA  -------------------------------------------------------------------------------------------
### Prepare data ###
data <- sites %>%
  select(id, plot, calciumcarbonatPerc, NtotalPerc, cnRatio, pH, sandPerc, siltPerc, clayPerc, phosphorous, potassium, magnesium, topsoilDepth, NtotalConc) %>%
  group_by(plot) %>%
  summarise(across(where(is.numeric), ~median(.x, na.rm = T))) %>%
  select(-plot)
### Calculate PCA ###
pca <- rda(X = decostand(data, method = "standardize"), scale = T)
biplot(pca, display = "species")
screeplot(pca, bstick = TRUE, type = "l", main = NULL)
### Make data frames ###
eigenvals <- pca %>%
  eigenvals() %>%
  summary() %>%
  as_tibble() %>%
  select(PC1:PC3) %>%
  bind_cols(c("Eigenvalues", "Proportion Explained", "Cumulative Proportion")) %>%
  rename("variables" = "...4") %>%
  mutate(across(where(is.numeric), as.double))
values <- pca %>%
  summary()
pca <- values$species[ ,1:3] %>%
  as_tibble() %>%
  bind_cols(c("calciumcarbonatPerc", "NtotalPerc", "cnRatio", "pH", "sandPerc", "siltPerc", "clayPerc", "phosphorous", "potassium", "magnesium", "topsoilDepth", "NtotalConc")) %>%
  rename(variables = "...4") %>%
  bind_rows(eigenvals)
data <- as_tibble(values$sites[,1:3])
### Add to sites ###
sites <- sites %>%
  select(id, plot, calciumcarbonatPerc, NtotalPerc, cnRatio, pH, sandPerc, siltPerc, clayPerc, phosphorous, potassium, magnesium, topsoilDepth, NtotalConc) %>%
  group_by(plot) %>%
  summarise(across(where(is.numeric), ~median(.x, na.rm = T))) %>%
  select(plot) %>%
  bind_cols(data) %>%
  left_join(sites, by = "plot") %>%
  select(-position, -calciumcarbonatPerc, -NtotalPerc, -cnRatio, -pH, -sandPerc, -siltPerc, -clayPerc, -phosphorous, -potassium, -magnesium, -topsoilDepth, -NtotalConc)
rm(list=setdiff(ls(), c("sites", "species", "traits", "pca")))

### b Climate PCA  -------------------------------------------------------------------------------------------

### * Temperature ####
data <- read_csv(here("data/raw/temperature/data/data_OBS_DEU_P1M_T2M.csv"), col_names = T, na = c("", "NA", "na"), col_types = 
                          cols(
                            .default = "?"
                          )) %>%
  rename(date = Zeitstempel, value = Wert, site = SDO_ID) %>%
  filter(date >= "2002-03-01") %>%
  mutate(site = factor(site)) %>%
  select(site, date, value) %>%
  mutate(season = floor_date(date, "season")) %>%
  mutate(year = year(season), season = factor(month(season))) %>%
  mutate(season = fct_recode(season, "spring" = "3", "summer" = "6", "autumn" = "9", "winter" = "12")) %>%
  group_by(year) %>%
  mutate(yearMean = round(mean(value), digits = 1)) %>%
  mutate(currentYear = if_else(season == "spring", 0, 1)) %>%
  mutate(currentYear = year + currentYear) %>%
  group_by(currentYear) %>%
  mutate(currentMean = round(mean(value), digits = 1)) %>%
  mutate(currentMean = if_else(season == "spring", currentMean, NA_real_)) %>%
  group_by(year, season, yearMean, currentYear, currentMean) %>%
  summarise(seasonMean = round(mean(value), digits = 1)) %>%
  pivot_wider(id_cols = c(year, yearMean, currentMean), names_from = season, values_from = seasonMean) %>%
  group_by(year) %>%
  summarise(across(where(is.numeric), ~max(., na.rm = T))) %>% #warnings because of lates year (summer, autumn, winter)
  mutate(year = factor(year))
sites <- sites %>%
  right_join(data %>% select(-currentMean),
             by = c("constructionYearF" = "year")) %>%
  rename("TempSpring_cYear" = "spring", "TempSummer_cYear" = "summer", "TempAutumn_cYear" = "autumn", "TempWinter_cYear" = "winter", "TempMean_cYear" = "yearMean") %>%
  right_join(data %>% select(-currentMean), 
             by = c("constructionYearFplus" = "year")) %>%
  rename("TempSpring_cYearPlus" = "spring", "TempSummer_cYearPlus" = "summer", "TempAutumn_cYearPlus" = "autumn", "TempWinter_cYearPlus" = "winter", "TempMean_cYearPlus" = "yearMean") %>%
  right_join(data %>% select(year, currentMean, spring), 
             by = c("surveyYearF" = "year")) %>%
  rename("TempMean_sYear" = "currentMean", "TempSpring_sYear" = "spring") %>%
  right_join(data %>% select(year, summer, autumn, winter), 
             by = c("surveyYearFminus" = "year")) %>%
  rename("TempSummer_sYear" = "summer", "TempAutumn_sYear" = "autumn", "TempWinter_sYear" = "winter")

### * Precipitation ####
data <- read_csv(here("data/raw/precipitation/data/data_OBS_DEU_P1M_RR.csv"), col_names = T, na = c("", "NA", "na"), col_types = 
                          cols(
                            .default = "?"
                          )) %>%
  rename(date = Zeitstempel, value = Wert, site = SDO_ID) %>%
  filter(date >= "2002-03-01") %>%
  mutate(site = factor(site)) %>%
  select(site, date, value) %>%
  mutate(season = floor_date(date, "season")) %>%
  mutate(year = year(season), season = factor(month(season))) %>%
  mutate(season = fct_recode(season, "spring" = "3", "summer" = "6", "autumn" = "9", "winter" = "12")) %>%
  group_by(site, year) %>%
  mutate(yearSum = round(sum(value), digits = 0)) %>%
  mutate(currentYear = if_else(season == "spring", 0, 1)) %>%
  mutate(currentYear = year + currentYear) %>%
  group_by(site, currentYear) %>%
  mutate(currentSum = round(sum(value), digits = 0)) %>%
  group_by(site, season, year, currentYear, yearSum, currentSum) %>%
  summarise(seasonSum = round(sum(value), digits = 0)) %>%
  group_by(year, season, currentYear) %>%
  summarise(seasonMean = round(mean(seasonSum), digits = 0),
            yearMean = round(mean(yearSum), digits = 0),
            currentMean = round(mean(currentSum), digits = 0)) %>%
  mutate(currentMean = if_else(season == "spring", currentMean, NA_real_)) %>%
  pivot_wider(id_cols = c(year, yearMean, currentMean), names_from = season, values_from = seasonMean) %>%
  group_by(year) %>%
  summarise(across(where(is.numeric), ~max(., na.rm = T))) %>% #warnings because of lates year (summer, autumn, winter)
  mutate(year = factor(year))
sites <- sites %>%
  right_join(data %>% select(-currentMean),
             by = c("constructionYearF" = "year")) %>%
  rename("PrecSpring_cYear" = "spring", "PrecSummer_cYear" = "summer", "PrecAutumn_cYear" = "autumn", "PrecWinter_cYear" = "winter", "PrecMean_cYear" = "yearMean") %>%
  right_join(data %>% select(-currentMean), 
             by = c("constructionYearFplus" = "year")) %>%
  rename("PrecSpring_cYearPlus" = "spring", "PrecSummer_cYearPlus" = "summer", "PrecAutumn_cYearPlus" = "autumn", "PrecWinter_cYearPlus" = "winter", "PrecMean_cYearPlus" = "yearMean") %>%
  right_join(data %>% select(year, currentMean, spring), 
             by = c("surveyYearF" = "year")) %>%
  rename("PrecMean_sYear" = "currentMean", "PrecSpring_sYear" = "spring") %>%
  right_join(data %>% select(year, summer, autumn, winter), 
             by = c("surveyYearFminus" = "year")) %>%
  rename("PrecSummer_sYear" = "summer", "PrecAutumn_sYear" = "autumn", "PrecWinter_sYear" = "winter")

### * Calculation ####
### Prepare data ###
data <- sites %>%
  select(id, plot, 
         TempMean_sYear, TempSpring_sYear, TempSummer_sYear, TempAutumn_sYear, TempWinter_sYear,
         PrecMean_sYear, PrecSpring_sYear, PrecSummer_sYear, PrecAutumn_sYear, PrecWinter_sYear)
data <- sites %>%
  select(id, plot, 
         TempMean_cYear, TempSpring_cYear, TempSummer_cYear, TempAutumn_cYear, TempWinter_cYear,
         TempMean_cYearPlus, TempSpring_cYearPlus, TempSummer_cYearPlus, TempAutumn_cYearPlus, TempWinter_cYearPlus,
         PrecMean_cYear, PrecSpring_cYear, PrecSummer_cYear, PrecAutumn_cYear, PrecWinter_cYear,
         PrecMean_cYearPlus, PrecSpring_cYearPlus, PrecSummer_cYearPlus, PrecAutumn_cYearPlus, PrecWinter_cYearPlus) %>%
  group_by(plot) %>%
  summarise(across(where(is.numeric), ~median(.x, na.rm = T))) %>%
  select(-plot)
### Calculate PCA ###
pca <- rda(X = decostand(data, method = "standardize"), scale = T)
biplot(pca, display = "species")
screeplot(pca, bstick = TRUE, type = "l", main = NULL)
### Make data frames ###
eigenvals <- pca %>%
  eigenvals() %>%
  summary() %>%
  as_tibble() %>%
  select(PC1:PC3) %>%
  bind_cols(c("Eigenvalues", "Proportion Explained", "Cumulative Proportion")) %>%
  rename("variables" = "...4") %>%
  mutate(across(where(is.numeric), as.double))
values <- pca %>%
  summary()
pca <- values$species[ ,1:3] %>%
  as_tibble() %>%
  bind_cols(c("calciumcarbonatPerc", "NtotalPerc", "cnRatio", "pH", "sandPerc", "siltPerc", "clayPerc", "phosphorous", "potassium", "magnesium", "topsoilDepth", "NtotalConc")) %>%
  rename(variables = "...4") %>%
  bind_rows(eigenvals)
data <- as_tibble(values$sites[,1:3])
### Add to sites ###
sites2 <- sites %>%
  select(id, plot, calciumcarbonatPerc, NtotalPerc, cnRatio, pH, sandPerc, siltPerc, clayPerc, phosphorous, potassium, magnesium, topsoilDepth, NtotalConc) %>%
  group_by(plot) %>%
  summarise(across(where(is.numeric), ~median(.x, na.rm = T))) %>%
  select(plot) %>%
  bind_cols(data) %>%
  left_join(sites, by = "plot") %>%
  select(-starts_with(Temp), -starts_with(Prec))
rm(permdisp, data, nmds, names, values, eigenvals)


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# C Save processed data ##############################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

### Data
write_csv(sites, here("data/processed/data_processed_sites.csv"))
write_csv(species, here("data/processed/data_processed_species.csv"))
write_csv(species17, here("data/processed/data_processed_species17.csv"))
write_csv(species18, here("data/processed/data_processed_species18.csv"))
write_csv(species19, here("data/processed/data_processed_species19.csv"))
write_csv(traits, here("data/processed/data_processed_traits.csv"))

### Tables
setwd(here("data/tables"))
write_csv(pca, "table_pca.csv")
