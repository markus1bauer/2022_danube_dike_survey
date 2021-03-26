# Prepare vegetation data ####
# Markus Bauer



### Packages ###
library(here)
library(tidyverse)
library(vegan)
library(FD) #dbFD()
library(naniar) #are_na()

### Start ###
#installr::updateR(browse_news = F, install_R = T, copy_packages = T, copy_Rprofile.site = T, keep_old_packages = T, update_packages = T, start_new_R = F, quit_R = T, print_R_versions = T, GUI = F)
rm(list = ls())
setwd(here("data/raw"))



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Load data ##############################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


### 1 Species #####################################################################################

species <- read_csv2("data_raw_species.csv", col_names = F, na = c("", "NA", "na"))
### bring table in long format ###
names <- species %>% 
  slice(1) %>%
  pivot_longer(-(X1:X4), names_to = "id", values_to = "name") %>%
  select(-(X1:X4)) %>%
  mutate(name = str_replace_all(name, " ", "_"))
#names[which(duplicated(names$name)),]
#names[which(duplicated(names$id)),]
species <- species %>%
  select(-(X1:X3)) %>%
  slice(-c(1:2)) %>%
  pivot_longer(-(X4), names_to = "id", values_to = "abu") %>%
  rename(plot = X4) %>%
  pivot_wider(names_from = plot, names_prefix = "X", values_from = abu) %>%
  mutate_all(funs(str_replace(., ",", "."))) %>%
  mutate(id = as_factor(id)) %>%
  mutate_if(is_character, as.numeric) %>%
  select(id, where(~ is.numeric(.x) && sum(.x, na.rm = T) !=0 ))
species <- left_join(names, species, by = "id") %>% 
  select(-id) %>%
  mutate(name = as_factor(name))
rm(names)
### Check that each species occurs at least one time ###
species <- species %>%
  group_by(name) %>%
  mutate(total = sum(c_across(starts_with("X")), na.rm = T)) %>%
  mutate(presence = if_else(total > 0, 1, 0)) %>%
  filter(presence == 1) %>%
  ungroup() %>%
  select(-total, -presence) %>%
  mutate(name = as.character(name)) %>%
  arrange(name) %>%
  mutate(name = as_factor(name)) %>%
### Check that each plot has at least one species ###
  pivot_longer(-name, names_to = "id", values_to = "value") %>%
  group_by(id) %>%
  mutate(sum = sum(value, na.rm = T)) %>%
  filter(sum > 0) %>%
  select(-sum) %>%
  pivot_wider(names_from = id, values_from = value)
rm(list = setdiff(ls(), c("species"))) 

### 2 Sites #####################################################################################

sites <- read_csv2("data_raw_sites.csv", col_names = T, na = c("", "NA", "na"), col_types = 
                     cols(
                       .default = col_double(),
                       id = col_factor(),
                       block = col_factor(),
                       location = col_factor(),
                       side = col_factor(),
                       exposition = col_factor(),
                       ageCategory = col_factor(),
                       surveyDate_2017 = col_date(),
                       surveyDate_2018 = col_date(),
                       surveyDate_2019 = col_date(),
                       HCl = col_factor(),
                       humusLevel = col_factor(),
                       cnLevel = col_factor(),
                       phosphorousClass = col_factor(),
                       potassiumClass = col_factor(),
                       magnesiumClass = col_character()
                     )) %>%
  select(id, block, location, RW, HW, side, exposition, constructionYear, starts_with("vegetationCov"), starts_with("vegetationHeight"), calciumcarbonatPerc, humusPerc, NtotalPerc, cnRatio, pH, sandPerc, siltPerc, clayPerc, phosphorous, phosphorousClass, potassium, magnesium, topsoilDepth, fmDepth, sceletonRatiov, fmDbd, ufcPerc, ufc) %>%
  pivot_longer(starts_with("vegetationCov") | starts_with("vegetationHeight"), names_to = "surveyYear", values_to ="value") %>%
  separate(surveyYear, c("vegetationMeasure", "surveyYear"), "_", extra = "drop", fill = "right") %>%
  mutate(surveyYear = as.numeric(surveyYear)) %>%
  pivot_wider(names_from = "vegetationMeasure", values_from = "value") %>%
  mutate(id = str_c(id, surveyYear, sep = "_"), .keep = "all") %>%
  mutate(id = paste0("X", id)) %>%
  mutate(plot = str_sub(id, start = 2, end = 3)) %>%
  mutate(position = str_sub(id, start = 5, end = 5)) %>%
### Remove plots with no species ###
  filter(id %in% colnames(species[-1])) %>%
  mutate(id = as_factor(id)) %>%
  mutate(plot = as_factor(plot)) %>%
  mutate(position = as_factor(position))



### 3 Traits #####################################################################################

traits <- read_csv2("data_raw_traits.csv", col_names = T, na = c("", "NA", "na"), col_types = 
                      cols(
                        .default = col_factor(),
                        l = col_double(),
                        t = col_double(),
                        k = col_double(),
                        f = col_double(),
                        r = col_double(),
                        n = col_double()
                      )) %>%
#Check congruency of traits and species table
#traits$name[which(!(traits$name %in% species$name))]
#species$name[which(!(species$name %in% traits$name))]
  left_join(species, traits, by = "name") %>%
  select(name, abb, family, group, rlg, rlb, t, f, n, targetHerb, targetGrass, targetArrhenatherion, ffh6510, ffh6210, table30, table33, table34, nitrogenIndicator, nitrogenIndicator2, ruderalIndicator) %>%
  mutate(name = as.character(name)) %>%
  arrange(name) %>%
  mutate(name = as_factor(name))



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
  ))

sites <- sites %>%
  mutate(conf.low = c(1:length(id))) %>%
  mutate(conf.high = c(1:length(id))) %>%
  mutate(fmMass = fmDepth * fmDbd * 10) %>%
  mutate(NtotalConc = fmMass * NtotalPerc / 100) %>%
  mutate(plotAge = surveyYear - constructionYear) %>%
  select(-fmDepth, -fmMass)


### 2 Coverages #####################################################################################

### a Graminoid covratio; graminoid, herb, and total coverage)-------------------------------------------------------------------------------------------
sites <- species %>%
  mutate(type = traits$family) %>%
  pivot_longer(names_to = "id", values_to = "n", cols = starts_with("X")) %>%
  group_by(id, type) %>%
  summarise(total = sum(n, na.rm = T)) %>%
  mutate(type = if_else(type == "Poaceae" | type == "Cyperaceae" | type == "Juncaceae", "graminoidCov", "herbCov")) %>%
  group_by(id, type) %>%
  summarise(total = sum(total, na.rm = T)) %>%
  spread(type, total) %>%
  mutate(graminoidCovratio = graminoidCov / (graminoidCov + herbCov)) %>%
  mutate(accumulatedCov = graminoidCov + herbCov) %>%
  ungroup() %>%
  mutate(graminoidCovratio = round(graminoidCovratio, 3), .keep = "unused") %>%
  mutate(accumulatedCov = round(accumulatedCov, 3), .keep = "unused") %>%
  right_join(sites, by = "id")

### b Target species' coverage -------------------------------------------------------------------------------------------
sites <- species %>%
  mutate(type = traits$target) %>%
  pivot_longer(names_to = "id", values_to = "n", cols = starts_with("X")) %>%
  group_by(id, type) %>%
  summarise(total = sum(n, na.rm = T)) %>%
  filter(type == "yes") %>%
  select(-type) %>%
  ungroup() %>%
  mutate(targetCov = round(total, 3), .keep = "unused") %>%
  right_join(sites, by = "id")

### c Target herb species' coverage -------------------------------------------------------------------------------------------
sites <- species %>%
  mutate(type = traits$targetHerb) %>%
  pivot_longer(names_to = "id", values_to = "n", cols = starts_with("X")) %>%
  group_by(id, type) %>%
  summarise(total = sum(n, na.rm = T)) %>%
  filter(type == "yes") %>%
  select(-type) %>%
  ungroup() %>%
  mutate(targetHerbCov = round(total, 3), .keep = "unused") %>%
  right_join(sites, by = "id")

### d Arrhenatherum species' cover ratio -------------------------------------------------------------------------------------------
sites <- species %>%
  mutate(type = traits$targetArrhenatherion) %>%
  pivot_longer(names_to = "id", values_to = "n", cols = starts_with("X")) %>%
  group_by(id, type) %>%
  summarise(total = sum(n, na.rm = T)) %>%
  filter(type == "yes") %>%
  select(-type) %>%
  ungroup() %>%
  mutate(arrhCov = round(total, 3), .keep = "unused") %>%
  right_join(sites, by = "id")

### e Lean indicator's coverage -------------------------------------------------------------------------------------------
sites <- species %>%
  mutate(type = traits$leanIndicator) %>%
  pivot_longer(names_to = "id", values_to = "n", cols = starts_with("X")) %>%
  group_by(id, type) %>%
  summarise(total = sum(n, na.rm = T)) %>%
  filter(type == "yes") %>%
  select(-type) %>%
  ungroup() %>%
  mutate(leanCov = round(total, 3), .keep = "unused") %>%
  right_join(sites, by = "id")

### f Nitrogen indicator's coverage -------------------------------------------------------------------------------------------
sites <- species %>%
  mutate(type = traits$nitrogenIndicator) %>%
  pivot_longer(names_to = "id", values_to = "n", cols = starts_with("X")) %>%
  group_by(id, type) %>%
  summarise(total = sum(n, na.rm = T)) %>%
  filter(type == "yes") %>%
  select(-type) %>%
  ungroup() %>%
  mutate(nitrogenCov = round(total, 3), .keep = "unused") %>%
  right_join(sites, by = "id")

### g Ruderal indicator's coverage -------------------------------------------------------------------------------------------
sites <- species %>%
  mutate(type = traits$ruderalIndicator) %>%
  pivot_longer(names_to = "id", values_to = "n", cols = starts_with("X")) %>%
  group_by(id, type) %>%
  summarise(total = sum(n, na.rm = T)) %>%
  filter(type == "yes") %>%
  select(-type) %>%
  ungroup() %>%
  mutate(ruderalCov = round(total, 3), .keep = "unused") %>%
  right_join(sites, by = "id")

### h Table 33 species' coverage -------------------------------------------------------------------------------------------
sites <- species %>%
  mutate(type = traits$table33) %>%
  pivot_longer(names_to = "id", values_to = "n", cols = starts_with("X")) %>%
  group_by(id, type) %>%
  summarise(total = sum(n, na.rm = T)) %>%
  mutate(type = if_else(type == "4" | type =="3" | type == "2", "table33Cov", "other")) %>%
  group_by(id, type) %>%
  summarise(total = sum(total, na.rm = T)) %>%
  spread(type, total) %>%
  select(id, table33Cov) %>%
  ungroup() %>%
  mutate(table33Cov = round(table33Cov, 1)) %>%
  right_join(sites, by = "id")


### 3 Species richness #####################################################################################

specRich <- left_join(species, traits, by = "name") %>%
  select(rlg, rlb, target, targetHerb, targetArrhenatherion, ffh6510, ffh6210, nitrogenIndicator, leanIndicator, table33, table34, starts_with("X"))

### a total species richness -------------------------------------------------------------------------------------------
specRich_all <- specRich %>%
  pivot_longer(names_to = "id", values_to = "n", cols = starts_with("X")) %>%
  group_by(id) %>%
  mutate(n = if_else(n > 0, 1, 0)) %>%
  summarise(total = sum(n, na.rm = T)) %>%
  group_by(id) %>%
  summarise(speciesRichness = sum(total, na.rm = T)) %>%
  ungroup()

### b red list Germany (species richness) -------------------------------------------------------------------------------------------
specRich_rlg <- specRich %>%
  pivot_longer(names_to = "id", values_to = "n", cols = starts_with("X")) %>%
  group_by(id, rlg) %>%
  mutate(n = if_else(n > 0, 1, 0)) %>%
  summarise(total = sum(n, na.rm = T)) %>%
  filter(rlg == "1" | rlg == "2" | rlg == "3" | rlg == "V") %>%
  group_by(id) %>%
  summarise(rlgRichness = sum(total, na.rm = T)) %>%
  ungroup()

### c red list Bavaria (species richness) -------------------------------------------------------------------------------------------
specRich_rlb <- specRich %>%
  pivot_longer(names_to = "id", values_to = "n", cols = starts_with("X")) %>%
  group_by(id, rlb) %>%
  mutate(n = if_else(n > 0, 1, 0)) %>%
  summarise(total = sum(n, na.rm = T)) %>%
  filter(rlb == "1" | rlb == "2" | rlb == "3" | rlb == "V") %>%
  group_by(id) %>%
  summarise(rlbRichness = sum(total, na.rm = T)) %>%
  ungroup()

### d target species (species richness) -------------------------------------------------------------------------------------------
specRich_target <- specRich %>%
  pivot_longer(names_to = "id", values_to = "n", cols = starts_with("X")) %>%
  group_by(id, target) %>%
  mutate(n = if_else(n > 0, 1, 0)) %>%
  summarise(total = sum(n, na.rm = T)) %>%
  filter(target == "yes") %>%
  group_by(id) %>%
  summarise(targetRichness = sum(total, na.rm = T)) %>%
  ungroup()

### e target herb species (species richness) -------------------------------------------------------------------------------------------
specRich_targetHerb <- specRich %>%
  pivot_longer(names_to = "id", values_to = "n", cols = starts_with("X")) %>%
  group_by(id, targetHerb) %>%
  mutate(n = if_else(n > 0, 1, 0)) %>%
  summarise(total = sum(n, na.rm = T)) %>%
  filter(targetHerb == "yes") %>%
  group_by(id) %>%
  summarise(targetHerbRichness = sum(total, na.rm = T)) %>%
  ungroup()

### f Arrhenatherion species (species richness) -------------------------------------------------------------------------------------------
specRich_arrh <- specRich %>%
  pivot_longer(names_to = "id", values_to = "n", cols = starts_with("X")) %>%
  group_by(id, targetArrhenatherion) %>%
  mutate(n = if_else(n > 0, 1, 0)) %>%
  summarise(total = sum(n, na.rm = T)) %>%
  filter(targetArrhenatherion == "yes") %>%
  group_by(id) %>%
  summarise(arrhRichness = sum(total, na.rm = T)) %>%
  ungroup()

### g ffh6510 species (species richness) -------------------------------------------------------------------------------------------
specRich_ffh6510 <- specRich %>%
  pivot_longer(names_to = "id", values_to = "n", cols = starts_with("X")) %>%
  group_by(id, ffh6510) %>%
  mutate(n = if_else(n > 0, 1, 0)) %>%
  summarise(total = sum(n, na.rm = T)) %>%
  filter(ffh6510 == "yes") %>%
  group_by(id) %>%
  summarise(ffh6510Richness = sum(total, na.rm = T)) %>%
  ungroup()

### h ffh6210 species (species richness) -------------------------------------------------------------------------------------------
specRich_ffh6210 <- specRich %>%
  pivot_longer(names_to = "id", values_to = "n", cols = starts_with("X")) %>%
  group_by(id, ffh6210) %>%
  mutate(n = if_else(n > 0, 1, 0)) %>%
  summarise(total = sum(n, na.rm = T)) %>%
  filter(ffh6210 == "yes") %>%
  group_by(id) %>%
  summarise(ffh6210Richness = sum(total, na.rm = T)) %>%
  ungroup()

### i leanIndicator species (species richness) -------------------------------------------------------------------------------------------
specRich_leanIndicator <- specRich %>%
  pivot_longer(names_to = "id", values_to = "n", cols = starts_with("X")) %>%
  group_by(id, leanIndicator) %>%
  mutate(n = if_else(n > 0, 1, 0)) %>%
  summarise(total = sum(n, na.rm = T)) %>%
  filter(leanIndicator == "yes") %>%
  group_by(id) %>%
  summarise(leanIndicatorRichness = sum(total, na.rm = T)) %>%
  ungroup()

### k table33 (species richness) -------------------------------------------------------------------------------------------
specRich_table33_2 <- specRich %>%
  pivot_longer(names_to = "id", values_to = "n", cols = starts_with("X")) %>%
  group_by(id, table33) %>%
  mutate(n = if_else(n > 0, 1, 0)) %>%
  summarise(total = sum(n, na.rm = T)) %>%
  filter(table33 == "2") %>%
  group_by(id) %>%
  summarise(table33_2Richness = sum(total, na.rm = T)) %>%
  ungroup()
specRich_table33_3 <- specRich %>%
  pivot_longer(names_to = "id", values_to = "n", cols = starts_with("X")) %>%
  group_by(id, table33) %>%
  mutate(n = if_else(n > 0, 1, 0)) %>%
  summarise(total = sum(n, na.rm = T)) %>%
  filter(table33 == "3") %>%
  group_by(id) %>%
  summarise(table33_3Richness = sum(total, na.rm = T)) %>%
  ungroup()
specRich_table33_4 <- specRich %>%
  pivot_longer(names_to = "id", values_to = "n", cols = starts_with("X")) %>%
  group_by(id, table33) %>%
  mutate(n = if_else(n > 0, 1, 0)) %>%
  summarise(total = sum(n, na.rm = T)) %>%
  filter(table33 == "4") %>%
  group_by(id) %>%
  summarise(table33_4Richness = sum(total, na.rm = T)) %>%
  ungroup()

### l table34 (species richness) -------------------------------------------------------------------------------------------
specRich_table34_2 <- specRich %>%
  pivot_longer(names_to = "id", values_to = "n", cols = starts_with("X")) %>%
  group_by(id, table34) %>%
  mutate(n = if_else(n > 0, 1, 0)) %>%
  summarise(total = sum(n, na.rm = T)) %>%
  filter(table34 == "2") %>%
  group_by(id) %>%
  summarise(table34_2Richness = sum(total, na.rm = T)) %>%
  ungroup()
specRich_table34_3 <- specRich %>%
  pivot_longer(names_to = "id", values_to = "n", cols = starts_with("X")) %>%
  group_by(id, table34) %>%
  mutate(n = if_else(n > 0, 1, 0)) %>%
  summarise(total = sum(n, na.rm = T)) %>%
  filter(table34 == "3") %>%
  group_by(id) %>%
  summarise(table34_3Richness = sum(total, na.rm = T)) %>%
  ungroup()

### m implement in sites data set -------------------------------------------------------------------------------------------
sites <- sites %>%
  right_join(specRich_all, by = "id") %>%
  right_join(specRich_rlg, by = "id") %>%
  right_join(specRich_rlb, by = "id") %>%
  right_join(specRich_target, by = "id") %>%
  right_join(specRich_targetHerb, by = "id") %>%
  right_join(specRich_arrh, by = "id") %>%
  right_join(specRich_ffh6510, by = "id") %>%
  right_join(specRich_ffh6210, by = "id") %>%
  right_join(specRich_leanIndicator, by = "id") %>%
  right_join(specRich_table33_2, by = "id") %>%
  right_join(specRich_table33_3, by = "id") %>%
  right_join(specRich_table33_4, by = "id") %>%
  right_join(specRich_table34_2, by = "id") %>%
  right_join(specRich_table34_3, by = "id") %>%
### Calcute the ratio of target species richness of total species richness
  mutate(targetRichratio = targetRichness / speciesRichness) %>%
  mutate(targetRichratio = round(targetRichratio, 3))
rm(list=setdiff(ls(), c("sites", "species", "traits")))


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


### 5 TBI #####################################################################################

tbi <- species %>%
  pivot_longer(-name, "id", "n") %>%
  pivot_wider(id, name)
tbi <- left_join(tbi, sites, by = "id") %>%
  select(id, exposition, surveyYear, Achillea_millefolium:Rubus_fruticosus_agg) %>%
  filter(surveyYear == 2017 | surveyYear == 2019) %>%
  #tbiStart <- sites %>% select(id, plot, surveyYear) %>% filter(surveyYear == 2017)
  #tbiEnd <- sites %>% select(id, plot, surveyYear) %>% filter(surveyYear == 2019)
  #anti_join(tbiStart, tbiEnd, by = "plot")
  #anti_join(tbiEnd, tbiStart, by = "plot")
  #rm(tbiStart, tbiEnd)
  filter(id != "X15_m_2017" &
           id != "X16_m_2017" &
           id != "X47_m_2017" &
           id != "X48_m_2017" &
           id != "X49_m_2017" &
           id != "X50_m_2017" &
           id != "X34_o_2019") %>%
  column_to_rownames(var = "id") %>%
  select(-surveyYear)

### a Abundance data --------------------------------------------------------------------------------------------
tbi2 <- select(tbi)
tbiAbu <- tbi2[,colSums(tbi2) > 0]

### b Presence-absence data--------------------------------------------------------------------------------------------
tbiPa <- tbiAbu
tbiPa[tbiPa > 0] = 1

### c Abundance Exposition--------------------------------------------------------------------------------------------
#tbi2 <- tbi %>%
  #filter(exposition == "north") %>%
 # select(-exposition)
#tbiAbuN <- tbi2[,colSums(tbi2) > 0]

tbiAbu <- rownames_to_column(tbiAbu, var = "id")
tbiPa <- rownames_to_column(tbiPa, var = "id")
#tbiAbuN <- rownames_to_column(tbiAbuN, var = "id")
rm(tbi, tbi2)

### 6 CWM of Ellenberg #####################################################################################

### a N value -------------------------------------------------------------------------------------------
Ntraits <- traits %>%
  select(name, n) %>%
  filter(n > 0) 
Nspecies <- semi_join(species, Ntraits, by = "name") %>%
  pivot_longer(-name, "site", "value") %>%
  pivot_wider(site, name) %>%
  column_to_rownames("site")
Ntraits <- column_to_rownames(Ntraits, "name")
### Calculate CWM ###
Nweighted <- dbFD(Ntraits, Nspecies, w.abun = T,
                        calc.FRic = F, calc.FDiv = F, corr = "sqrt")

### b F value -------------------------------------------------------------------------------------------
Ftraits <- traits %>%
  select(name, f) %>%
  filter(f > 0) 
Fspecies <- semi_join(species, Ftraits, by = "name") %>%
  pivot_longer(-name, "site", "value") %>%
  pivot_wider(site, name) %>%
  column_to_rownames("site")
Ftraits <- column_to_rownames(Ftraits, "name")
### Calculate CWM ###
Fweighted <- dbFD(Ftraits, Fspecies, w.abun = T,
                  calc.FRic = F, calc.FDiv = F, corr = "sqrt")

### c T value -------------------------------------------------------------------------------------------
Ttraits <- traits %>%
  select(name, t) %>%
  filter(t > 0) 
Tspecies <- semi_join(species, Ttraits, by = "name") %>%
  pivot_longer(-name, "site", "value") %>%
  pivot_wider(site, name) %>%
  column_to_rownames("site")
Ttraits <- column_to_rownames(Ttraits, "name")
### Calculate CWM ###
Tweighted <- dbFD(Ttraits, Tspecies, w.abun = T,
                  calc.FRic = F, calc.FDiv = F, corr = "sqrt")

### d implement in sites data set -------------------------------------------------------------------------------------------
sites$cwmAbuN <- round(as.numeric(as.character(Nweighted$CWM$n)), 3)
sites$cwmAbuF <- round(as.numeric(as.character(Fweighted$CWM$f)), 3)
sites$cwmAbuT <- round(as.numeric(as.character(Tweighted$CWM$t)), 3)
rm(list=setdiff(ls(), c("sites", "species", "traits", "tbiAbu", "tbiPa")))


### 7 Functional plant traits #####################################################################################

#### * read LEDA data #### 
dataSLA <- data.table::fread("data_raw_LEDA_20210223_sla.txt", 
                             sep = ";",
                             dec = ".",
                             skip = 3,
                             header = T,
                             select = c("SBS name", "single value [mm^2/mg]"),
                             ) %>%
  as_tibble() %>%
  rename(name = "SBS name") %>%
  rename(sla = "single value [mm^2/mg]") %>%
  mutate(name = as_factor(str_replace_all(name, " ", "_")))

dataSM <- data.table::fread("data_raw_LEDA_20210223_seedmass.txt", 
                            sep = ";",
                            dec = ".",
                            skip = 3,
                            header = T,
                            select = c("SBS name", "single value [mg]"),
                            ) %>%
  as_tibble() %>%
  rename(name = "SBS name") %>%
  rename(seedmass = "single value [mg]") %>%
  mutate(name = as_factor(str_replace_all(name, " ", "_")))

dataH <- data.table::fread("data_raw_LEDA_20210223_canopy_height.txt", 
                           sep = ";",
                           dec = ".",
                           skip = 3,
                           header = T,
                           select = c("SBS name", "single value [m]")
                             ) %>%
  as_tibble() %>%
  rename(name = "SBS name") %>%
  rename(height = "single value [m]") %>%
  mutate(name = as_factor(str_replace_all(name, " ", "_")))


data <- full_join(dataSLA, dataSM, by = "name") %>%
  full_join(dataH, by = "name")
rm(dataSLA, dataSM, dataH)
##### Synonyms ###
#traits$name[which(!(traits$name %in% data$name))]
#data %>%
  #group_by(name) %>%
  #summarise(across(where(is.double), ~median(.x, na.rm = T))) %>%
  #filter(str_detect(name, "incana"))
traits <- data %>%
  mutate(name = fct_recode(name, Centaurea_stoebe = "Centaurea_stoebe_s.lat.")) %>%
  mutate(name = fct_recode(name, Carex_praecox_ssp_praecox = "Carex_praecox")) %>%
  mutate(name = fct_recode(name, Cerastium_fontanum_ssp_vulgare = "Cerastium_fontanum")) %>%
  mutate(name = fct_recode(name, Cerastium_fontanum_ssp_vulgare = "Cerastium_fontanum_s._vulgare")) %>%
  mutate(name = fct_recode(name, Cota_tinctoria = "Anthemis_tinctoria")) %>%
  mutate(name = fct_recode(name, Cyanus_segetum = "Centaurea_cyanus")) %>%
  mutate(name = fct_recode(name, Euphorbia_verrucosa = "Euphorbia_brittingeri")) %>%
  mutate(name = fct_recode(name, Helictotrichon_pubescens = "Avenula_pubescens")) %>%
  mutate(name = fct_recode(name, Hypochaeris_radicata = "Hypochoeris_radicata")) %>%
  mutate(name = fct_recode(name, Jacobaea_vulgaris = "Senecio_vulgaris")) %>%
  mutate(name = fct_recode(name, Medicago_falcata = "Medicago_sativa_s._falcata")) %>%
  mutate(name = fct_recode(name, Ononis_spinosa_ssp_procurrens = "Ononis_repens")) %>%
  mutate(name = fct_recode(name, Persicaria_amphibia = "Polygonum_amphibium")) %>%
  mutate(name = fct_recode(name, Pilosella_caespitosa = "Hieracium_caespitosum")) %>%
  mutate(name = fct_recode(name, Pilosella_officinarum = "Hieracium_pilosella")) %>%
  mutate(name = fct_recode(name, Pilosella_piloselloides = "Hieracium_piloselloides")) %>%
  mutate(name = fct_recode(name, Plantago_major_ssp_intermedia = "Plantago_major_subsp._intermedia")) %>%
  mutate(name = fct_recode(name, Rubus_fruticosus_agg = "Rubus_fruticosus")) %>%
  mutate(name = fct_recode(name, Rubus_fruticosus_agg = "Rubus_fruticosus_ag._L.")) %>%
  mutate(name = fct_recode(name, Securigera_varia = "Coronilla_varia")) %>%
  mutate(name = fct_recode(name, "Silene_flos-cuculi" = "Lychnis_flos-cuculi")) %>%
  mutate(name = fct_recode(name, Taraxacum_campylodes = "Taraxacum_officinale")) %>%
  mutate(name = fct_recode(name, Taraxacum_campylodes = "Taraxacum_Sec._Ruderalia")) %>%
  mutate(name = fct_recode(name, Tripleurospermum_maritimum = "Matricaria_maritima")) %>%
  mutate(name = fct_recode(name, Vicia_sativa_ssp_nigra= "Vicia_sativa_s._nigra")) %>%  
  group_by(name) %>%
  summarise(across(where(is.double), ~median(.x, na.rm = T))) %>%
  right_join(traits, by = "name")
### check completeness of LEDA ###
#(test2 <- traits %>%
  #select(name, sla, seedmass, height) %>%
  #filter(!complete.cases(.)))

### * read TRY data ####
data <- data.table::fread("data_raw_TRY_20210306_13996.txt", 
              header = T, 
              sep = "\t", 
              dec = ".", 
              quote = "") %>%
  as_tibble() %>%
  rename(name = "AccSpeciesName") %>%
  rename(value = "StdValue") %>%
  rename(trait = "TraitID") %>%
  select(name, value, trait) %>%
  mutate(name = as_factor(str_replace_all(name, " ", "_"))) %>%
  drop_na %>%
  mutate(trait = str_replace(trait, "26", "seedmass")) %>%
  mutate(trait = str_replace(trait, "3106", "height")) %>%
  mutate(trait = str_replace(trait, "3107", "height")) %>%
  mutate(trait = str_replace(trait, "3115", "sla")) %>%
  mutate(trait = str_replace(trait, "3116", "sla")) %>%
  mutate(trait = str_replace(trait, "3117", "sla"))
##### Synonyms ###
#test2$name[which(!(test2$name %in% data$name))]
#data %>%
  #group_by(name) %>%
  #summarise(across(where(is.double), ~median(.x, na.rm = T))) %>%
  #filter(str_detect(name, "Equisetum"))
traits <- data %>%
  mutate(name = fct_recode(name, Carex_praecox_ssp_praecox = "Carex_praecox")) %>%
  mutate(name = fct_recode(name, Plantago_major_ssp_intermedia = "Plantago_major_subsp._intermedia")) %>%
  mutate(name = fct_recode(name, Ranunculus_serpens_ssp_nemorosus = "Ranunculus_serpens_subsp._nemorosus")) %>%
  mutate(trait = as_factor(trait)) %>%
  group_by(name, trait) %>%
  summarise(across(where(is.double), ~median(.x, na.rm = T))) %>%
  spread(trait, value) %>%
  right_join(traits, by = "name") %>%
  mutate(sla = coalesce(sla.x, sla.y), .keep = "unused") %>%
  mutate(seedmass = coalesce(seedmass.x, seedmass.y), .keep = "unused") %>%
  mutate(height = coalesce(height.x, height.y), .keep = "unused")
### check completeness of LEDA + TRY ###
#(test2 <- traits %>%
  #select(name, sla, seedmass, height) %>%
  #filter(!complete.cases(.)))

### * read GrooT data ####
data <- read.csv("data_raw_GrooT.csv", header = T, na.strings = c("", "NA")) %>%
  filter(traitName == "Specific_root_length" |
           traitName == "Root_length_density_volume" |
           traitName == "Root_mass_fraction" |
           traitName == "Lateral_spread") %>%
  mutate(name = as_factor(str_c(genusTNRS, speciesTNRS, sep = "_"))) %>%
  rename(trait = traitName, value = medianSpecies) %>%
  select(name, trait, value) %>%
  pivot_wider(names_from = "trait", values_from = "value") %>%
  rename(rmf = Root_mass_fraction, rld = Root_length_density_volume, srl = Specific_root_length, lateral = Lateral_spread)
##### Find Synonyms ###
#traits$name[which(!(traits$name %in% data$name))]
#data %>%
#group_by(name) %>%
#summarise(across(where(is.double), ~median(.x, na.rm = T))) %>%
#filter(str_detect(name, "Angelica"))
traits <- data %>%
  mutate(name = fct_recode(name, Cota_tinctoria = "Anthemis_tinctoria")) %>%
  mutate(name = fct_recode(name, Carex_praecox_ssp_praecox = "Carex_praecox")) %>%
  mutate(name = fct_recode(name, Cerastium_fontanum_ssp_vulgare = "Cerastium_fontanum")) %>%
  mutate(name = fct_recode(name, Plantago_major_ssp_major = "Plantago_major")) %>%
  mutate(name = fct_recode(name, Jacobaea_vulgaris = "Senecio_jacobaea")) %>%
  mutate(name = fct_recode(name, Ranunculus_serpens_ssp_nemorosus = "Ranunculus_nemorosus")) %>%
  mutate(name = fct_recode(name, "Silene_flos-cuculi" = "Lychnis_flos-cuculi")) %>%
  mutate(name = fct_recode(name, Silene_latifolia_ssp_alba = "Silene_latifolia")) %>%
  mutate(name = fct_recode(name, Vicia_sativa_ssp_nigra = "Vicia_sativa")) %>%
  group_by(name) %>%
  summarise(across(where(is.double), ~median(.x, na.rm = T))) %>%
  right_join(traits, by = "name")
### check completeness of Roots ###
(test2 <- traits %>%
    select(name, rmf, rld, srl, lateral) %>%
    filter(!complete.cases(rmf, rld, srl, lateral)))

### * check completeness of traits ####
test <- traits %>%
  select(t, n, f, sla, seedmass, height, rmf, rld, srl, lateral)
vis_miss(test, cluster = F, sort_miss = T)
#gg_miss_var(test)
gg_miss_case(test, order_cases = F)
(test2 <- traits %>%
    select(name, sla, seedmass, height, rmf, rld, srl, lateral) %>%
    filter(!complete.cases(sla, seedmass, height, rmf, rld, srl, lateral)))
rm(test, test2, data)

### * prepare data frames ####
species <- species %>%
  mutate(name = as.character(name)) %>%
  arrange(name) %>%
  mutate(name = as_factor(name))
traits <- traits %>%
  mutate(name = as.character(name)) %>%
  arrange(name) %>%
  mutate(name = as_factor(name))
herbCount <- traits %>%
  filter(group != "tree" & group != "shrub") %>%
  left_join(species, by = "name") %>%
  count() %>%
  pull()
undefinedSpeciesCount <- traits %>%
  filter(group != "tree" & group != "shrub") %>%
  filter(str_detect(name, "_spec")) %>%
  left_join(species, by = "name") %>%
  count() %>%
  pull()
traitsLHS <- traits %>%
  filter(group != "tree" & group != "shrub") %>%
  select(name, sla, seedmass, height) %>%
  drop_na()
traitsSLA <- traits %>%
  filter(group != "tree" & group != "shrub") %>%
  select(name, sla) %>%
  drop_na()  
traitsSM <- traits %>%
  filter(group != "tree" & group != "shrub") %>%
  select(name, seedmass) %>%
  drop_na()
traitsH <- traits %>%
  filter(group != "tree" & group != "shrub") %>%
  select(name, height) %>%
  drop_na()
traitsSRL <- traits %>%
  filter(group != "tree" & group != "shrub") %>%
  select(name, srl) %>%
  drop_na()
traitsRLD <- traits %>%
  filter(group != "tree" & group != "shrub") %>%
  select(name, rld) %>%
  drop_na()
traitsRMF <- traits %>%
  filter(group != "tree" & group != "shrub") %>%
  select(name, rmf) %>%
  drop_na()
traitsAll <- traits %>%
  filter(group != "tree" & group != "shrub") %>%
  select(name, sla, seedmass, height, srl) %>%
  drop_na()


### 8 CWM and FDis of functional plant traits #####################################################################################

### a LHS -------------------------------------------------------------------------------------------
Tspecies <- semi_join(species, traitsLHS, by = "name")
Ttraits <- semi_join(traitsLHS, Tspecies, by = "name")
Tspecies <- Tspecies %>%
  pivot_longer(-name, "site", "value") %>%
  pivot_wider(site, name) %>%
  column_to_rownames("site")
Ttraits <- column_to_rownames(Ttraits, "name")
log_Ttraits <- log(Ttraits)
TdiversityAbu <- dbFD(log_Ttraits, Tspecies, w.abun = T,
               calc.FRic = F, calc.FDiv = F, corr = "cailliez")
sites$fdisAbuLHS <- TdiversityAbu$FDis

### b SLA -------------------------------------------------------------------------------------------
Tspecies <- semi_join(species, traitsSLA, by = "name")
Ttraits <- semi_join(traitsSLA, Tspecies, by = "name")
Tspecies <- Tspecies %>%
  pivot_longer(-name, "site", "value") %>%
  pivot_wider(site, name) %>%
  column_to_rownames("site")
Ttraits <- column_to_rownames(Ttraits, "name")
log_Ttraits <- log(Ttraits)
TdiversityAbu <- dbFD(log_Ttraits, Tspecies, w.abun = T, 
                   calc.FRic = F, calc.FDiv = F, corr = "sqrt");
sites$fdisAbuSla <- TdiversityAbu$FDis
sites$cwmAbuSla <- exp(as.numeric(as.character(TdiversityAbu$CWM$sla)))

### c Seed mass -------------------------------------------------------------------------------------------
Tspecies <- semi_join(species, traitsSM, by = "name")
Ttraits <- semi_join(traitsSM, Tspecies, by = "name")
Tspecies <- Tspecies %>%
  pivot_longer(-name, "site", "value") %>%
  pivot_wider(site, name) %>%
  column_to_rownames("site")
Ttraits <- column_to_rownames(Ttraits, "name")
log_Ttraits <- log(Ttraits)
TdiversityAbu <- dbFD(log_Ttraits, Tspecies, w.abun = T, 
                      calc.FRic = F, calc.FDiv = F, corr = "sqrt");
sites$fdisAbuSeedmass <- TdiversityAbu$FDis
sites$cwmAbuSeedmass <- exp(as.numeric(as.character(TdiversityAbu$CWM$seedmass)))

### d Canopy height -------------------------------------------------------------------------------------------
Tspecies <- semi_join(species, traitsH, by = "name")
Ttraits <- semi_join(traitsH, Tspecies, by = "name")
Tspecies <- Tspecies %>%
  pivot_longer(-name, "site", "value") %>%
  pivot_wider(site, name) %>%
  column_to_rownames("site")
Ttraits <- column_to_rownames(Ttraits, "name")
log_Ttraits <- log(Ttraits)
TdiversityAbu <- dbFD(log_Ttraits, Tspecies, w.abun = T, 
                      calc.FRic = F, calc.FDiv = F, corr = "sqrt");
sites$fdisAbuHeight <- TdiversityAbu$FDis
sites$cwmAbuHeight <- exp(as.numeric(as.character(TdiversityAbu$CWM$height)))
rm(TdiversityAbu, log_Ttraits, Ttraits, Tspecies, traitsLHS, traitsSLA, traitsSM, traitsH)
### e Specific root length -------------------------------------------------------------------------------------------
Tspecies <- semi_join(species, traitsSRL, by = "name")
Ttraits <- semi_join(traitsSRL, Tspecies, by = "name")
Tspecies <- Tspecies %>%
  pivot_longer(-name, "site", "value") %>%
  pivot_wider(site, name) %>%
  column_to_rownames("site")
Ttraits <- column_to_rownames(Ttraits, "name")
log_Ttraits <- log(Ttraits)
TdiversityAbu <- dbFD(log_Ttraits, Tspecies, w.abun = T, 
                      calc.FRic = F, calc.FDiv = F, corr = "sqrt");
sites$fdisAbuSrl <- TdiversityAbu$FDis
sites$cwmAbuSrl <- TdiversityAbu$CWM$srl %>%
  as.character() %>% as.numeric() %>% exp()
length(traitsSRL$name) / (herbCount - undefinedSpeciesCount)

### f Root mass fraction -------------------------------------------------------------------------------------------
Tspecies <- semi_join(species, traitsRMF, by = "name")
Ttraits <- semi_join(traitsRMF, Tspecies, by = "name")
Tspecies <- Tspecies %>%
  pivot_longer(-name, "site", "value") %>%
  pivot_wider(site, name) %>%
  column_to_rownames("site")
Ttraits <- column_to_rownames(Ttraits, "name")
log_Ttraits <- log(Ttraits)
TdiversityAbu <- dbFD(log_Ttraits, Tspecies, w.abun = T, 
                      calc.FRic = F, calc.FDiv = F, corr = "sqrt");
sites$fdisAbuRmf <- TdiversityAbu$FDis
sites$cwmAbuRmf <- TdiversityAbu$CWM$rmf %>%
  as.character() %>% as.numeric() %>% exp()
length(traitsRMF$name) / (herbCount - undefinedSpeciesCount)

### g All -------------------------------------------------------------------------------------------
Tspecies <- semi_join(species, traitsAll, by = "name")
Ttraits <- semi_join(traitsAll, Tspecies, by = "name")
Tspecies <- Tspecies %>%
  pivot_longer(-name, "site", "value") %>%
  pivot_wider(site, name) %>%
  column_to_rownames("site")
Ttraits <- column_to_rownames(Ttraits, "name")
log_Ttraits <- log(Ttraits)
TdiversityAbu <- dbFD(log_Ttraits, Tspecies, w.abun = T,
                      calc.FRic = F, calc.FDiv = F, corr = "cailliez")
sites$fdisAbuAll <- TdiversityAbu$FDis
length(traitsAll$name) / (herbCount - undefinedSpeciesCount)
rm(list=setdiff(ls(), c("sites", "species", "traits", "tbiPa", "tbiAbu")))


### 9 Ordination #####################################################################################

### a NMDS -------------------------------------------------------------------------------------------
nmds <- species %>%
  pivot_longer(-name, "id", "value") %>%
  pivot_wider(id, name) %>%
  arrange(id) %>%
  column_to_rownames("id")
set.seed(1)
(ordi <- metaMDS(nmds, try = 99, previous.best = T, na.rm = T))
sites <- ordi %>%
  scores() %>%
  as.data.frame() %>%
  rownames_to_column(var = "id") %>%
  as_tibble() %>%
  select(id, NMDS1, NMDS2) %>%
  right_join(sites, by = "id")

### b PERMDISP -------------------------------------------------------------------------------------------
(permdisp <- betadisper(d = vegdist(nmds), group = sites$plot))
sites <- mutate(sites, permdisp = as.numeric(permdisp$distances))
rm(permdisp, ordi, nmds)
  
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# C Save processed data ##############################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


setwd(here("data/processed"))
write_csv2(sites, "data_processed_sites.csv")
write_csv2(species, "data_processed_species.csv")
write_csv2(traits, "data_processed_traits.csv")
write_csv2(tbiAbu, "data_processed_tbiAbu.csv")
write_csv2(tbiPa, "data_processed_tbiPa.csv")
