# Prepare vegetation data ####
# Markus Bauer
# Citation: Markus Bauer, Jakob Huber, Johannes Kollmann
# DOI


### Packages ###
library(tidyverse)
library(vegan)
library(FD) #dbFD()
library(naniar) #are_na()

### Start ###
#installr::updateR(browse_news = F, install_R = T, copy_packages = T, copy_Rprofile.site = T, keep_old_packages = T, update_packages = T, start_new_R = F, quit_R = T, print_R_versions = T, GUI = F)
rm(list = ls())
setwd("Z:/Documents/0_Donaudeiche/3_Aufnahmen_und_Ergebnisse/2022_Danube_old_dikes/data/raw")



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Load data ##############################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


### 1 Species #####################################################################################

species <- read_csv2("data_raw_species.csv", col_names = F, na = c("", "NA", "na"))
### bring table in long format ###
names <- species %>% 
  slice(1) %>%
  pivot_longer(-(X1:X4), names_to = "id", values_to = "name") %>%
  select(-(X1:X4))
names$name <- str_replace_all(names$name, " ", "_")
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
  select(-total, -presence)
 

### 2 Sites #####################################################################################

sites <- read_csv2("data_raw_sites.csv", col_names = T, na = c("", "NA", "na"), col_types = 
                     cols(
                       .default = col_double(),
                       id = col_factor(),
                       pair = col_factor(),
                       location = col_factor(),
                       side = col_factor(),
                       exposition = col_factor(),
                       ageCategory = col_factor(),
                       HCl = col_factor(),
                       humusLevel = col_factor(),
                       cnLevel = col_factor(),
                       phosphorousClass = col_factor(),
                       potassiumClass = col_factor(),
                       magnesiumClass = col_character()
                     )        
)
sites <- sites %>%
  select(id, starts_with("survey"), pair, location, RW, HW, side, exposition, constructionYear, calciumcarbonatPerc, humusPerc, NtotalPerc, cnRatio, pH, sandPerc, siltPerc, clayPerc, phosphorous, phosphorousClass, potassium, magnesium, topsoilDepth, fmDepth, sceletonRatiov, fmDbd, ufcPerc, ufc) %>%
  pivot_longer(starts_with("survey"), names_to = "survey", values_to ="surveyYear") %>%
  select(-survey) %>%
  mutate(id = str_c(id, surveyYear, sep = "_"), .keep = "all") %>%
  mutate(id = paste0("X", id)) %>%
  mutate(plot = str_sub(id, start = 2, end = 3)) %>%
  mutate(position = str_sub(id, start = 5, end = 5))
### Remove plots with no species ###
ids <- colnames(species[-1])
sites <- sites %>%
  filter(id %in% ids) %>%
  mutate(id = as_factor(id)) %>%
  mutate(plot = as_factor(plot)) %>%
  mutate(position = as_factor(position))
rm(ids)


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
                      )
                   )
#Check congruency of traits and species table
#traits$name[which(!(traits$name %in% species$name))]
#species$name[which(!(species$name %in% traits$name))]
traits <- left_join(species, traits, by = "name") %>%
  select(name, abb, family, rlg, rlb, t, f, n, targetHerb, targetGrass, targetArrhenatherion, ffh6510, ffh6210, table30, table33, table34, nitrogenIndicator, nitrogenIndicator2, ruderalIndicator)

### Check missing values ###
#n_miss(sites)
#n_miss(species)
#miss_var_summary(species, order = T)
#vis_miss(species, cluster = F)



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Create variables ##############################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


### 1 Create simple variables #####################################################################################

sites <- sites %>%
  mutate(conf.low = c(1:length(id))) %>%
  mutate(conf.high = c(1:length(id))) %>%
  mutate(fmMass = fmDepth * fmDbd * 10) %>%
  mutate(NtotalConc = fmMass * NtotalPerc / 100) %>%
  select(-fmDepth, -fmMass)

traits <- traits %>%
  mutate(leanIndicator = if_else(
    !(is.na(table30)) | !(is.na(table33)) | !(is.na(table34)), "yes", "no"
  )) %>%
  mutate(target = if_else(
    targetHerb == "yes" | targetGrass == "yes", "yes", "no"
  ))


### 2 Coverages #####################################################################################

### a Graminoid covratio; graminoid, herb, and total coverage)-------------------------------------------------------------------------------------------
data <- species %>%
  mutate(type = traits$family) %>%
  pivot_longer(names_to = "id", values_to = "n", cols = starts_with("X")) %>%
  group_by(id, type) %>%
  summarise(total = sum(n)) %>%
  mutate(type = if_else(type == "Poaceae" | type == "Cyperaceae" | type == "Juncaceae", "graminoidCov", "herbCov")) %>%
  group_by(id, type) %>%
  summarise(total = sum(total)) %>%
  spread(type, total) %>%
  mutate(graminoidCovratio = graminoidCov / (graminoidCov + herbCov)) %>%
  mutate(accumulatedCov = graminoidCov + herbCov) %>%
  ungroup() %>%
  mutate(graminoidCovratio = round(graminoidCovratio, 3), .keep = "unused") %>%
  mutate(accumulatedCov = round(accumulatedCov, 3), .keep = "unused")
sites <- left_join(sites, data, by = "id")

### b Target species' coverage -------------------------------------------------------------------------------------------
data <- species %>%
  mutate(type = traits$target) %>%
  pivot_longer(names_to = "id", values_to = "n", cols = starts_with("X")) %>%
  group_by(id, type) %>%
  summarise(total = sum(n)) %>%
  filter(type == "yes") %>%
  select(-type) %>%
  ungroup() %>%
  mutate(targetCov = round(total, 3), .keep = "unused")
sites <- left_join(sites, data, by = "id")

### c Target herb species' coverage -------------------------------------------------------------------------------------------
data <- species %>%
  mutate(type = traits$targetHerb) %>%
  pivot_longer(names_to = "id", values_to = "n", cols = starts_with("X")) %>%
  group_by(id, type) %>%
  summarise(total = sum(n)) %>%
  filter(type == "yes") %>%
  select(-type) %>%
  ungroup() %>%
  mutate(targetHerbCov = round(total, 3), .keep = "unused")
sites <- left_join(sites, data, by = "id")

### d Arrhenatherum species' cover ratio -------------------------------------------------------------------------------------------
data <- species %>%
  mutate(type = traits$targetArrhenatherion) %>%
  pivot_longer(names_to = "id", values_to = "n", cols = starts_with("X")) %>%
  group_by(id, type) %>%
  summarise(total = sum(n)) %>%
  filter(type == "yes") %>%
  select(-type) %>%
  ungroup() %>%
  mutate(arrhCov = round(total, 3), .keep = "unused")
sites <- left_join(sites, data, by = "id")

### e Lean indicator's coverage -------------------------------------------------------------------------------------------
data <- species %>%
  mutate(type = traits$leanIndicator) %>%
  pivot_longer(names_to = "id", values_to = "n", cols = starts_with("X")) %>%
  group_by(id, type) %>%
  summarise(total = sum(n)) %>%
  filter(type == "yes") %>%
  select(-type) %>%
  ungroup() %>%
  mutate(leanCov = round(total, 3), .keep = "unused")
sites <- left_join(sites, data, by = "id")

### f Nitrogen indicator's coverage -------------------------------------------------------------------------------------------
data <- species %>%
  mutate(type = traits$nitrogenIndicator) %>%
  pivot_longer(names_to = "id", values_to = "n", cols = starts_with("X")) %>%
  group_by(id, type) %>%
  summarise(total = sum(n)) %>%
  filter(type == "yes") %>%
  select(-type) %>%
  ungroup() %>%
  mutate(nitrogenCov = round(total, 3), .keep = "unused")
sites <- left_join(sites, data, by = "id")

### g Ruderal indicator's coverage -------------------------------------------------------------------------------------------
data <- species %>%
  mutate(type = traits$ruderalIndicator) %>%
  pivot_longer(names_to = "id", values_to = "n", cols = starts_with("X")) %>%
  group_by(id, type) %>%
  summarise(total = sum(n)) %>%
  filter(type == "yes") %>%
  select(-type) %>%
  ungroup() %>%
  mutate(ruderalCov = round(total, 3), .keep = "unused")
sites <- left_join(sites, data, by = "id")

### h Table 33 species' coverage -------------------------------------------------------------------------------------------
data <- species %>%
  mutate(type = traits$table33) %>%
  pivot_longer(names_to = "id", values_to = "n", cols = starts_with("X")) %>%
  group_by(id, type) %>%
  summarise(total = sum(n)) %>%
  mutate(type = if_else(type == "4" | type =="3" | type == "2", "table33Cov", "other")) %>%
  group_by(id, type) %>%
  summarise(total = sum(total)) %>%
  spread(type, total) %>%
  select(id, table33Cov) %>%
  ungroup() %>%
  mutate(table33Cov = round(table33Cov, 1))
sites <- left_join(sites, data, by = "id")
rm(data)


### 3 Species richness #####################################################################################

specRich <- left_join(species, traits, by = "name") %>%
  select(rlg, rlb, target, targetHerb, targetArrhenatherion, ffh6510, ffh6210, nitrogenIndicator, leanIndicator, table33, table34, starts_with("X"))

### a total species richness -------------------------------------------------------------------------------------------
specRich_all <- specRich %>%
  pivot_longer(names_to = "id", values_to = "n", cols = starts_with("X")) %>%
  group_by(id) %>%
  mutate(n = if_else(n > 0, 1, 0)) %>%
  summarise(total = sum(n)) %>%
  group_by(id) %>%
  summarise(speciesRichness = sum(total)) %>%
  ungroup()

### b red list Germany (species richness) -------------------------------------------------------------------------------------------
specRich_rlg <- specRich %>%
  pivot_longer(names_to = "id", values_to = "n", cols = starts_with("X")) %>%
  group_by(id, rlg) %>%
  mutate(n = if_else(n > 0, 1, 0)) %>%
  summarise(total = sum(n)) %>%
  filter(rlg == "1" | rlg == "2" | rlg == "3" | rlg == "V") %>%
  group_by(id) %>%
  summarise(rlgRichness = sum(total)) %>%
  ungroup()

### c red list Bavaria (species richness) -------------------------------------------------------------------------------------------
specRich_rlb <- specRich %>%
  pivot_longer(names_to = "id", values_to = "n", cols = starts_with("X")) %>%
  group_by(id, rlb) %>%
  mutate(n = if_else(n > 0, 1, 0)) %>%
  summarise(total = sum(n)) %>%
  filter(rlb == "1" | rlb == "2" | rlb == "3" | rlb == "V") %>%
  group_by(id) %>%
  summarise(rlbRichness = sum(total)) %>%
  ungroup()

### d target species (species richness) -------------------------------------------------------------------------------------------
specRich_target <- specRich %>%
  pivot_longer(names_to = "id", values_to = "n", cols = starts_with("X")) %>%
  group_by(id, target) %>%
  mutate(n = if_else(n > 0, 1, 0)) %>%
  summarise(total = sum(n)) %>%
  filter(target == "yes") %>%
  group_by(id) %>%
  summarise(targetRichness = sum(total)) %>%
  ungroup()

### e target herb species (species richness) -------------------------------------------------------------------------------------------
specRich_targetHerb <- specRich %>%
  pivot_longer(names_to = "id", values_to = "n", cols = starts_with("X")) %>%
  group_by(id, targetHerb) %>%
  mutate(n = if_else(n > 0, 1, 0)) %>%
  summarise(total = sum(n)) %>%
  filter(targetHerb == "yes") %>%
  group_by(id) %>%
  summarise(targetHerbRichness = sum(total)) %>%
  ungroup()

### f Arrhenatherion species (species richness) -------------------------------------------------------------------------------------------
specRich_arrh <- specRich %>%
  pivot_longer(names_to = "id", values_to = "n", cols = starts_with("X")) %>%
  group_by(id, targetArrhenatherion) %>%
  mutate(n = if_else(n > 0, 1, 0)) %>%
  summarise(total = sum(n)) %>%
  filter(targetArrhenatherion == "yes") %>%
  group_by(id) %>%
  summarise(arrhRichness = sum(total)) %>%
  ungroup()

### g ffh6510 species (species richness) -------------------------------------------------------------------------------------------
specRich_ffh6510 <- specRich %>%
  pivot_longer(names_to = "id", values_to = "n", cols = starts_with("X")) %>%
  group_by(id, ffh6510) %>%
  mutate(n = if_else(n > 0, 1, 0)) %>%
  summarise(total = sum(n)) %>%
  filter(ffh6510 == "yes") %>%
  group_by(id) %>%
  summarise(ffh6510Richness = sum(total)) %>%
  ungroup()

### h ffh6210 species (species richness) -------------------------------------------------------------------------------------------
specRich_ffh6210 <- specRich %>%
  pivot_longer(names_to = "id", values_to = "n", cols = starts_with("X")) %>%
  group_by(id, ffh6210) %>%
  mutate(n = if_else(n > 0, 1, 0)) %>%
  summarise(total = sum(n)) %>%
  filter(ffh6210 == "yes") %>%
  group_by(id) %>%
  summarise(ffh6210Richness = sum(total)) %>%
  ungroup()

### i leanIndicator species (species richness) -------------------------------------------------------------------------------------------
specRich_leanIndicator <- specRich %>%
  pivot_longer(names_to = "id", values_to = "n", cols = starts_with("X")) %>%
  group_by(id, leanIndicator) %>%
  mutate(n = if_else(n > 0, 1, 0)) %>%
  summarise(total = sum(n)) %>%
  filter(leanIndicator == "yes") %>%
  group_by(id) %>%
  summarise(leanIndicatorRichness = sum(total)) %>%
  ungroup()

### k table33 (species richness) -------------------------------------------------------------------------------------------
specRich_table33_2 <- specRich %>%
  pivot_longer(names_to = "id", values_to = "n", cols = starts_with("X")) %>%
  group_by(id, table33) %>%
  mutate(n = if_else(n > 0, 1, 0)) %>%
  summarise(total = sum(n)) %>%
  filter(table33 == "2") %>%
  group_by(id) %>%
  summarise(table33_2Richness = sum(total)) %>%
  ungroup()
specRich_table33_3 <- specRich %>%
  pivot_longer(names_to = "id", values_to = "n", cols = starts_with("X")) %>%
  group_by(id, table33) %>%
  mutate(n = if_else(n > 0, 1, 0)) %>%
  summarise(total = sum(n)) %>%
  filter(table33 == "3") %>%
  group_by(id) %>%
  summarise(table33_3Richness = sum(total)) %>%
  ungroup()
specRich_table33_4 <- specRich %>%
  pivot_longer(names_to = "id", values_to = "n", cols = starts_with("X")) %>%
  group_by(id, table33) %>%
  mutate(n = if_else(n > 0, 1, 0)) %>%
  summarise(total = sum(n)) %>%
  filter(table33 == "4") %>%
  group_by(id) %>%
  summarise(table33_4Richness = sum(total)) %>%
  ungroup()

### l table34 (species richness) -------------------------------------------------------------------------------------------
specRich_table34_2 <- specRich %>%
  pivot_longer(names_to = "id", values_to = "n", cols = starts_with("X")) %>%
  group_by(id, table34) %>%
  mutate(n = if_else(n > 0, 1, 0)) %>%
  summarise(total = sum(n)) %>%
  filter(table34 == "2") %>%
  group_by(id) %>%
  summarise(table34_2Richness = sum(total)) %>%
  ungroup()
specRich_table34_3 <- specRich %>%
  pivot_longer(names_to = "id", values_to = "n", cols = starts_with("X")) %>%
  group_by(id, table34) %>%
  mutate(n = if_else(n > 0, 1, 0)) %>%
  summarise(total = sum(n)) %>%
  filter(table34 == "3") %>%
  group_by(id) %>%
  summarise(table34_3Richness = sum(total)) %>%
  ungroup()

### m implement in sites data set -------------------------------------------------------------------------------------------
sites <- left_join(sites, specRich_all, by = "id")
sites <- left_join(sites, specRich_rlg, by = "id")
sites <- left_join(sites, specRich_rlb, by = "id")
sites <- left_join(sites, specRich_target, by = "id")
sites <- left_join(sites, specRich_targetHerb, by = "id")
sites <- left_join(sites, specRich_arrh, by = "id")
sites <- left_join(sites, specRich_ffh6510, by = "id")
sites <- left_join(sites, specRich_ffh6210, by = "id")
sites <- left_join(sites, specRich_leanIndicator, by = "id")
sites <- left_join(sites, specRich_table33_2, by = "id")
sites <- left_join(sites, specRich_table33_3, by = "id")
sites <- left_join(sites, specRich_table33_4, by = "id")
sites <- left_join(sites, specRich_table34_2, by = "id")
sites <- left_join(sites, specRich_table34_3, by = "id")
sites <- sites %>%
  mutate(targetRichratio = targetRichness / speciesRichness) %>%
  mutate(targetRichratio = round(targetRichratio, 3))
rm(list=setdiff(ls(), c("sites", "species", "traits")))


### 4 Biotope types #####################################################################################
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
sites <- left_join(sites, biotopetypes, by = "id")  

#test <- left_join(sites, biotopetypes, by = "id") %>%
 # mutate_if(is.factor, as.character) %>%
  #select(id, biotopeType.x, biotopeType.y, arrhCov, arrhRichness, leanCov) %>%
  #mutate(similar = if_else(biotopeType.x == biotopeType.y, "yes", "no")) %>%
  #filter(similar == "no") %>%
  #select(-similar)
  #rm(test)
sites <- select(sites, -targetHerbCov, -arrhCov, -leanCov, -nitrogenCov, -table33Cov, -targetHerbRichness, -arrhRichness, -leanIndicatorRichness, -table33_2Richness, -table33_3Richness, table33_4Richness, table34_2Richness, table34_3Richness)
rm(biotopetypes)  


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
tbi2 <- tbi %>%
  filter(exposition == "north") %>%
  select(-exposition)
tbiAbuN <- tbi2[,colSums(tbi2) > 0]
tbi2 <- tbi %>%
  filter(exposition == "south") %>%
  select(-exposition)
tbiAbuS <- tbi2[,colSums(tbi2) > 0]
tbi2 <- tbi %>%
  filter(exposition == "east") %>%
  select(-exposition)
tbiAbuE <- tbi2[,colSums(tbi2) > 0]
tbi2 <- tbi %>%
  filter(exposition == "west") %>%
  select(-exposition)
tbiAbuW <- tbi2[,colSums(tbi2) > 0]
rm(tbi,tbi2)

tbiAbu <- rownames_to_column(tbiAbu, var = "id")
tbiPa <- rownames_to_column(tbiPa, var = "id")
tbiAbuN <- rownames_to_column(tbiAbuN, var = "id")
tbiAbuS <- rownames_to_column(tbiAbuS, var = "id")
tbiAbuE <- rownames_to_column(tbiAbuE, var = "id")
tbiAbuW <- rownames_to_column(tbiAbuW, var = "id")


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
rm(list=setdiff(ls(), c("sites", "species", "traits")))


### 7 Functional plant traits #####################################################################################

#### * make names equal #### 
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


data <- full_join(dataSLA, dataSM, by = "name")
data <- full_join(data, dataH, by = "name")
rm(dataSLA, dataSM, dataH)

#traits$name[which(!(traits$name %in% data$name))]
#data %>%
  #group_by(name) %>%
  #summarise(across(where(is.double), ~median(.x, na.rm = T))) %>%
  #filter(str_detect(name, "incana"))
##### Synonyme ###
data$name <- fct_recode(data$name, Centaurea_stoebe = "Centaurea_stoebe_s.lat.")
data$name <- fct_recode(data$name, Carex_praecox_ssp_praecox = "Carex_praecox")
data$name <- fct_recode(data$name, Cerastium_fontanum_ssp_vulgare = "Cerastium_fontanum")
data$name <- fct_recode(data$name, Cerastium_fontanum_ssp_vulgare = "Cerastium_fontanum_s._vulgare")
data$name <- fct_recode(data$name, Cota_tinctoria = "Anthemis_tinctoria")
data$name <- fct_recode(data$name, Cyanus_segetum = "Centaurea_cyanus")
data$name <- fct_recode(data$name, Euphorbia_verrucosa = "Euphorbia_brittingeri")
data$name <- fct_recode(data$name, Helictotrichon_pubescens = "Avenula_pubescens")
data$name <- fct_recode(data$name, Hypochaeris_radicata = "Hypochoeris_radicata")
data$name <- fct_recode(data$name, Jacobaea_vulgaris = "Senecio_vulgaris")
data$name <- fct_recode(data$name, Medicago_falcata = "Medicago_sativa_s._falcata")
data$name <- fct_recode(data$name, Ononis_spinosa_ssp_procurrens = "Ononis_repens")
data$name <- fct_recode(data$name, Persicaria_amphibia = "Polygonum_amphibium")
data$name <- fct_recode(data$name, Pilosella_caespitosa = "Hieracium_caespitosum")
data$name <- fct_recode(data$name, Pilosella_officinarum = "Hieracium_pilosella")
data$name <- fct_recode(data$name, Pilosella_piloselloides = "Hieracium_piloselloides")
data$name <- fct_recode(data$name, Plantago_major_ssp_intermedia = "Plantago_major_subsp._intermedia")
data$name <- fct_recode(data$name, Rubus_fruticosus_agg = "Rubus_fruticosus")
data$name <- fct_recode(data$name, Rubus_fruticosus_agg = "Rubus_fruticosus_ag._L.")
data$name <- fct_recode(data$name, Securigera_varia = "Coronilla_varia")
data$name <- fct_recode(data$name, "Silene_flos-cuculi" = "Lychnis_flos-cuculi")
data$name <- fct_recode(data$name, Taraxacum_campylodes = "Taraxacum_officinale")
data$name <- fct_recode(data$name, Taraxacum_campylodes = "Taraxacum_Sec._Ruderalia")
data$name <- fct_recode(data$name, Tripleurospermum_maritimum = "Matricaria_maritima")
data$name <- fct_recode(data$name, Vicia_sativa_ssp_nigra= "Vicia_sativa_s._nigra")
##### Take value of similar species ###
#data$name <- fct_recode(data$name, Ranunculus_serpens_ssp_nemorosus = "Ranunculus_polyanthemos")
#data$name <- fct_recode(data$name, Potentilla_grandiflora = "Potentilla_cinerea")

#### * add missing values ####
####  SLA (missing values) ###
####SLA von Cerabolini et al. (2010)
#####SLA von Pierce et al. (2007)
#####SLA von Pipenbaher et al. (2007)
#dataSLA <- dataSLA %>% add_row(name = "Euphorbia_brittingeri", value = 16.2) #Bei Pipenbaher als E. verrucosa gef?hrt
#dataSLA <- dataSLA %>% add_row(name = "Rhinanthus_aristatus", value = 15.9)
#data$name <- fct_recode(data$name, Carex_praecox_ssp_curvata = )

#### seed mass (missing values) ###
######Seedmass-Daten: Eigenannahme von Conradi & Kollmann (2016)
#####Seedmass-Daten: Hintze et al. (2013)
#dataSM <- dataSM %>% add_row(name = "Euphorbia_brittingeri", value = 2.635) # Bei Hintze als E. verrucosa gef?hrt
#dataSM <- dataSM %>% add_row(name = "Potentilla_cinerea", value = 0.131) # Bei Hintze als P. incana gef?hrt

#### canopy height (missing values) ###
#####Height von Pipenbaher et al. (2007)
#dataH <- dataH %>% add_row(name = "Rhinanthus_aristatus", value = 0.190)
#####Height von J?ger (2011)

#### * put values to traits data frame #### 
data <- data %>%
  group_by(name) %>%
  summarise(across(where(is.double), ~median(.x, na.rm = T)))
traits <- left_join(traits, data, by = "name")
rm(data)

### * Check completeness ####
test <- traits %>%
  select(name, t, n, f, sla, seedmass, height)
n_miss(test$sla);pct_complete(test$sla)
n_miss(test$seedmass);pct_complete(test$seedmass)
n_miss(test$height);pct_complete(test$height)
vis_miss(test, cluster = F, sort_miss = T)
gg_miss_var(test)
gg_miss_case(test, order_cases = F)
test2 <- test %>%
  select(-t, -n, -f) %>%
  filter(!complete.cases(.))
rm(test, test2)

### * Prepare data frames ####
traitsLHS <- traits %>%
  select(name, sla, seedmass, height) %>%
  drop_na()
traitsSLA <- traits %>%
  select(name, sla) %>%
  drop_na()  
traitsSM <- traits %>%
  select(name, seedmass) %>%
  drop_na()
traitsH <- traits %>%
  select(name, height) %>%
  drop_na()


### 8 CWM and FDis of functional plant traits #####################################################################################

### a All -------------------------------------------------------------------------------------------
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

### b Seed mass -------------------------------------------------------------------------------------------
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

### c Canopy height -------------------------------------------------------------------------------------------
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

traits <- select(traits, -targetArrhenatherion, -table30, -table33, -table34, -nitrogenIndicator, -nitrogenIndicator2, leanIndicator)



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# C Save processed data ##############################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


setwd("Z:/Documents/0_Donaudeiche/3_Aufnahmen_und_Ergebnisse/2022_Danube_old_dikes/data/processed")
write_csv2(sites, "data_processed_sites.csv")
write_csv2(species, "data_processed_species.csv")
write_csv2(traits, "data_processed_traits.csv")
write_csv2(tbiAbu, "data_processed_tbiAbu.csv")
write_csv2(tbiPa, "data_processed_tbiPa.csv")
