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
library(vegan) #metaMDS()
library(FD) #dbFD()
library(adespatial)
#remotes::install_github("larsito/tempo")
library(tempo) #calc_sync()


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
  select(-position, -starts_with("surveyDate_"), -starts_with("topsoilDepth_"), -cnLevel, -ends_with("Class"), -starts_with("sceleton")) %>%
  relocate(plot, .after = id) %>%
  relocate(c("locationAbb", "locationYear"), .after = location) %>%
  relocate(c("surveyYear", "surveyYearF", "surveyYearFminus"), .after = riverkm) %>%
  relocate(c("constructionYearF", "constructionYearFplus"), .after = constructionYear)


### 2 Species #####################################################################################

species <- data.table::fread("20211102_data_raw_species.csv", 
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
  select(name, all_of(sites$id)) %>%
  mutate(total = sum(c_across(starts_with("X")), na.rm = T),
         presence = if_else(total > 0, 1, 0),
         name = factor(name)) %>%
  filter(presence == 1) %>%
  ungroup() %>%
  select(name, sort(tidyselect::peek_vars()), -total, -presence) %>%
  mutate(across(where(is.numeric), ~replace(., is.na(.), 0)))

### Create list with species names and their frequency ###
specieslist <- species %>%
  mutate_if(is.numeric, ~1 * (. != 0)) %>%
  mutate(sum = rowSums(across(where(is.numeric)), na.rm = T), .keep = "unused") %>%
  group_by(name) %>%
  summarise(sum = sum(sum))
#write_csv(specieslist, "specieslist_2022xxxx.csv")


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
  select(-ssp) %>%
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
  summarise(across(where(is.double), ~ sum(.x, na.rm = T))) %>%
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
  select(name, starts_with("X01")) %>%
  filter(if_any(starts_with("X"), ~ . > 0)) %>%
  print(n = 70)

### Check missing data ###
miss_var_summary(sites, order = T)
vis_miss(sites, cluster = F, sort_miss = T)
vis_miss(traits, cluster = F, sort_miss = T)

rm(list = ls(pattern = "[^species|traits|sites]"))



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Create variables ##############################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


## 1 Create simple variables #####################################################################################

traits <- traits %>%
  mutate(leanIndicator = if_else(
    !(is.na(table30)) | !(is.na(table33)) | !(is.na(table34)), "yes", "no"
    ),
    target = if_else(
      targetHerb == "yes" | targetGrass == "yes", "yes", "no"
      ),
    ruderal = if_else(
      sociology >= 3300 & sociology < 3700, "yes", "no"
      ),
    targetEllenberg = if_else(
      sociology >= 5300 & sociology < 5400, "dry_grassland", if_else(
        sociology >= 5400 & sociology < 6000, "hay_meadow", if_else(
          sociology >= 5100 & sociology < 5200, "nardus_grassland", if_else(
            sociology >= 5200 & sociology < 5300, "sand_grasland", "no"
            )))))

sites <- sites %>%
  mutate(conf.low = c(1:length(id)),
         conf.high = c(1:length(id)),
         fmMass = fmDepth * fmDbd * 10,
         fmMass = round(fmMass, 3),
         NtotalConc = fmMass * NtotalPerc / 100,
         plotAge = surveyYear - constructionYear,
         ageCategory = if_else(surveyYear > 2016, "young", "old")) %>%
  select(-fmDepth, -fmMass)


## 2 Coverages #####################################################################################

cover <- left_join(species, traits, by = "name") %>%
  select(name, family, target, targetHerb, targetArrhenatherion, leanIndicator, nitrogenIndicator, ruderalIndicator, table33, starts_with("X")) %>%
  pivot_longer(names_to = "id", values_to = "n", cols = starts_with("X")) %>%
  group_by(id)

### * graminoid, herb, and total coverage) ####
cover_total_and_graminoid <- cover %>%
  group_by(id, family) %>%
  summarise(total = sum(n, na.rm = T), .groups = "keep") %>%
  mutate(type = if_else(family == "Poaceae" | family == "Cyperaceae" | family == "Juncaceae", "graminoidCov", "herbCov")) %>%
  group_by(id, type) %>%
  summarise(total = sum(total, na.rm = T), .groups = "keep") %>%
  spread(type, total) %>%
  mutate(accumulatedCov = graminoidCov + herbCov) %>%
  ungroup()

### * Target species' coverage ####
cover_target <- cover %>%
  filter(target == "yes") %>%
  summarise(targetCov = sum(n, na.rm = T)) %>%
  mutate(targetCov = round(targetCov, 1)) %>%
  ungroup()

### * Target herb species' coverage ####
cover_targetHerb <- cover %>%
  filter(targetHerb == "yes") %>%
  summarise(targetHerbCov = sum(n, na.rm = T)) %>%
  mutate(targetHerbCov = round(targetHerbCov, 1)) %>%
  ungroup()

### * Arrhenatherum species' cover ratio ####
cover_targetArrhenatherion <- cover %>%
  filter(targetArrhenatherion == "yes") %>%
  summarise(arrhCov = sum(n, na.rm = T)) %>%
  mutate(arrhCov = round(arrhCov, 1)) %>%
  ungroup()

### * Lean indicator's coverage ####
cover_leanIndicator <- cover %>%
  filter(leanIndicator == "yes") %>%
  summarise(leanCov = sum(n, na.rm = T)) %>%
  mutate(leanCov = round(leanCov, 1)) %>%
  ungroup()

### * Nitrogen indicator's coverage ####
cover_nitrogenIndicator <- cover %>%
  filter(nitrogenIndicator == "yes") %>%
  summarise(nitrogenCov = sum(n, na.rm = T)) %>%
  mutate(nitrogenCov = round(nitrogenCov, 1)) %>%
  ungroup()

### * Ruderal indicator's coverage ####
cover_ruderalIndicator <- cover %>%
  filter(ruderalIndicator == "yes") %>%
  summarise(ruderalCov = sum(n, na.rm = T)) %>%
  mutate(ruderalCov = round(ruderalCov, 1)) %>%
  ungroup()
  
### * Table 33 species' coverage ####
cover_table33 <- cover %>%
  mutate(table33 = if_else(table33 == "4" | table33 == "3" | table33 == "2", "table33Cov", "other")) %>%
  filter(table33 == "table33Cov") %>%
  summarise(table33Cov = sum(n, na.rm = T)) %>%
  mutate(table33Cov = round(table33Cov, 1)) %>%
  ungroup()

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
         graminoidCovratio = graminoidCov / accumulatedCov)

rm(list = ls(pattern = "[^species|traits|sites]"))


## 3 Alpha diversity #####################################################################################

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
  mutate(targetRichratio = targetRichness / speciesRichness)

### b Species eveness and shannon ----------------------------------------------------------------------

data <- species  %>%
  mutate(across(where(is.numeric), ~replace(., is.na(.), 0))) %>%
  pivot_longer(-name, names_to = "id", values_to = "value") %>%
  pivot_wider(names_from = "name", values_from = "value") %>%
  column_to_rownames("id") %>%
  diversity(index = "shannon") %>%
  as_tibble(rownames = NA) %>%
  rownames_to_column(var = "id") %>%
  mutate(id = factor(id)) %>%
  rename(shannon = value)
sites <- sites %>%
  left_join(data, by = "id") %>%
  mutate(eveness = shannon / log(speciesRichness))

rm(list = ls(pattern = "[^species|traits|sites]"))


## 4 Biotope types #####################################################################################

### a Calculate types -------------------------------------------------------------------------------------------

biotopetypes <- sites %>%
  select(id, table33_2Richness, table33_3Richness, table33_4Richness, table33Cov, table34_2Richness, table34_3Richness, targetRichness, targetHerbRichness, arrhRichness, targetCov, leanCov, arrhCov, targetHerbCov, nitrogenCov) %>%
  mutate(table33Rich_proof = if_else(
    table33_2Richness >= 2 | table33_2Richness + table33_3Richness >= 3 | table33_2Richness + table33_3Richness + table33_4Richness >= 4, "yes", "no"
    ),
    table33Cov_proof = if_else(
      table33Cov >= 25, "yes", "no"
      ),
    table33_proof = if_else(
      table33Rich_proof == "yes" & table33Cov_proof == "yes", "yes", "no"
      ),
    table34Rich_proof = if_else(
      table34_2Richness >= 2 | table34_2Richness + table34_3Richness >= 3, "yes", "no"
      ),
    G312_GT6210_type = if_else(
      table33_proof == "yes" & table34Rich_proof == "yes", "yes", "no"
      ),
    GE_proof = if_else(
      targetHerbCov >= 12.5 & targetRichness >= 20 & nitrogenCov < 25, "yes", "no"
      ),
    G214_GE6510_type = if_else(
      GE_proof == "yes" & arrhRichness >= 1 & arrhCov > 0.5 & leanCov >= 25 & targetHerbCov >= 12.5, "yes", "no"
      ),
    G212_LR6510_type = if_else(
      GE_proof == "yes" & arrhRichness >= 1 & arrhCov > 0.5 & leanCov < 25, "yes", "no"
      ),
    G214_GE00BK_type = if_else(
      GE_proof == "yes" & arrhRichness == 0 & leanCov >= 25 & targetHerbCov >= 12.5, "yes", if_else(
        GE_proof == "yes" & arrhCov >= 0.5 & leanCov >= 25 & targetHerbCov >= 12.5, "yes", "no"
        )),
    G213_GE00BK_type = if_else(
      GE_proof == "yes" & arrhRichness == 0 & leanCov >= 25 & targetHerbCov < 12.5, "yes", if_else(
        GE_proof == "yes" & arrhCov >= 0.5 & leanCov >= 25 & targetHerbCov < 12.5, "yes", "no"
      )),
    G213_type = if_else(
      leanCov >= 25, "yes", "no"
      ),
    G212_type = if_else(
      leanCov >= 1 & leanCov < 25 & targetHerbCov >= 12.5 & targetHerbRichness >= 10, "yes", "no"
      ),
    G211_type = if_else(
      leanCov >= 1 & leanCov < 25 & targetHerbCov >= 1 & targetHerbRichness >= 5, "yes", "no"
      ),
    biotopeType = if_else(
      G312_GT6210_type == "yes", "G312-GT6210", if_else(
        G214_GE6510_type == "yes", "G214-GE6510", if_else(
          G212_LR6510_type == "yes", "G212-LR6510", if_else(
            G214_GE00BK_type == "yes", "G214-GE00BK", if_else(
              G213_GE00BK_type == "yes", "G213-GE00BK", if_else(
                G213_type == "yes", "G213", if_else(
                  G212_type == "yes", "G212", if_else(
                    G211_type == "yes", "G211", "other"
                    )))))))),
    biotopeType = as_factor(biotopeType)) %>%
  select(id, biotopeType, -ends_with("proof"), -starts_with("table")) %>%
  mutate(ffh6510 = str_match(biotopeType, "6510"),
         ffh6210 = str_match(biotopeType, "6210"),
         baykompv = as_factor(str_sub(biotopeType, start = 1, end = 4))) %>%
  unite(ffh, ffh6510, ffh6210, sep = "") %>%
  mutate(ffh = str_replace(ffh, "NA", ""),
         ffh = as_factor(str_replace(ffh, "NA", "non-FFH")),
         biotopeType = as_factor(biotopeType),
         biotopePoints = if_else(
           biotopeType == "G312-GT6210", 13, if_else(
             biotopeType == "G214-GE6510", 12, if_else(
               biotopeType == "G212-LR6510", 9, if_else(
                 biotopeType == "G214-GE00BK", 12, if_else(
                   biotopeType == "G213-GE00BK", 9, if_else(
                     biotopeType == "G213", 8, if_else(
                       biotopeType == "G212", 8, if_else(
                         biotopeType == "G211", 6, 0
                         )))))))),
         min8 = as_factor(if_else(biotopePoints >= 8, "yes", "no")),
         min9 = as_factor(if_else(biotopePoints >= 9, "yes", "no")))
sites <- left_join(sites, biotopetypes, by = "id") %>%  
  select(-targetHerbCov, -arrhCov, -leanCov, -nitrogenCov, -table33Cov, -targetHerbRichness, -arrhRichness, -leanIndicatorRichness, -ffh6510Richness, -ffh6210Richness, -table33_2Richness, -table33_3Richness, -table33_4Richness, -table34_2Richness, -table34_3Richness)
traits <- traits %>%
  select(-targetArrhenatherion, -table30, -table33, -table34, -nitrogenIndicator, -nitrogenIndicator2, leanIndicator)

### b Calculate constance -------------------------------------------------------------------------------------------

data <- sites %>%
  select(id, plot, surveyYear, ffh) %>%
  group_by(plot) %>%
  mutate(count = n()) %>%
  filter(count == max(count)) %>%
  pivot_wider(id_cols = -id, names_from = "surveyYear", values_from = "ffh") %>%  #group_by(plot) %>%
  rename(x17 = "2017", x18 = "2018", x19 = "2019") %>%
  mutate(changeType = ifelse((x17 == "non-FFH" & x18 != "non-FFH" & x19 != "non-FFH"), "better", ifelse(
    (x17 != "non-FFH" & x18 == "non-FFH" & x19 != "non-FFH") | (x17 == "non-FFH" & x18 != "non-FFH" & x19 == "non-FFH") | (x17 == "non-FFH" & x18 == "non-FFH" & x19 != "non-FFH") | (x17 != "non-FFH" & x18 != "non-FFH" & x19 == "non-FFH"), "change", ifelse((x17 == "non-FFH" & x18 == "non-FFH" & x19 == "non-FFH"), "non-FFH", ifelse(
        x17 == "6510" & x18 == "6510" & x19 == "6510", "FFH6510", ifelse(
          x17 == "6210" & x18 == "6210" & x19 == "6210", "FFH6210", ifelse(
            x17 != "non-FFH" & x18 != "non-FFH" & x19 != "non-FFH", "any-FFH", "worse"
            ))))))) %>%
  select(plot, changeType)
sites <- left_join(sites, data, by = "plot")

rm(list = ls(pattern = "[^species|traits|sites]"))


## 5 Beta diversity #####################################################################################

### * Prepare data ####
data_sites <- sites %>%
  filter(accumulatedCov > 0)
data_species <- species %>%
  pivot_longer(-name, "id", "value") %>%
  pivot_wider(id, name) %>%
  arrange(id) %>%
  semi_join(data_sites, by = "id") %>%
  mutate(across(where(is.numeric), ~replace(., is.na(.), 0)))

### a NMDS -------------------------------------------------------------------------------------------

data_species_nmds <- data_species %>%
  column_to_rownames(var = "id")
### Calculate NMDS ###
set.seed(1)
(nmds <- metaMDS(data_species_nmds, dist = "bray", binary = F,
                 try = 50, previous.best = T, na.rm = T))
### Add to sites ###
data <- nmds %>%
  scores() %>%
  as.data.frame() %>%
  rownames_to_column(var = "id") %>%
  as_tibble() %>%
  select(id, NMDS1, NMDS2)
sites <- left_join(sites, data, by = "id")

### b PERMDISP -------------------------------------------------------------------------------------------

### Presence-Absence data ###
data <- betadisper(d = vegdist(data_species_nmds, method = "bray", binary = T),
                   group = data_sites$plot)
data <- enframe(data$distances) %>%
  rename(id = name, 
         permdispPresabs = value)
sites <- left_join(sites, data, by = "id")

### Abundance data ###
### Presence-Absence data ###
data <- betadisper(d = vegdist(data_species_nmds, method = "bray", binary = F),
                   group = data_sites$plot)
data <- enframe(data$distances) %>%
  rename(id = name, 
         permdispAbu = value)
sites <- left_join(sites, data, by = "id")

rm(list = ls(pattern = "[^species|traits|sites|data_sites|data_species]"))

### c dbMEM -------------------------------------------------------------------------------------------

source('https://raw.githubusercontent.com/zdealveindy/anadat-r/master/scripts/NumEcolR2/quickMEM.R')
### * 2017 ####
data_sites_dbMEM <- data_sites %>%
  filter(surveyYear == 2017) %>%
  select(id, longitude, latitude)
data_species_dbMEM <- data_species %>%
  semi_join(data_sites_dbMEM, by = "id") %>%
  column_to_rownames(var = "id") %>%
  decostand("hellinger")
data_sites_dbMEM <- data_sites_dbMEM %>%
  column_to_rownames("id")
m <- quickMEM(data_species_dbMEM, data_sites_dbMEM, 
              alpha = 0.05, 
              detrend = F,
              method = "fwd",
              rangexy = T,
              perm.max = 999) #R2adj of minimum (final) model = 0.096
m$RDA_test # p = 0.001
m$RDA_axes_test #1 sig axes
m$RDA # PC1 = 0.093
dbMEMred <- dbMEMred %>%
  rownames_to_column(var = "id") %>%
  select(id, MEM1) %>%
  rename(MEM1_2017 = MEM1)
sites <- left_join(sites, dbMEMred, by = "id")

### * 2018 ####
data_sites_dbMEM <- data_sites %>%
  filter(surveyYear == 2018) %>%
  select(id, longitude, latitude)
data_species_dbMEM <- data_species %>%
  semi_join(data_sites_dbMEM, by = "id") %>%
  column_to_rownames(var = "id") %>%
  decostand("hellinger")
data_sites_dbMEM <- data_sites_dbMEM %>%
  column_to_rownames("id")
m <- quickMEM(data_species_dbMEM, data_sites_dbMEM, 
              alpha = 0.05, 
              detrend = F,
              method = "fwd",
              rangexy = T,
              perm.max = 999) #R2adj of minimum (final) model = 0.086
m$RDA_test # p = 0.001
m$RDA_axes_test # 2 sig. axes
m$RDA # PC1 = 0.090
dbMEMred <- dbMEMred %>%
  rownames_to_column(var = "id") %>%
  select(id, MEM1, MEM2) %>%
  rename(MEM1_2018 = MEM1, MEM2_2018 = MEM2)
sites <- left_join(sites, dbMEMred, by = "id")

### * 2019 ####
data_sites_dbMEM <- data_sites %>%
  filter(surveyYear == 2019) %>%
  select(id, longitude, latitude)
data_species_dbMEM <- data_species %>%
  semi_join(data_sites_dbMEM, by = "id") %>%
  column_to_rownames(var = "id") %>%
  decostand("hellinger")
data_sites_dbMEM <- data_sites_dbMEM %>%
  column_to_rownames("id")
m <- quickMEM(data_species_dbMEM, data_sites_dbMEM, 
              alpha = 0.05, 
              detrend = F,
              method = "fwd",
              rangexy = T,
              perm.max = 999) #R2adj of minimum (final) model = 0.093 
m$RDA_test # p = 0.001
m$RDA_axes_test # 2sig axes
m$RDA # PC1 = 0.065, PC2 = 0.060
dbMEMred <- dbMEMred %>%
  rownames_to_column(var = "id") %>%
  select(id, MEM1, MEM2) %>%
  rename(MEM1_2019 = MEM1, MEM2_2019 = MEM2)
sites <- left_join(sites, dbMEMred, by = "id")

### * 2021 ####
data_sites_dbMEM <- data_sites %>%
  filter(surveyYear == 2021) %>%
  select(id, longitude, latitude)
data_species_dbMEM <- data_species %>%
  semi_join(data_sites_dbMEM, by = "id") %>%
  column_to_rownames(var = "id") %>%
  decostand("hellinger")
data_sites_dbMEM <- data_sites_dbMEM %>%
  column_to_rownames("id")
m <- quickMEM(data_species_dbMEM, data_sites_dbMEM, 
              alpha = 0.05, 
              detrend = F,
              method = "fwd",
              rangexy = T,
              perm.max = 999) #R2adj of minimum (final) model = 0.064 
m$RDA_test # p = 0.002
m$RDA_axes_test # 2 sig axes
m$RDA # PC1 = 0.096, PC2 = 0.085
dbMEMred <- dbMEMred %>%
  rownames_to_column(var = "id") %>%
  select(id, MEM1, MEM2) %>%
  rename(MEM1_2021 = MEM1, MEM2_2021 = MEM2)
sites <- left_join(sites, dbMEMred, by = "id")

#### d LCBD (Local Contributions to Beta Diversity) -------------------------------------------------------------------------------------------

#### * 2017 ####
data_species_lcbd <- data_species %>%
  filter(str_detect(id, "2017")) %>%
  column_to_rownames(var = "id")
data <- beta.div(data_species_lcbd, method = "hellinger", nperm = 9999)
data_2017 <- data$LCBD
row.names(data_species_lcbd[which(p.adjust(data$p.LCBD, "holm") <= 0.05),]) #X62_m_2017

#### * 2018 ####
data_species_lcbd <- data_species %>%
  filter(str_detect(id, "2018")) %>%
  column_to_rownames(var = "id")
data <- beta.div(data_species_lcbd, method = "hellinger", nperm = 9999)
data_2018 <- data$LCBD
row.names(data_species_lcbd[which(p.adjust(data$p.LCBD, "holm") <= 0.05),]) # none

#### * 2019 ####
data_species_lcbd <- data_species %>%
  filter(str_detect(id, "2019")) %>%
  column_to_rownames(var = "id")
data <- beta.div(data_species_lcbd, method = "hellinger", nperm = 9999)
data_2019 <- data$LCBD
row.names(data_species_lcbd[which(p.adjust(data$p.LCBD, "holm") <= 0.05),]) # none

#### * 2021 ####
data_species_lcbd <- data_species %>%
  filter(str_detect(id, "2021")) %>%
  column_to_rownames(var = "id")
data <- beta.div(data_species_lcbd, method = "hellinger", nperm = 9999)
data_2021 <- data$LCBD
row.names(data_species_lcbd[which(p.adjust(data$p.LCBD, "holm") <= 0.05),]) # none

#### * combine datasets ####
data <- c(data_2017, data_2018, data_2019, data_2021)
data <- data_sites %>%
  mutate(lcbd = data) %>%
  select(id, lcbd)
sites <- sites %>%
  left_join(data, by = "id")

rm(list = ls(pattern = "[^species|traits|sites]"))

### e Synchrony -------------------------------------------------------------------------------------------

data_sites <- sites %>%
  select(id, plot, vegetationCov) %>%
  filter(vegetationCov > 0) %>%
  add_count(plot) %>%
  filter(n == max(n))
data_species <- species %>%
  select(where(~!all(is.na(.x)))) %>%
  mutate(across(where(is.numeric), ~replace(., is.na(.), 0))) %>%
  pivot_longer(-name, names_to = "id", values_to = "value") %>%
  pivot_wider(id, name) %>%
  mutate(year = str_match(id, "\\d{4}"),
         plot = str_match(id, "\\d{2}"),
         year = factor(year),
         plot = factor(plot)) %>%
  arrange(id) %>%
  semi_join(data_sites, by = "id") %>%
  column_to_rownames(var = "id") %>%
  select(plot, year, tidyselect::peek_vars())
data <- data_species %>%
  split(data_species$plot, drop = T) %>%
  map(~ (.x %>% select(-plot, -year))) %>%
  map(~ (.x %>% select(where(~!all(is.na(.x))))))
sync_indices <- map(data, calc_sync)
data <- do.call("rbind", sync_indices) %>% #map() does not work because 'plot' is missing
  rownames_to_column("plot") %>%
  as_tibble() %>%
  mutate(plot = str_extract(plot, "\\d\\d")) %>%
  select(plot, syn_total, syn_trend, syn_detrend) #log_varrrat_t3 does not allow missing years
sites <- left_join(sites, data, by = "plot")

rm(list = ls(pattern = "[^species|traits|sites]"))


## 6 Environmental variables #####################################################################################

### a Soil PCA  -------------------------------------------------------------------------------------------

### Prepare data ###
data <- sites %>%
  select(id, plot, calciumcarbonatPerc, humusPerc, NtotalPerc, cnRatio, pH, sandPerc, siltPerc, clayPerc, phosphorus, potassium, magnesium, topsoilDepth, NtotalConc) %>%
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
### create summary table ###
pcaSoil <- values$species[ ,1:3] %>%
  as_tibble() %>%
  bind_cols(c("calciumcarbonatPerc", "humusPerc", "NtotalPerc", "cnRatio", "pH", "sandPerc", "siltPerc", "clayPerc", "phosphorous", "potassium", "magnesium", "topsoilDepth", "NtotalConc")) %>%
  rename(variables = "...4") %>%
  bind_rows(eigenvals)
data <- as_tibble(values$sites[,1:3]) %>%
  rename(PC1soil = PC1,
         PC2soil = PC2,
         PC3soil = PC3)
### Add to sites ###
data <- sites %>%
  group_by(plot) %>%
  slice(1) %>%
  ungroup() %>%
  select(plot) %>%
  bind_cols(data)
sites <- left_join(sites, data, by = "plot") %>%
  select(-calciumcarbonatPerc, -humusPerc, -cnRatio, -pH, -sandPerc, -siltPerc, -clayPerc, -phosphorus, -potassium, -magnesium, -NtotalConc, -topsoilDepth,
         -HCl, -C550Perc, -CorgPerc, -humusLevel, -ufcPerc, -ufc, -fmDbd)

rm(list = ls(pattern = "[^species|traits|sites|species2017|species2018|species2019|species2021|pcaSoil]"))

### b Climate PCA  -------------------------------------------------------------------------------------------

### * Temperature ####
data <- read_csv(here("data/raw/temperature/data/data_OBS_DEU_P1M_T2M.csv"), col_names = T, na = c("", "NA", "na"), col_types = 
                   cols(
                     .default = "?"
                   )) %>%
  rename(date = Zeitstempel, value = Wert, site = SDO_ID) %>%
  select(site, date, value) %>%
  filter(date >= "2002-03-01") %>%
  mutate(site = factor(site),
         season = floor_date(date, "season"),
         year = year(season), 
         season = month(season),
         season = factor(season),
         season = fct_recode(season, "spring" = "3", "summer" = "6", "autumn" = "9", "winter" = "12")) %>%
  group_by(year) %>%
  mutate(yearMean = mean(value),
         yearMean = round(yearMean, digits = 1),
         currentYear = if_else(season == "spring", 0, 1),
         currentYear = year + currentYear) %>% #current year of 2021 is e.g. from summer 2020 to spring 2021 = climate for surveyYear
  group_by(currentYear) %>%
  mutate(currentMean = mean(value),
         currentMean = round(currentMean, digits = 1),
         currentMean = if_else(season == "spring", currentMean, NA_real_)) %>%
  group_by(year, season, yearMean, currentYear, currentMean) %>%
  summarise(seasonMean = round(mean(value), digits = 1), .groups = "keep") %>%
  pivot_wider(id_cols = c(year, yearMean, currentMean), names_from = season, values_from = seasonMean) %>%
  group_by(year) %>%
  summarise(across(where(is.numeric), ~ max(., na.rm = T)), .groups = "keep") %>% #warnings because of lates year (summer, autumn, winter), can be ignored
  mutate(year = factor(year))

sites <- sites %>%
  left_join(data %>% select(-currentMean), by = c("constructionYearF" = "year")) %>%
  rename("tempSpring_constructionYear" = "spring", 
         "tempSummer_constructionYear" = "summer", 
         "tempAutumn_constructionYear" = "autumn", 
         "tempWinter_constructionYear" = "winter", 
         "tempMean_constructionYear" = "yearMean") %>%
  left_join(data %>% select(-currentMean), by = c("constructionYearFplus" = "year")) %>%
  rename("tempSpring_constructionYearPlus" = "spring", 
         "tempSummer_constructionYearPlus" = "summer", 
         "tempAutumn_constructionYearPlus" = "autumn", 
         "tempWinter_constructionYearPlus" = "winter", 
         "tempMean_constructionYearPlus" = "yearMean") %>%
  left_join(data %>% select(year, currentMean, spring), by = c("surveyYearF" = "year")) %>%
  rename("tempMean_surveyYear" = "currentMean", 
         "tempSpring_surveyYear" = "spring") %>%
  left_join(data %>% select(year, summer, autumn, winter), by = c("surveyYearFminus" = "year")) %>%
  rename("tempSummer_surveyYear" = "summer", 
         "tempAutumn_surveyYear" = "autumn", 
         "tempWinter_surveyYear" = "winter")
rm(list = ls(pattern = "[^species|traits|sites|species2017|species2018|species2019|species2021|pcaSoil|pcaSurveyYear|pcaConstructionYear]"))

### * Precipitation ####
data <- read_csv(here("data/raw/precipitation/data/data_OBS_DEU_P1M_RR.csv"), col_names = T, na = c("", "NA", "na"), col_types = 
                   cols(
                     .default = "?"
                   )) %>%
  rename(date = Zeitstempel, value = Wert, site = SDO_ID) %>%
  select(site, date, value) %>%
  filter(date >= "2002-03-01") %>%
  mutate(site = factor(site),
         season = floor_date(date, "season"),
         year = year(season), 
         season = month(season),
         season = factor(season),
         season = fct_recode(season, "spring" = "3", "summer" = "6", "autumn" = "9", "winter" = "12")) %>%
  group_by(site, year) %>%
  mutate(yearSum = round(sum(value), digits = 0),
         currentYear = if_else(season == "spring", 0, 1),
         currentYear = year + currentYear) %>% #current year of 2021 is e.g. from summer 2020 to spring 2021 = climate for surveyYear
  group_by(site, currentYear) %>%
  mutate(currentSum = sum(value),
         currentSum = round(currentSum, 0)) %>%
  group_by(site, season, year, currentYear, yearSum, currentSum) %>%
  summarise(seasonSum = round(sum(value), digits = 0),
            .groups = "keep") %>%
  group_by(year, season, currentYear) %>%
  summarise(seasonMean = round(mean(seasonSum), digits = 0),
            yearMean = round(mean(yearSum), digits = 0),
            currentMean = round(mean(currentSum), digits = 0),
            .groups = "keep") %>%
  mutate(currentMean = if_else(season == "spring", currentMean, NA_real_)) %>%
  pivot_wider(id_cols = c(year, yearMean, currentMean), names_from = season, values_from = seasonMean) %>%
  group_by(year) %>%
  summarise(across(where(is.numeric), ~ max(., na.rm = T)),
            .groups = "keep") %>% #warnings because of lates year (summer, autumn, winter), can be ignored
  mutate(year = factor(year))
sites <- sites %>%
  left_join(data %>% select(-currentMean), by = c("constructionYearF" = "year")) %>%
  rename("precSpring_constructionYear" = "spring", 
         "precSummer_constructionYear" = "summer", 
         "precAutumn_constructionYear" = "autumn", 
         "precWinter_constructionYear" = "winter", 
         "precMean_constructionYear" = "yearMean") %>%
  left_join(data %>% select(-currentMean), by = c("constructionYearFplus" = "year")) %>%
  rename("precSpring_constructionYearPlus" = "spring", 
         "precSummer_constructionYearPlus" = "summer", 
         "precAutumn_constructionYearPlus" = "autumn", 
         "precWinter_constructionYearPlus" = "winter", 
         "precMean_constructionYearPlus" = "yearMean") %>%
  left_join(data %>% select(year, currentMean, spring), by = c("surveyYearF" = "year")) %>%
  rename("precMean_surveyYear" = "currentMean", 
         "precSpring_surveyYear" = "spring") %>%
  left_join(data %>% select(year, summer, autumn, winter), by = c("surveyYearFminus" = "year")) %>%
  rename("precSummer_surveyYear" = "summer", 
         "precAutumn_surveyYear" = "autumn", 
         "precWinter_surveyYear" = "winter")
rm(list = ls(pattern = "[^species|traits|sites|species2017|species2018|species2019|species2021|pcaSoil|pcaSurveyYear|pcaConstructionYear]"))

### * Calculation surveyYear ####
### Prepare data ###
data <- sites %>%
  arrange(id) %>%
  select(tempMean_surveyYear, tempSpring_surveyYear, tempSummer_surveyYear, tempAutumn_surveyYear, tempWinter_surveyYear,
         precMean_surveyYear, precSpring_surveyYear, precSummer_surveyYear, precAutumn_surveyYear, precWinter_surveyYear)
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
### create summary table ###
pcaSurveyYear <- values$species[ ,1:3] %>%
  as_tibble() %>%
  bind_cols(c("tempMean_surveyYear", "tempSpring_surveyYear", "tempSummer_surveyYear", "tempAutumn_surveyYear", "tempWinter_surveyYear",
              "precMean_surveyYear", "precSpring_surveyYear", "precSummer_surveyYear", "precAutumn_surveyYear", "precWinter_surveyYear")) %>%
  rename(variables = "...4") %>%
  bind_rows(eigenvals)
data <- as_tibble(values$sites[,1:3]) %>%
  rename(PC1surveyYear = PC1,
         PC2surveyYear = PC2,
         PC3surveyYear = PC3)
### Add to sites ###
sites <- sites %>%
  bind_cols(data) %>%
  select(-tempSpring_surveyYear, -tempSummer_surveyYear, -tempAutumn_surveyYear, -tempWinter_surveyYear,
         -precSpring_surveyYear, -precSummer_surveyYear, -precAutumn_surveyYear, -precWinter_surveyYear)

### * Calculation constructionYear ####
### Prepare data ###
data <- sites %>%
  arrange(id) %>%
  select(plot,
         tempMean_constructionYear, tempSpring_constructionYear, tempSummer_constructionYear, tempAutumn_constructionYear, tempWinter_constructionYear,
         tempMean_constructionYearPlus, tempSpring_constructionYearPlus, tempSummer_constructionYearPlus, tempAutumn_constructionYearPlus, tempWinter_constructionYearPlus,
         precMean_constructionYear, precSpring_constructionYear, precSummer_constructionYear, precAutumn_constructionYear, precWinter_constructionYear,
         precMean_constructionYearPlus, precSpring_constructionYearPlus, precSummer_constructionYearPlus, precAutumn_constructionYearPlus, precWinter_constructionYearPlus) %>%
  group_by(plot) %>%
  summarise(across(where(is.numeric), ~ median(.x, na.rm = T))) %>%
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
### create summary table ###
pcaConstuctionYear <- values$species[ ,1:3] %>%
  as_tibble() %>%
  bind_cols(c("tempMean_constructionYear", "tempSpring_constructionYear", "tempSummer_constructionYear", "tempAutumn_constructionYear", "tempWinter_constructionYear",
              "tempMean_constructionYearPlus", "tempSpring_constructionYearPlus", "tempSummer_constructionYearPlus", "tempAutumn_constructionYearPlus", "tempWinter_constructionYearPlus",
              "precMean_constructionYear", "precSpring_constructionYear", "precSummer_constructionYear", "precAutumn_constructionYear", "precWinter_constructionYear",
              "precMean_constructionYearPlus", "precSpring_constructionYearPlus", "precSummer_constructionYearPlus", "precAutumn_constructionYearPlus", "precWinter_constructionYearPlus")) %>%
  rename(variables = "...4") %>%
  bind_rows(eigenvals)
data <- as_tibble(values$sites[,1:3]) %>%
  rename(PC1constructionYear = PC1,
         PC2constructionYear = PC2,
         PC3constructionYear = PC3)
### Add to sites ###
data <- sites %>%
  group_by(plot) %>%
  slice(1) %>%
  ungroup() %>%
  select(plot) %>%
  bind_cols(data)
sites <- left_join(sites, data, by = "plot") %>%
  select(-tempSpring_constructionYear, -tempSummer_constructionYear, -tempAutumn_constructionYear, -tempWinter_constructionYear,
         -tempSpring_constructionYearPlus, -tempSummer_constructionYearPlus, -tempAutumn_constructionYearPlus, -tempWinter_constructionYearPlus,
         -precSpring_constructionYear, -precSummer_constructionYear, -precAutumn_constructionYear, -precWinter_constructionYear,
         -precSpring_constructionYearPlus, -precSummer_constructionYearPlus, -precAutumn_constructionYearPlus, -precWinter_constructionYearPlus)

rm(list = ls(pattern = "[^species|traits|sites|pcaSoil|pcaSurveyYear|pcaConstructionYear]"))


## 7 TBI: Temporal Beta diversity Index #####################################################################################

### * Prepare data ####
data_sites <- sites %>%
  select(id, plot, vegetationCov) %>%
  filter(vegetationCov > 0) %>%
  add_count(plot) %>%
  filter(n == max(n))
data_species <- species %>%
  select(where(~!all(is.na(.x)))) %>%
  mutate(across(where(is.numeric), ~replace(., is.na(.), 0))) %>%
  pivot_longer(-name, names_to = "id", values_to = "value") %>%
  pivot_wider(id, name) %>%
  mutate(year = str_match(id, "\\d{4}"),
         plot = str_match(id, "\\d{2}"),
         year = factor(year),
         plot = factor(plot)) %>%
  arrange(id) %>%
  semi_join(data_sites, by = "id") %>%
  select(plot, year, tidyselect::peek_vars(), -id)

### Separate each year in several tibbles ###

for(i in unique(data_species$year)) {
  nam <- paste("species", i, sep = "")
  
  assign(nam, data_species %>%
           filter(year == i) %>%
           select(-year) %>%
           column_to_rownames(var = "plot")
  )
}

### a Calculate TBI Presence -------------------------------------------------------------------------------------------

#### * 2017 vs. 2018 ####
res1718 <- TBI(species2017, species2018, method = "sorensen", 
               nperm = 9999, test.t.perm = T, clock = T)
res1718$BCD.summary #B = 0.223, C = 0.155, D = 0.378 (58.9% vs. 41.0%)
res1718$t.test_B.C # p.perm = 0.0058
tbi1718 <- as_tibble(res1718$BCD.mat) %>%
  mutate(comparison = "1718")
#### Test plot
plot(res1718, type = "BC")

#### * 2018 vs. 2019 ####
res1819 <- TBI(species2018, species2019, method = "sorensen", 
               nperm = 9999, test.t.perm = T, clock = T)
res1819$BCD.summary #B = 0.118, C = 0.214, D = 0.332 (35.6% vs. 64.3%)
res1819$t.test_B.C # p.perm = 1e-04
tbi1819 <- as_tibble(res1819$BCD.mat)  %>%
  mutate(comparison = "1819")
#### Test plot
plot(res1819, type = "BC")

#### * 2019 vs. 2021 ####
res1921 <- TBI(species2019, species2021, method = "sorensen", 
               nperm = 9999, test.t.perm = T, clock = T)
res1921$BCD.summary #B = 0.249, C = 0.140, D = 0.390 (63.8% vs. 36.1%)
res1921$t.test_B.C # p.perm = 1e-04
tbi1921 <- as_tibble(res1921$BCD.mat) %>%
  mutate(comparison = "1921")
#### Test plot
plot(res1921, type = "BC")

#### * 2017 vs. 2019 ####
res1719 <- TBI(species2017, species2019, method = "sorensen", 
               nperm = 9999, test.t.perm = T, clock = T)
res1719$BCD.summary #B = 0.186, C = 0.212, D = 0.399 (46.7% vs. 53.2%)
res1719$t.test_B.C # p.perm = 0.273
tbi1719 <- as_tibble(res1719$BCD.mat) %>%
  mutate(comparison = "1719")
#### Test plot
plot(res1719, type = "BC")

#### * 2017 vs. 2021 ####
res1721 <- TBI(species2017, species2021, method = "sorensen", 
               nperm = 9999, test.t.perm = T, clock = T)
res1721$BCD.summary #B = 0.184, C = 0.450, D = 0.590 (59.0% vs. 40.9%)
res1721$t.test_B.C # p.perm = 0.0021
tbi1721 <- as_tibble(res1721$BCD.mat) %>%
  mutate(comparison = "1721")
#### Test plot
plot(res1721, type = "BC")

#### * Combine datasets ####
data_presence <- bind_rows(tbi1718, tbi1819, tbi1921, tbi1719, tbi1721) %>%
  mutate(presabu = "presence")

### b Calculate TBI Abundance -------------------------------------------------------------------------------------------

#### * 2017 vs. 2018 ####
res1718 <- TBI(species2017, species2018, method = "%diff", 
               nperm = 9999, test.t.perm = T, clock = T)
res1718$BCD.summary #B = 0.213, C = 0.260, D = 0.473 (45.0% vs. 54.9%)
res1718$t.test_B.C # p.perm = 0.1756
tbi1718 <- as_tibble(res1718$BCD.mat) %>%
  mutate(comparison = "1718")
#### Test plot
plot(res1718, type = "BC")

#### * 2018 vs. 2019 ####
res1819 <- TBI(species2018, species2019, method = "%diff", 
               nperm = 9999, test.t.perm = T, clock = T)
res1819$BCD.summary #B = 0.167, C = 0.302, D = 0.470 (35.7% vs. 64.2%)
res1819$t.test_B.C # p.perm = 1e-04
tbi1819 <- as_tibble(res1819$BCD.mat) %>%
  mutate(comparison = "1819") 
#### Test plot
plot(res1819, type = "BC")

#### * 2019 vs. 2021 ####
res1921 <- TBI(species2019, species2021, method = "%diff", 
               nperm = 9999, test.t.perm = T, clock = T)
res1921$BCD.summary #B = 0.331, C = 0.168, D = 0.499 (66.3% vs. 33.6%)
res1921$t.test_B.C # p.perm = 1e-04
tbi1921 <- as_tibble(res1921$BCD.mat) %>%
  mutate(comparison = "1921")
#### Test plot
plot(res1921, type = "BC")

#### * 2017 vs. 2019 ####
res1719 <- TBI(species2017, species2019, method = "%diff", 
               nperm = 9999, test.t.perm = T, clock = T)
res1719$BCD.summary #B = 0.210, C = 0.390, D = 0.601 (35.0% vs. 64.9%)
res1719$t.test_B.C # p.perm = 1e-04
tbi1719 <- as_tibble(res1719$BCD.mat) %>%
  mutate(comparison = "1719")
#### Test plot
plot(res1719, type = "BC")

#### * 2017 vs. 2021 ####
res1721 <- TBI(species2017, species2021, method = "%diff", 
               nperm = 9999, test.t.perm = T, clock = T)
res1721$BCD.summary #B = 0.301, C = 0.319, D = 0.620 (48.5% vs. 51.4%)
res1721$t.test_B.C # p.perm = 0.598
tbi1721 <- as_tibble(res1721$BCD.mat) %>%
  mutate(comparison = "1721")
#### Test plot
plot(res1721, type = "BC")

#### * Combine datasets ####
data_abundance <- bind_rows(tbi1718, tbi1819, tbi1921, tbi1719, tbi1721) %>%
  mutate(presabu = "abundance")
plot <- data_sites %>%
  filter(str_detect(id, "2017")) %>%
  pull(plot)
data <- add_row(data_presence, data_abundance) %>%
  mutate(plot = rep(plot, length(data_abundance$comparison) * 2 / 38))
tbi <- sites %>%
  filter(surveyYearF == "2017") %>%
  left_join(data, by = "plot") %>%
  rename(B = "B/(2A+B+C)", C = "C/(2A+B+C)", D = "D=(B+C)/(2A+B+C)", change = Change) %>%
  mutate(change = C - B) %>%
  select(id, plot, exposition, side, block, location, locationAbb, locationYear, longitude, latitude, riverkm, distanceRiver, constructionYear, botanist, conf.low, conf.high, PC1soil, PC2soil, PC3soil, PC1constructionYear, PC2constructionYear, PC3constructionYear, B, C, D, comparison, presabu) %>%
  mutate(across(c(PC1soil, PC2soil, PC3soil, 
                  PC1constructionYear, PC2constructionYear, PC3constructionYear,
                  B, C, D), 
                ~ round(.x, digits = 4))) %>%
  mutate(across(c(distanceRiver), 
                ~ round(.x, digits = 1)))
  
rm(list = ls(pattern = "[^species|traits|sites|pcaSoil|pcaSurveyYear|pcaConstructionYear|tbi]"))


## 8 Rounding ##############################################################################################

sites <- sites %>%
  mutate(across(c(NMDS1, NMDS2, permdispPresabs, permdispAbu, lcbd, 
                  syn_total, syn_trend, syn_detrend, 
                  PC1soil, PC2soil, PC3soil, PC1surveyYear, PC2surveyYear, PC3surveyYear, 
                  PC1constructionYear, PC2constructionYear, PC3constructionYear), 
                ~ round(.x, digits = 4))) %>%
  mutate(across(c(NtotalPerc, targetCovratio, graminoidCovratio, targetRichratio, shannon, eveness), 
                ~ round(.x, digits = 3))) %>%
  mutate(across(c(distanceRiver, accumulatedCov), 
                ~ round(.x, digits = 1)))


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# C Save processed data ################################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


### Data ###
write_csv(sites, here("data/processed/data_processed_sites.csv"))
write_csv(species, here("data/processed/data_processed_species.csv"))
write_csv(traits, here("data/processed/data_processed_traits.csv"))
write_csv(tbi, here("data/processed/data_processed_tbi.csv"))

### Tables ###
write_csv(pcaSoil, here("outputs/tables/table_pca_soil.csv"))
write_csv(pcaSurveyYear, here("outputs/tables/table_pca_survey_year.csv"))
write_csv(pcaConstuctionYear, here("outputs/tables/table_pca_construction_year.csv"))
