# Beta diversity on dike grasslands
# Prepare species, sites, and traits data ####
# Markus Bauer
# 2022-01-11



### Packages ###
library(here)
library(tidyverse)
library(naniar) # are_na
library(lubridate) # modify dates
library(vegan) # metaMDS
library(FD) # dbFD
library(adespatial)
#remotes::install_github(file.path("larsito", "tempo"))
library(tempo) # calc_sync
#remotes::install_github(file.path("inbo", "checklist"))

### Start ###
#installr::updateR(browse_news = FALSE, install_R = TRUE, copy_packages = TRUE, copy_Rprofile.site = TRUE, keep_old_packages = TRUE, update_packages = TRUE, start_new_R = FALSE, quit_R = TRUE, print_R_versions = TRUE, GUI = TRUE)
#checklist::setup_source()
#checklist::check_source()
#renv::status()
#renv::snapshot()
#renv::restore()

rm(list = ls())
setwd(here("data", "raw"))



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Load data ##########################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


### 1 Sites ############################################################

sites <- read_csv("data_raw_sites.csv",
  col_names = TRUE,
  na = c("", "NA", "na"), col_types =
    cols(
      .default = "?",
      id = "f",
      block = "f",
      survey_date.2017 = col_date(format = "%Y-%m-%d"),
      survey_date.2018 = col_date(format = "%Y-%m-%d"),
      survey_date.2019 = col_date(format = "%Y-%m-%d"),
      survey_date.2021 = col_date(format = "%Y-%m-%d")
    )) %>%
  mutate(across(
    starts_with("vegetation_cover") |
      starts_with("vegetation_height") |
      starts_with("opensoil_cover") |
      starts_with("moss_cover") |
      starts_with("litter_cover") |
      starts_with("botanist"),
    ~ as.character(.x)
  )) %>%
  pivot_longer(
    starts_with("vegetation_cover.") |
      starts_with("vegetation_height.") |
      starts_with("opensoil_cover.") |
      starts_with("moss_cover.") |
      starts_with("litter_cover.") |
      starts_with("botanist."),
    names_to = c("x", "survey_year"),
    names_sep = "\\.",
    values_to = "n"
    ) %>%
  pivot_wider(names_from = x, values_from = n) %>%
  mutate(across(
    c(survey_year, vegetation_cover, vegetation_height,
      opensoil_cover, moss_cover, litter_cover),
    ~ as.numeric(.x)
  )) %>%
  mutate(
    survey_year_factor = factor(survey_year),
    survey_year_factor_minus = factor(survey_year - 1),
    construction_year_factor = factor(construction_year),
    construction_year_factor_plus = factor(construction_year + 1)
  ) %>%
  mutate(
    id = str_c(id, survey_year, sep = "_"),
    .keep = "all",
    id = paste0("X", id),
    plot = str_sub(id, start = 2, end = 3),
    position = str_sub(id, start = 5, end = 5),
    location_abb = str_sub(location, 1, 3),
    location_abb = str_to_upper(location_abb),
    location_abb = factor(location_abb,
      levels = unique(location_abb[order(construction_year)])
    ),
    location_construction_year = str_c(location_abb, construction_year,
                                       sep = "-")
  ) %>%
  select(
    -position, -starts_with("survey_date_"), -starts_with("topsoil_depth_"),
    -cn_level, -ends_with("_class")
  ) %>%
  relocate(plot,
           .after = id) %>%
  relocate(c("location_abb", "location_construction_year"),
           .after = location) %>%
  relocate(c("survey_year", "survey_year_factor", "survey_year_factor_minus"),
           .after = river_km) %>%
  relocate(c("construction_year_factor", "construction_year_factor_plus"),
           .after = construction_year)

sites_splot <- read_delim(here("data", "raw", "splot", "sPlotOpen_header.txt"),
                         col_names = TRUE, na = c("", "NA", "na"),
                         col_types = cols(
                           .default = "?"
                           ))


### 2 Species #########################################################

species <- data.table::fread("data_raw_species.csv",
  sep = ",",
  dec = ".",
  skip = 0,
  header = TRUE,
  na.strings = c("", "NA", "na"),
  colClasses = list(
    character = "name"
  )
) %>%
  ### Check that each species occurs at least one time ###
  group_by(name) %>%
  arrange(name) %>%
  select(name, all_of(sites$id)) %>%
  mutate(
    total = sum(c_across(starts_with("X")), na.rm = TRUE),
    presence = if_else(total > 0, 1, 0),
    name = factor(name)
  ) %>%
  filter(presence == 1) %>%
  ungroup() %>%
  select(name, sort(tidyselect::peek_vars()), -total, -presence) %>%
  mutate(across(where(is.numeric), ~ replace(., is.na(.), 0)))

species_splot <- read_delim(here("data", "raw", "splot",
                        "sPlotOpen_DT.txt"),
                   col_names = TRUE, na = c("", "NA", "na"), col_types =
                     cols(
                       .default = "?"
                     )) %>%
  filter(Abundance_scale == "CoverPerc")


### Create list with species names and their frequency ###
specieslist <- species %>%
  mutate_if(is.numeric, ~ 1 * (. != 0)) %>%
  mutate(
    sum = rowSums(across(where(is.numeric)), na.rm = TRUE),
    .keep = "unused"
  ) %>%
  group_by(name) %>%
  summarise(sum = sum(sum))
# write_csv(specieslist, "specieslist_2022xxxx.csv")


### 3 Traits ##########################################################

traits <- read_csv("data_raw_traits.csv",
  col_names = TRUE, na = c("", "NA", "na"),
  col_types =
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
    )
) %>%
  separate(name, c("genus", "species", "ssp", "subspecies"), "_",
    remove = FALSE, extra = "drop", fill = "right"
  ) %>%
  mutate(
    genus = str_sub(genus, 1, 4),
    species = str_sub(species, 1, 4),
    subspecies = str_sub(subspecies, 1, 4),
    name = factor(name)
  ) %>%
  unite(abb, genus, species, subspecies, sep = "") %>%
  mutate(
    abb = str_replace(abb, "NA", ""),
    abb = as_factor(abb)
  ) %>%
  select(-ssp) %>%
  arrange(name)

### Check congruency of traits and species table ###
traits[duplicated(traits$abb), ]
# traits$name[which(!(traits$name %in% species$name))]
species$name[which(!(species$name %in% traits$name))]

### Combine with species table ###
traits <- traits %>%
  semi_join(species, by = "name")


### 4 Temperature #######################################################

sites_temperature <- read_csv(here("data", "raw", "temperature", "data",
                             "data_OBS_DEU_P1M_T2M.csv"),
                 col_names = TRUE, na = c("", "NA", "na"), col_types =
                   cols(
                     .default = "?"
                   )) %>%
  rename(date = Zeitstempel, value = Wert, site = SDO_ID) %>%
  select(site, date, value)

  
### 5 Precipitation #######################################################

sites_precipitation <-  read_csv(here("data", "raw", "precipitation", "data",
                                "data_OBS_DEU_P1M_RR.csv"),
           col_names = TRUE, na = c("", "NA", "na"), col_types =
             cols(
               .default = "?"
             )
  ) %>%
  rename(date = Zeitstempel, value = Wert, site = SDO_ID) %>%
  select(site, date, value)
  

### 6 Check data frames ################################################

### Check typos ###
sites %>%
  filter(!str_detect(id, "_seeded$")) %>%
  janitor::tabyl(vegetation_cover)
# sites %>% filter(vegetation_cover == 17)
species %>%
  select(-name) %>%
  unlist() %>%
  janitor::tabyl()
species %>% # Check special typos
  pivot_longer(-name, names_to = "id", values_to = "value") %>%
  filter(value == 8)

### Compare vegetation_cover and accumulated_cover ###
species %>%
  summarise(across(where(is.double), ~ sum(.x, na.rm = TRUE))) %>%
  pivot_longer(cols = everything(), names_to = "id", values_to = "value") %>%
  mutate(id = factor(id)) %>%
  full_join(sites, by = "id") %>%
  mutate(diff = (value - vegetation_cover)) %>%
  select(id, survey_year, value, vegetation_cover, diff) %>%
  filter(diff > 50 | diff < -30) %>%
  arrange(diff) %>%
  print(n = 100)

### Check plots over time ###
species %>%
  select(name, starts_with("X01")) %>%
  filter(if_any(starts_with("X"), ~ . > 0)) %>%
  print(n = 70)

### Check missing data ###
miss_var_summary(sites, order = TRUE)
vis_miss(sites, cluster = FALSE, sort_miss = TRUE)
vis_miss(traits, cluster = FALSE, sort_miss = TRUE)


rm(list = setdiff(ls(), c("sites", "sites_splot",
                          "sites_precipitation", "sites_temperature",
                          "species", "species_splot",
                          "traits")))



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Create variables ###########################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


## 1 Create simple variables ###################################################

traits <- traits %>%
  mutate(
    lean_indicator = if_else(
      !(is.na(table30)) | !(is.na(table33)) | !(is.na(table34)), "yes", "no"
    ),
    target = if_else(
      target_herb == "yes" | target_grass == "yes", "yes", "no"
    ),
    ruderal = if_else(
      sociology >= 3300 & sociology < 3700, "yes", "no"
    ),
    target_ellenberg = if_else(
      sociology >= 5300 & sociology < 5400, "dry_grassland", if_else(
        sociology >= 5400 & sociology < 6000, "hay_meadow", if_else(
          sociology >= 5100 & sociology < 5200, "nardus_grassland", if_else(
            sociology >= 5200 & sociology < 5300, "sand_grassland", "no"
          )
        )
      )
    )
  )

sites <- sites %>%
  mutate(
    n_total_concentration = finematerial_depth * finematerial_density * 10 *
      n_total_ratio / 100,
    plot_age = survey_year - construction_year
  )


## 2 Coverages #################################################################

cover <- species %>%
  left_join(traits, by = "name") %>%
  select(
    name, family, target, target_herb, target_arrhenatherion, lean_indicator,
    nitrogen_indicator, ruderal, table33, target_ellenberg,
    starts_with("X")
  ) %>%
  pivot_longer(names_to = "id", values_to = "n", cols = starts_with("X")) %>%
  group_by(id)

### * graminoid, herb, and total coverage ####
cover_total_and_graminoid <- cover %>%
  group_by(id, family) %>%
  summarise(total = sum(n, na.rm = TRUE), .groups = "keep") %>%
  mutate(type = if_else(family == "Poaceae" |
                          family == "Cyperaceae" |
                          family == "Juncaceae",
    "graminoid_cover", "herb_cover"
  )) %>%
  group_by(id, type) %>%
  summarise(total = sum(total, na.rm = TRUE), .groups = "keep") %>%
  spread(type, total) %>%
  mutate(accumulated_cover = graminoid_cover + herb_cover) %>%
  ungroup()

### * German biotope mapping: Target species' coverage ####
cover_target <- cover %>%
  filter(target == "yes") %>%
  summarise(target_cover = sum(n, na.rm = TRUE)) %>%
  mutate(target_cover = round(target_cover, 1)) %>%
  ungroup()

### * German biotope mapping: Target herb species' coverage ####
cover_target_herb <- cover %>%
  filter(target_herb == "yes") %>%
  summarise(target_herb_cover = sum(n, na.rm = TRUE)) %>%
  mutate(target_herb_cover = round(target_herb_cover, 1)) %>%
  ungroup()

### * German biotope mapping: Arrhenatherum species' cover ratio ####
cover_target_arrhenatherion <- cover %>%
  filter(target_arrhenatherion == "yes") %>%
  summarise(arrh_cover = sum(n, na.rm = TRUE)) %>%
  mutate(arrh_cover = round(arrh_cover, 1)) %>%
  ungroup()

### * German biotope mapping: Lean indicator's coverage ####
cover_lean_indicator <- cover %>%
  filter(lean_indicator == "yes") %>%
  summarise(lean_cover = sum(n, na.rm = TRUE)) %>%
  mutate(lean_cover = round(lean_cover, 1)) %>%
  ungroup()

### * German biotope mapping: Nitrogen indicator's coverage ####
cover_nitrogen_indicator <- cover %>%
  filter(nitrogen_indicator == "yes") %>%
  summarise(nitrogen_cover = sum(n, na.rm = TRUE)) %>%
  mutate(nitrogen_cover = round(nitrogen_cover, 1)) %>%
  ungroup()

### * Ruderal's coverage ####
cover_ruderal <- cover %>%
  filter(ruderal == "yes") %>%
  summarise(ruderal_cover = sum(n, na.rm = TRUE)) %>%
  mutate(ruderal_cover = round(ruderal_cover, 1)) %>%
  ungroup()

### * German biotope mapping: Table 33 species' coverage ####
cover_table33 <- cover %>%
  mutate(table33 = if_else(table33 == "4" |
    table33 == "3" |
    table33 == "2",
  "table33_cover", "other"
  )) %>%
  filter(table33 == "table33_cover") %>%
  summarise(table33_cover = sum(n, na.rm = TRUE)) %>%
  mutate(table33_cover = round(table33_cover, 1)) %>%
  ungroup()

### * Ellenberg target species' coverage ####
cover_ellenberg <- cover %>%
  mutate(ellenberg = if_else(target_ellenberg == "dry_grassland" |
                          target_ellenberg == "hay_meadow",
                        "ellenberg_cover", "other"
  )) %>%
  filter(ellenberg == "ellenberg_cover") %>%
  summarise(ellenberg_cover = sum(n, na.rm = TRUE)) %>%
  mutate(ellenberg_cover = round(ellenberg_cover, 1)) %>%
  ungroup()

### * implement in sites data set ####
sites <- sites %>%
  right_join(cover_total_and_graminoid, by = "id") %>%
  right_join(cover_target, by = "id") %>%
  right_join(cover_target_herb, by = "id") %>%
  right_join(cover_target_arrhenatherion, by = "id") %>%
  right_join(cover_lean_indicator, by = "id") %>%
  right_join(cover_nitrogen_indicator, by = "id") %>%
  right_join(cover_ruderal, by = "id") %>%
  right_join(cover_table33, by = "id") %>%
  right_join(cover_ellenberg, by = "id") %>%
  ### Calcute the ratio of target species richness of total species richness
  mutate(
    target_cover_ratio = target_cover / accumulated_cover,
    graminoid_cover_ratio = graminoid_cover / accumulated_cover
  )

rm(list = setdiff(ls(), c("sites", "sites_splot",
                          "sites_precipitation", "sites_temperature",
                          "species", "species_splot",
                          "traits")))


## 3 Alpha diversity ###################################################

### a Species richness -------------------------------------------------

richness <- left_join(species, traits, by = "name") %>%
  select(
    name, rlg, rlb, target, target_herb, target_arrhenatherion,
    ffh6510, ffh6210, nitrogen_indicator, lean_indicator, table33, table34,
    starts_with("X")
  ) %>%
  pivot_longer(names_to = "id", values_to = "n", cols = starts_with("X")) %>%
  mutate(n = if_else(n > 0, 1, 0)) %>%
  group_by(id)

### * total species richness ####
richness_total <- richness %>%
  summarise(species_richness = sum(n, na.rm = TRUE)) %>%
  ungroup()

### * red list Germany (species richness) ####
richness_rlg <- richness %>%
  filter(rlg == "1" | rlg == "2" | rlg == "3" | rlg == "V") %>%
  summarise(rlg_richness = sum(n, na.rm = TRUE)) %>%
  ungroup()

### * red list Bavaria (species richness) ####
richness_rlb <- richness %>%
  filter(rlb == "1" | rlb == "2" | rlb == "3" | rlb == "V") %>%
  summarise(rlb_richness = sum(n, na.rm = TRUE)) %>%
  ungroup()

### * German biotope mapping: target species (species richness) ####
richness_target <- richness %>%
  filter(target == "yes") %>%
  summarise(target_richness = sum(n, na.rm = TRUE)) %>%
  ungroup()

### * German biotope mapping: target herb species (species richness) ####
richness_target_herb <- richness %>%
  filter(target_herb == "yes") %>%
  summarise(target_herb_richness = sum(n, na.rm = TRUE)) %>%
  ungroup()

### * German biotope mapping: Arrhenatherion species (species richness) ####
richness_arrh <- richness %>%
  filter(target_arrhenatherion == "yes") %>%
  summarise(arrh_richness = sum(n, na.rm = TRUE)) %>%
  ungroup()

### * ffh6510 species (species richness) ####
richness_ffh6510 <- richness %>%
  filter(ffh6510 == "yes") %>%
  summarise(ffh6510_richness = sum(n, na.rm = TRUE)) %>%
  ungroup()

### * ffh6210 species (species richness) ####
richness_ffh6210 <- richness %>%
  filter(ffh6210 == "yes") %>%
  summarise(ffh6210_richness = sum(n, na.rm = TRUE)) %>%
  ungroup()

### * German biotope mapping: lean_indicator species (species richness) ####
richness_lean_indicator <- richness %>%
  filter(lean_indicator == "yes") %>%
  summarise(lean_indicator_richness = sum(n, na.rm = TRUE)) %>%
  ungroup()

### * German biotope mapping: table33 (species richness) ####
richness_table33_2 <- richness %>%
  filter(table33 == "2") %>%
  summarise(table33_2_richness = sum(n, na.rm = TRUE)) %>%
  ungroup()
richness_table33_3 <- richness %>%
  filter(table33 == "3") %>%
  summarise(table33_3_richness = sum(n, na.rm = TRUE)) %>%
  ungroup()
richness_table33_4 <- richness %>%
  filter(table33 == "4") %>%
  summarise(table33_4_richness = sum(n, na.rm = TRUE)) %>%
  ungroup()

### * German biotope mapping: table34 (species richness) ####
richness_table34_2 <- richness %>%
  filter(table34 == "2") %>%
  summarise(table34_2_richness = sum(n, na.rm = TRUE)) %>%
  ungroup()
richness_table34_3 <- richness %>%
  filter(table34 == "3") %>%
  summarise(table34_3_richness = sum(n, na.rm = TRUE)) %>%
  ungroup()

### * implement in sites data set (species richness) ####
sites <- sites %>%
  right_join(richness_total, by = "id") %>%
  right_join(richness_rlg, by = "id") %>%
  right_join(richness_rlb, by = "id") %>%
  right_join(richness_target, by = "id") %>%
  right_join(richness_target_herb, by = "id") %>%
  right_join(richness_arrh, by = "id") %>%
  right_join(richness_ffh6510, by = "id") %>%
  right_join(richness_ffh6210, by = "id") %>%
  right_join(richness_lean_indicator, by = "id") %>%
  right_join(richness_table33_2, by = "id") %>%
  right_join(richness_table33_3, by = "id") %>%
  right_join(richness_table33_4, by = "id") %>%
  right_join(richness_table34_2, by = "id") %>%
  right_join(richness_table34_3, by = "id") %>%
  ### Calcute the ratio of target species richness of total species richness
  mutate(target_richness_ratio = target_richness / species_richness)

### b Species eveness and shannon --------------------------------------

data <- species %>%
  mutate(across(where(is.numeric), ~ replace(., is.na(.), 0))) %>%
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
  mutate(eveness = shannon / log(species_richness))

rm(list = setdiff(ls(), c("sites", "sites_splot",
                          "sites_precipitation", "sites_temperature",
                          "species", "species_splot",
                          "traits")))


## 4 Biotope types #####################################################

### a Calculate types --------------------------------------------------

biotopetypes <- sites %>%
  select(
    id, table33_2_richness, table33_3_richness, table33_4_richness,
    table33_cover, table34_2_richness, table34_3_richness, target_richness,
    target_herb_richness, arrh_richness, target_cover, lean_cover, arrh_cover,
    target_herb_cover, nitrogen_cover
  ) %>%
  mutate(
    table33Rich_proof = if_else(
      table33_2_richness >= 2 | table33_2_richness + table33_3_richness >= 3 |
        table33_2_richness + table33_3_richness + table33_4_richness >= 4,
      "yes", "no"
    ),
    table33_cover_proof = if_else(
      table33_cover >= 25,
      "yes", "no"
    ),
    table33_proof = if_else(
      table33Rich_proof == "yes" & table33_cover_proof == "yes",
      "yes", "no"
    ),
    table34Rich_proof = if_else(
      table34_2_richness >= 2 | table34_2_richness + table34_3_richness >= 3,
      "yes", "no"
    ),
    G312_GT6210_type = if_else(
      table33_proof == "yes" & table34Rich_proof == "yes",
      "yes", "no"
    ),
    GE_proof = if_else(
      target_herb_cover >= 12.5 & target_richness >= 20 & nitrogen_cover < 25,
      "yes", "no"
    ),
    G214_GE6510_type = if_else(
      GE_proof == "yes" &
        arrh_richness >= 1 &
        arrh_cover > 0.5 &
        lean_cover >= 25 &
        target_herb_cover >= 12.5,
      "yes", "no"
    ),
    G212_LR6510_type = if_else(
      GE_proof == "yes" &
        arrh_richness >= 1 &
        arrh_cover > 0.5 &
        lean_cover < 25,
      "yes", "no"
    ),
    G214_GE00BK_type = if_else(
      GE_proof == "yes" &
        arrh_richness == 0 &
        lean_cover >= 25 &
        target_herb_cover >= 12.5,
      "yes", if_else(
        GE_proof == "yes" &
          arrh_cover >= 0.5 &
          lean_cover >= 25 &
          target_herb_cover >= 12.5,
        "yes", "no"
      )
    ),
    G213_GE00BK_type = if_else(
      GE_proof == "yes" &
        arrh_richness == 0 &
        lean_cover >= 25 &
        target_herb_cover < 12.5,
      "yes", if_else(
        GE_proof == "yes" &
          arrh_cover >= 0.5 &
          lean_cover >= 25 &
          target_herb_cover < 12.5,
        "yes", "no"
      )
    ),
    G213_type = if_else(
      lean_cover >= 25,
      "yes", "no"
    ),
    G212_type = if_else(
      lean_cover >= 1 &
        lean_cover < 25 &
        target_herb_cover >= 12.5 &
        target_herb_richness >= 10,
      "yes", "no"
    ),
    G211_type = if_else(
      lean_cover >= 1 &
        lean_cover < 25 &
        target_herb_cover >= 1 &
        target_herb_richness >= 5,
      "yes", "no"
    ),
    biotopeType = if_else(
      G312_GT6210_type == "yes", "G312-GT6210", if_else(
        G214_GE6510_type == "yes", "G214-GE6510", if_else(
          G212_LR6510_type == "yes", "G212-LR6510", if_else(
            G214_GE00BK_type == "yes", "G214-GE00BK", if_else(
              G213_GE00BK_type == "yes", "G213-GE00BK", if_else(
                G213_type == "yes", "G213", if_else(
                  G212_type == "yes", "G212", if_else(
                    G211_type == "yes",
                    "G211", "other"
                  )
                )
              )
            )
          )
        )
      )
    ),
    biotopeType = as_factor(biotopeType)
  ) %>%
  select(id, biotopeType, -ends_with("proof"), -starts_with("table")) %>%
  mutate(
    ffh6510 = str_match(biotopeType, "6510"),
    ffh6210 = str_match(biotopeType, "6210"),
    baykompv = as_factor(str_sub(biotopeType, start = 1, end = 4))
  ) %>%
  unite(ffh, ffh6510, ffh6210, sep = "") %>%
  mutate(
    ffh = str_replace(ffh, "NA", ""),
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
                    biotopeType == "G211",
                    6, 0
                  )
                )
              )
            )
          )
        )
      )
    ),
    min8 = as_factor(if_else(biotopePoints >= 8, "yes", "no")),
    min9 = as_factor(if_else(biotopePoints >= 9, "yes", "no"))
  )
sites <- sites %>%
  left_join(biotopetypes, by = "id") %>%
  select(
    -target_herb_cover, -arrh_cover, -lean_cover, -nitrogen_cover, -table33_cover,
    -target_herb_richness, -arrh_richness, -lean_indicator_richness,
    -ffh6510_richness, -ffh6210_richness, -table33_2_richness,
    -table33_3_richness, -table33_4_richness, -table34_2_richness,
    -table34_3_richness
  )

traits <- traits %>%
  select(
    -target_arrhenatherion, -table30, -table33, -table34,
    -nitrogen_indicator, -nitrogen_indicator2, -lean_indicator,
    -grazing_indicator, -rlb, -target_herb, -target_grass,
    -sociology, -legal
  )

### b Calculate constancy of ffh evaluation -----------------------------------

data <- sites %>%
  select(id, plot, survey_year, ffh) %>%
  group_by(plot) %>%
  mutate(count = n()) %>%
  filter(count == max(count)) %>%
  pivot_wider(id_cols = -id,
              names_from = "survey_year",
              values_from = "ffh") %>% # group_by(plot) %>%
  rename(x17 = "2017", x18 = "2018", x19 = "2019") %>%
  mutate(
    changeType = ifelse(
      (x17 == "non-FFH" & x18 != "non-FFH" & x19 != "non-FFH"),
      "better", ifelse(
        (x17 != "non-FFH" & x18 == "non-FFH" & x19 != "non-FFH") |
          (x17 == "non-FFH" & x18 != "non-FFH" & x19 == "non-FFH") |
          (x17 == "non-FFH" & x18 == "non-FFH" & x19 != "non-FFH") |
          (x17 != "non-FFH" & x18 != "non-FFH" & x19 == "non-FFH"),
        "change", ifelse(
          (x17 == "non-FFH" & x18 == "non-FFH" & x19 == "non-FFH"),
          "non-FFH", ifelse(
            x17 == "6510" & x18 == "6510" & x19 == "6510",
            "FFH6510", ifelse(
              x17 == "6210" & x18 == "6210" & x19 == "6210",
              "FFH6210", ifelse(
                x17 != "non-FFH" & x18 != "non-FFH" & x19 != "non-FFH",
                "any-FFH", "worse"
              )
            )
          )
        )
      )
      )
    ) %>%
  select(plot, changeType)
sites <- sites %>%
  left_join(data, by = "plot")

rm(list = setdiff(ls(), c("sites", "sites_splot",
                          "sites_precipitation", "sites_temperature",
                          "species", "species_splot",
                          "traits")))


## 5 Beta diversity ###################################################


### a dbMEM (41 plots) ------------------------------------------------

### * Prepare data ####
source("https://raw.githubusercontent.com/zdealveindy/anadat-r/master/scripts/NumEcolR2/quickMEM.R")
data_sites <- sites %>%
  filter(accumulated_cover > 0) %>%
  add_count(plot) %>%
  filter(n == max(n))
data_species <- species %>%
  pivot_longer(-name, names_to = "id", values_to = "value") %>%
  pivot_wider(id, names_from = "name", values_from = "value") %>%
  arrange(id) %>%
  semi_join(data_sites, by = "id") %>%
  mutate(across(where(is.numeric), ~ replace(., is.na(.), 0)))

### * 2017 ####
data_sites_dbmem <- data_sites %>%
  filter(survey_year == 2017) %>%
  select(id, longitude, latitude)
data_species_dbmem <- data_species %>%
  semi_join(data_sites_dbmem, by = "id") %>%
  column_to_rownames(var = "id") %>%
  decostand("hellinger")
i <- (colSums(data_species_dbmem, na.rm = T) != 0)
data_species_dbmem <- data_species_dbmem[, i]
data_sites_dbmem <- data_sites_dbmem %>%
  column_to_rownames("id")
m <- quickMEM(data_species_dbmem, data_sites_dbmem,
  alpha = 0.05,
  method = "fwd",
  rangexy = TRUE,
  perm.max = 999
)
# --> undetrended data, R2adj of minimum (final) model = 0.034
m$RDA_test # p = 0.002
m$RDA_axes_test # RDA1 sig axes
m$RDA # PC1 = 0.114
dbmem_reduced <- m$dbMEM %>%
  rownames_to_column(var = "id") %>%
  select(id, MEM1) %>%
  rename(mem1_2017 = MEM1)
sites <- left_join(sites, dbmem_reduced, by = "id")

### * 2018 ####
data_sites_dbmem <- data_sites %>%
  filter(survey_year == 2018) %>%
  select(id, longitude, latitude)
data_species_dbmem <- data_species %>%
  semi_join(data_sites_dbmem, by = "id") %>%
  column_to_rownames(var = "id") %>%
  decostand("hellinger")
i <- (colSums(data_species_dbmem, na.rm = T) != 0)
data_species_dbmem <- data_species_dbmem[, i]
data_sites_dbmem <- data_sites_dbmem %>%
  column_to_rownames("id")
m1 <- quickMEM(data_species_dbmem, data_sites_dbmem,
  alpha = 0.05,
  method = "fwd",
  rangexy = TRUE,
  perm.max = 999
)
# --> undetrended data, R2adj of minimum (final) model = 0.036
m1$RDA_test # p = 0.006
m1$RDA_axes_test #  RDA2 sig. axes
m1$RDA # PC1 = 0.109, PC2 = 0.082
dbmem_reduced <- m1$dbMEM %>%
  rownames_to_column(var = "id") %>%
  select(id, MEM2) %>%
  rename(mem2_2018 = MEM2)
sites <- left_join(sites, dbmem_reduced, by = "id")

### * 2019 ####
data_sites_dbmem <- data_sites %>%
  filter(survey_year == 2019) %>%
  select(id, longitude, latitude)
data_species_dbmem <- data_species %>%
  semi_join(data_sites_dbmem, by = "id") %>%
  column_to_rownames(var = "id") %>%
  decostand("hellinger")
i <- (colSums(data_species_dbmem, na.rm = T) != 0)
data_species_dbmem <- data_species_dbmem[, i]
data_sites_dbmem <- data_sites_dbmem %>%
  column_to_rownames("id")
m <- quickMEM(data_species_dbmem, data_sites_dbmem,
  alpha = 0.05,
  method = "fwd",
  detrend = FALSE,
  rangexy = TRUE,
  perm.max = 999
)
# --> not detrended, R2adj of minimum (final) model = 0.056
m$RDA_test # p = 0.001
m$RDA_axes_test # RDA1 sig axes
m$RDA # PC1 = 0.080
dbmem_reduced <- m$dbMEM %>%
  rownames_to_column(var = "id") %>%
  select(id, MEM1) %>%
  rename(mem1_2019 = MEM1)
sites <- left_join(sites, dbmem_reduced, by = "id")

### * 2021 ####
data_sites_dbmem <- data_sites %>%
  filter(survey_year == 2021) %>%
  select(id, longitude, latitude)
data_species_dbmem <- data_species %>%
  semi_join(data_sites_dbmem, by = "id") %>%
  column_to_rownames(var = "id") %>%
  decostand("hellinger")
i <- (colSums(data_species_dbmem, na.rm = T) != 0)
data_species_dbmem <- data_species_dbmem[, i]
data_sites_dbmem <- data_sites_dbmem %>%
  column_to_rownames("id")
m <- quickMEM(data_species_dbmem, data_sites_dbmem,
  alpha = 0.05,
  method = "fwd",
  detrend = FALSE,
  rangexy = TRUE,
  perm.max = 999
)
# --> R2adj of minimum (final) model = 0.064
m$RDA_test # p = 0.001
m$RDA_axes_test # RDA1, RDA2 sig axes
m$RDA # PC1 = 0.096, PC2 = 0.085
dbmem_reduced <- m$dbMEM %>%
  rownames_to_column(var = "id") %>%
  select(id, MEM1, MEM2) %>%
  rename(mem1_2021 = MEM1, mem2_2021 = MEM2)
sites <- left_join(sites, dbmem_reduced, by = "id")

### b LCBD (Local Contributions to Beta Diversity) --------------------

#### * 2017 ####
data_species_lcbd <- data_species %>%
  filter(str_detect(id, "2017")) %>%
  column_to_rownames(var = "id")
data <- beta.div(data_species_lcbd, method = "hellinger", nperm = 9999)
data_2017 <- data$LCBD
row.names(data_species_lcbd[which(p.adjust(data$p.LCBD, "holm") <= 0.05), ])
#--> X62_m_2017

#### * 2018 ####
data_species_lcbd <- data_species %>%
  filter(str_detect(id, "2018")) %>%
  column_to_rownames(var = "id")
data <- beta.div(data_species_lcbd, method = "hellinger", nperm = 9999)
data_2018 <- data$LCBD
row.names(data_species_lcbd[which(p.adjust(data$p.LCBD, "holm") <= 0.05), ])
#--> none

#### * 2019 ####
data_species_lcbd <- data_species %>%
  filter(str_detect(id, "2019")) %>%
  column_to_rownames(var = "id")
data <- beta.div(data_species_lcbd, method = "hellinger", nperm = 9999)
data_2019 <- data$LCBD
row.names(data_species_lcbd[which(p.adjust(data$p.LCBD, "holm") <= 0.05), ])
#--> none

#### * 2021 ####
data_species_lcbd <- data_species %>%
  filter(str_detect(id, "2021")) %>%
  column_to_rownames(var = "id")
data <- beta.div(data_species_lcbd, method = "hellinger", nperm = 9999)
data_2021 <- data$LCBD
row.names(data_species_lcbd[which(p.adjust(data$p.LCBD, "holm") <= 0.05), ])
#--> none

#### * combine datasets ####
data <- c(data_2017, data_2018, data_2019, data_2021)
data <- data_sites %>%
  mutate(lcbd = data) %>%
  select(id, lcbd)
sites <- sites %>%
  left_join(data, by = "id")

rm(list = setdiff(ls(), c("sites", "sites_splot",
                          "sites_precipitation", "sites_temperature",
                          "species", "species_splot",
                          "traits")))

### c Synchrony --------------------------------------------------------

data_sites <- sites %>%
  select(id, plot, vegetation_cover) %>%
  filter(vegetation_cover > 0) %>%
  add_count(plot) %>%
  filter(n == max(n))
data_species <- species %>%
  select(where(~ !all(is.na(.x)))) %>%
  mutate(across(where(is.numeric), ~ replace(., is.na(.), 0))) %>%
  pivot_longer(-name, names_to = "id", values_to = "value") %>%
  pivot_wider(id, name) %>%
  mutate(
    year = str_match(id, "\\d{4}"),
    plot = str_match(id, "\\d{2}"),
    year = factor(year),
    plot = factor(plot)
  ) %>%
  arrange(id) %>%
  semi_join(data_sites, by = "id") %>%
  column_to_rownames(var = "id") %>%
  select(plot, year, tidyselect::peek_vars())
data <- data_species %>%
  split(data_species$plot, drop = TRUE) %>%
  map(~ (.x %>% select(-plot, -year))) %>%
  map(~ (.x %>% select(where(~ !all(is.na(.x))))))
sync_indices <- map(data, calc_sync)
data <- do.call("rbind", sync_indices) %>%
  # map() does not work because 'plot' is missing
  rownames_to_column("plot") %>%
  as_tibble() %>%
  mutate(plot = str_extract(plot, "\\d\\d")) %>%
  select(plot, syn_total, syn_trend, syn_detrend)
#--> log_varrrat_t3 does not allow missing years
sites <- left_join(sites, data, by = "plot")

rm(list = setdiff(ls(), c("sites", "sites_splot",
                          "sites_precipitation", "sites_temperature",
                          "species", "species_splot",
                          "traits")))



## 6 Environmental variables ###########################################

### a Soil PCA  --------------------------------------------------------

### Prepare data ###
data <- sites %>%
  filter(accumulated_cover > 0) %>%
  add_count(plot) %>%
  filter(n == max(n) & survey_year == 2017) %>%
  select(
    plot, ph,
    calciumcarbonat, humus, n_total_ratio, n_total_concentration, cn_ratio,
    sand, silt, clay,
    phosphorus, potassium, magnesium,
    topsoil_depth
  ) %>%
  column_to_rownames(var = "plot")
### Calculate PCA ###
pca <- rda(data, scale = TRUE)
summary(pca, axes = 0)
screeplot(pca, bstick = TRUE, npcs = length(pca$CA$eig))
biplot(pca, display = "species", scaling = 2)
### create summary table ###
eigenvals <- pca %>%
  eigenvals() %>%
  summary() %>%
  as.data.frame() %>%
  rownames_to_column(var = "variables") %>%
  tibble() %>%
  select(PC1:PC4, variables)
data <- pca %>%
  summary() %>%
  magrittr::extract2("species") %>%
  as.data.frame() %>%
  rownames_to_column(var = "variables") %>%
  tibble() %>%
  select(PC1:PC4, variables)
pca_soil <- data %>%
  bind_rows(eigenvals)
### Add to sites ###
data <- pca %>%
  summary() %>%
  magrittr::extract2("sites") %>%
  as.data.frame() %>%
  rownames_to_column(var = "plot") %>%
  tibble() %>%
  select(plot, PC1:PC4) %>%
  rename(
    pc1_soil = PC1,
    pc2_soil = PC2,
    pc3_soil = PC3,
    pc4_soil = PC4
  )
sites <- sites %>%
  left_join(data, by = "plot")

rm(list = setdiff(ls(), c("sites", "sites_splot",
                          "sites_precipitation", "sites_temperature",
                          "species", "species_splot",
                          "traits", "pca_soil")))


### b Climate PCA  -----------------------------------------------------

### * Temperature ####
data <- sites_temperature %>%
  filter(date >= "2002-03-01") %>%
  mutate(
    site = factor(site),
    season = floor_date(date, "season"),
    year = year(season),
    season = month(season),
    season = factor(season),
    season = fct_recode(season,
      "spring" = "3",
      "summer" = "6",
      "autumn" = "9",
      "winter" = "12"
    )
  ) %>%
  group_by(year) %>%
  mutate(
    yearMean = mean(value),
    yearMean = round(yearMean, digits = 1),
    currentYear = if_else(season == "spring", 0, 1),
    currentYear = year + currentYear
  ) %>%
  #--> current year of 2021 is e.g. from summer 2020 to spring 2021 = climate for survey_year
  group_by(currentYear) %>%
  mutate(
    currentMean = mean(value),
    currentMean = round(currentMean, digits = 1),
    currentMean = if_else(season == "spring", currentMean, NA_real_)
  ) %>%
  group_by(year, season, yearMean, currentYear, currentMean) %>%
  summarise(seasonMean = round(mean(value), digits = 1), .groups = "keep") %>%
  pivot_wider(
    id_cols = c(year, yearMean, currentMean),
    names_from = season, values_from = seasonMean
  ) %>%
  group_by(year) %>%
  summarise(across(where(is.numeric), ~ max(., na.rm = TRUE)),
            .groups = "keep") %>%
  # warnings because of lates year (summer, autumn, winter), can be ignored
  mutate(year = factor(year))

sites <- sites %>%
  left_join(data %>% select(-currentMean),
    by = c("construction_year_factor" = "year")
  ) %>%
  rename(
    "tempSpring_construction_year" = "spring",
    "tempSummer_construction_year" = "summer",
    "tempAutumn_construction_year" = "autumn",
    "tempWinter_construction_year" = "winter",
    "tempMean_construction_year" = "yearMean"
  ) %>%
  left_join(data %>% select(-currentMean),
    by = c("construction_year_factor_plus" = "year")
  ) %>%
  rename(
    "tempSpring_construction_year_plus" = "spring",
    "tempSummer_construction_year_plus" = "summer",
    "tempAutumn_construction_year_plus" = "autumn",
    "tempWinter_construction_year_plus" = "winter",
    "tempMean_construction_year_plus" = "yearMean"
  ) %>%
  left_join(data %>% select(year, currentMean, spring),
    by = c("survey_year_factor" = "year")
  ) %>%
  rename(
    "tempMean_survey_year" = "currentMean",
    "tempSpring_survey_year" = "spring"
  ) %>%
  left_join(data %>% select(year, summer, autumn, winter),
    by = c("survey_year_factor_minus" = "year")
  ) %>%
  rename(
    "tempSummer_survey_year" = "summer",
    "tempAutumn_survey_year" = "autumn",
    "tempWinter_survey_year" = "winter"
  )

### * Precipitation ####
data <- sites_precipitation %>%
  filter(date >= "2002-03-01") %>%
  mutate(
    site = factor(site),
    season = floor_date(date, "season"),
    year = year(season),
    season = month(season),
    season = factor(season),
    season = fct_recode(season,
      "spring" = "3",
      "summer" = "6",
      "autumn" = "9",
      "winter" = "12"
    )
  ) %>%
  group_by(site, year) %>%
  mutate(
    yearSum = round(sum(value), digits = 0),
    currentYear = if_else(season == "spring", 0, 1),
    currentYear = year + currentYear
  ) %>%
  # current year of 2021 is e.g. from summer 2020 to spring 2021 = climate for survey_year
  group_by(site, currentYear) %>%
  mutate(
    currentSum = sum(value),
    currentSum = round(currentSum, 0)
  ) %>%
  group_by(site, season, year, currentYear, yearSum, currentSum) %>%
  summarise(
    seasonSum = round(sum(value), digits = 0),
    .groups = "keep"
  ) %>%
  group_by(year, season, currentYear) %>%
  summarise(
    seasonMean = round(mean(seasonSum), digits = 0),
    yearMean = round(mean(yearSum), digits = 0),
    currentMean = round(mean(currentSum), digits = 0),
    .groups = "keep"
  ) %>%
  mutate(currentMean = if_else(season == "spring", currentMean, NA_real_)) %>%
  pivot_wider(id_cols = c(year, yearMean, currentMean),
              names_from = season, values_from = seasonMean) %>%
  group_by(year) %>%
  summarise(across(where(is.numeric), ~ max(., na.rm = TRUE)),
    .groups = "keep"
  ) %>%
  # warnings because of lates year (summer, autumn, winter), can be ignored
  mutate(year = factor(year))
sites <- sites %>%
  left_join(data %>% select(-currentMean),
    by = c("construction_year_factor" = "year")
  ) %>%
  rename(
    "precSpring_construction_year" = "spring",
    "precSummer_construction_year" = "summer",
    "precAutumn_construction_year" = "autumn",
    "precWinter_construction_year" = "winter",
    "precMean_construction_year" = "yearMean"
  ) %>%
  left_join(data %>% select(-currentMean),
    by = c("construction_year_factor_plus" = "year")
  ) %>%
  rename(
    "precSpring_construction_year_plus" = "spring",
    "precSummer_construction_year_plus" = "summer",
    "precAutumn_construction_year_plus" = "autumn",
    "precWinter_construction_year_plus" = "winter",
    "precMean_construction_year_plus" = "yearMean"
  ) %>%
  left_join(data %>% select(year, currentMean, spring),
    by = c("survey_year_factor" = "year")
  ) %>%
  rename(
    "precMean_survey_year" = "currentMean",
    "precSpring_survey_year" = "spring"
  ) %>%
  left_join(data %>% select(year, summer, autumn, winter),
    by = c("survey_year_factor_minus" = "year")
  ) %>%
  rename(
    "precSummer_survey_year" = "summer",
    "precAutumn_survey_year" = "autumn",
    "precWinter_survey_year" = "winter"
  )

### * Calculation survey_year ####
### Prepare data ###
data <- sites %>%
  arrange(id) %>%
  select(
    tempMean_survey_year, tempSpring_survey_year, tempSummer_survey_year,
    tempAutumn_survey_year, tempWinter_survey_year,
    precMean_survey_year, precSpring_survey_year, precSummer_survey_year,
    precAutumn_survey_year, precWinter_survey_year
  )
### Calculate PCA ###
pca <- rda(X = decostand(data, method = "standardize"), scale = TRUE)
summary(pca, axes = 0)
screeplot(pca, bstick = TRUE, npcs = length(pca$CA$eig))
biplot(pca, display = "species", scaling = 2)
### Make data frames ###
eigenvals <- pca %>%
  eigenvals() %>%
  summary() %>%
  as.data.frame() %>%
  rownames_to_column(var = "variables") %>%
  tibble() %>%
  select(PC1:PC3, variables)
data <- pca %>%
  summary() %>%
  magrittr::extract2("species") %>%
  as.data.frame() %>%
  rownames_to_column(var = "variables") %>%
  tibble() %>%
  select(PC1:PC3, variables)
pca_survey_year <- data %>%
  bind_rows(eigenvals)
### create summary table ###
data <- pca %>%
  summary() %>%
  magrittr::extract2("sites") %>%
  as.data.frame() %>%
  rownames_to_column(var = "plot") %>%
  tibble() %>%
  select(plot, PC1:PC3) %>%
  rename(
    pc1_survey_year = PC1,
    pc2_survey_year = PC2,
    pc3_survey_year = PC3
  )
### Add to sites ###
sites <- sites %>%
  left_join(data, by = "plot")

### * Calculation construction_year ####
### Prepare data ###
data <- sites %>%
  arrange(id) %>%
  select(
    plot,
    tempMean_construction_year, tempSpring_construction_year,
    tempSummer_construction_year, tempAutumn_construction_year,
    tempWinter_construction_year,
    tempMean_construction_year_plus, tempSpring_construction_year_plus,
    tempSummer_construction_year_plus, tempAutumn_construction_year_plus,
    tempWinter_construction_year_plus,
    precMean_construction_year, precSpring_construction_year,
    precSummer_construction_year, precAutumn_construction_year,
    precWinter_construction_year,
    precMean_construction_year_plus, precSpring_construction_year_plus,
    precSummer_construction_year_plus, precAutumn_construction_year_plus,
    precWinter_construction_year_plus
  ) %>%
  group_by(plot) %>%
  summarise(across(where(is.numeric), ~ median(.x, na.rm = TRUE))) %>%
  select(-plot)
### Calculate PCA ###
pca <- rda(X = decostand(data, method = "standardize"), scale = TRUE)
summary(pca, axes = 0)
screeplot(pca, bstick = TRUE, npcs = length(pca$CA$eig))
biplot(pca, display = "species", scaling = 2)
### Make data frames ###
eigenvals <- pca %>%
  eigenvals() %>%
  summary() %>%
  as.data.frame() %>%
  rownames_to_column(var = "variables") %>%
  tibble() %>%
  select(PC1:PC3, variables)
data <- pca %>%
  summary() %>%
  magrittr::extract2("species") %>%
  as.data.frame() %>%
  rownames_to_column(var = "variables") %>%
  tibble() %>%
  select(PC1:PC3, variables)
pca_survey_year <- data %>%
  bind_rows(eigenvals)
### create summary table ###
data <- pca %>%
  summary() %>%
  magrittr::extract2("sites") %>%
  as.data.frame() %>%
  rownames_to_column(var = "plot") %>%
  tibble() %>%
  select(plot, PC1:PC3) %>%
  rename(
    pc1_construction_year = PC1,
    pc2_construction_year = PC2,
    pc3_construction_year = PC3
  )
### Add to sites ###
sites <- sites %>%
  left_join(data, by = "plot")


rm(list = setdiff(ls(), c("sites", "sites_splot", "species", "species_splot",
                          "traits", "pca_soil", "pca_construction_year",
                          "pca_survey_year")))


## 7 sPlot open data ##########################################################

data_sites <- sites_splot %>%
  filter(
    #Hay meadow: EUNIS2007 code E2.2 Chytry et al. 2020 Appl Veg Sci
    (ESY == "E22" |
       #Dry grassland: EUNIS2007 code E1.2a Chytry et al. 2020 Appl Veg Sci
       ESY == "E12a") &
      Releve_area >= 10 &
      Releve_area <= 40 &
      Longitude > 10.89845 & # West: Augsburg
      Longitude < 13.46434 & # East: Passau
      Latitude > 	47.85298 & # South: Rosenheim
      Latitude < 49.45095 & # North: Nuernberg
      Elevation < 700
  ) %>%
  rename_with(tolower) %>%
  rename(id = plotobservationid, survey_year = date_of_recording,
         plotSize = releve_area, reference = country) %>%
  mutate(
    id = paste0("X", id),
    reference = str_replace(reference, "Germany", "reference"),
    survey_year = year(survey_year)
  ) %>%
  select(id, survey_year, longitude, latitude, elevation, plotSize,
         reference, esy)


data_species <- species_splot %>%
  rename(id = PlotObservationID, name = Species,
         abundance = Original_abundance) %>%
  mutate(id = paste0("X", id)) %>%
  semi_join(data_sites, by = "id") %>%
  select(id, name, abundance) %>%
  pivot_wider(names_from = "id",
              values_from = "abundance",
              values_fn = sum) %>%
  mutate(
    name = str_replace(name, " ", "_"),
    name = str_replace(name, "Helianthemum_ovatum", "Helianthemum_nummularium"),
    name = str_replace(name, "Galium_album", "Galium_mollugo"),
    name = str_replace(name, "Taraxacum", "Taraxacum_campylodes"),
    name = str_replace(name, "Cerastium_fontanum", "Cerastium_fontanum_ssp_vulgare"),
    name = str_replace(name, "Leucanthemum_ircutianum", "Leucanthemum_vulgare"),
    name = str_replace(name, "Tragopogon_orientalis", "Tragopogon_pratensis"),
    name = factor(name)
  ) %>%
  group_by(name) %>%
  summarise(across(everything(), ~ sum(.x, na.rm = TRUE)))


data_species$name[which(!(data_species$name %in% traits$name))]

species_splot %>%
  select(Species) %>%
  filter(Species == "Galium mollugo")

write_csv(
  data_sites,
  here("data", "processed", "data_processed_sites_splot.csv")
)
write_csv(
  data_species,
  here("data", "processed", "data_processed_species_splot.csv")
)

rm(list = setdiff(ls(), c("sites", "species", "traits", "pca_soil")))


## 8 TBI: Temporal Beta diversity Index ################################

### * Prepare data ####
data_sites <- sites %>%
  ### Choose only plots which were surveyed in each year
  filter(accumulated_cover > 0) %>%
  add_count(plot) %>%
  filter(n == max(n)) %>%
  select(id, plot)
data_species <- species %>%
  select(where(~ !all(is.na(.x)))) %>%
  mutate(across(where(is.numeric), ~ replace(., is.na(.), 0))) %>%
  pivot_longer(-name, names_to = "id", values_to = "value") %>%
  pivot_wider(id, name) %>%
  mutate(
    year = str_match(id, "\\d{4}"),
    plot = str_match(id, "\\d{2}"),
    year = factor(year),
    plot = factor(plot)
  ) %>%
  arrange(id) %>%
  semi_join(data_sites, by = "id") %>%
  select(plot, year, tidyselect::peek_vars(), -id)

### Separate each year in several tibbles ###

for (i in unique(data_species$year)) {
  nam <- paste("species", i, sep = "")

  assign(nam, data_species %>%
    filter(year == i) %>%
    select(-year) %>%
    column_to_rownames(var = "plot"))
}

### a Calculate TBI Presence -------------------------------------------

#### * 2017 vs. 2018 ####
res1718 <- TBI(species2017, species2018,
  method = "sorensen",
  nperm = 9999, test.t.perm = TRUE, clock = TRUE
)
res1718$BCD.summary # B = 0.223, C = 0.155, D = 0.378 (58.9% vs. 41.0%)
res1718$t.test_B.C # p.perm = 0.0058
tbi1718 <- as_tibble(res1718$BCD.mat) %>%
  mutate(comparison = "1718")
#### Test plot
plot(res1718, type = "BC")

#### * 2018 vs. 2019 ####
res1819 <- TBI(species2018, species2019,
  method = "sorensen",
  nperm = 9999, test.t.perm = TRUE, clock = TRUE
)
res1819$BCD.summary # B = 0.118, C = 0.214, D = 0.332 (35.6% vs. 64.3%)
res1819$t.test_B.C # p.perm = 1e-04
tbi1819 <- as_tibble(res1819$BCD.mat) %>%
  mutate(comparison = "1819")
#### Test plot
plot(res1819, type = "BC")

#### * 2019 vs. 2021 ####
res1921 <- TBI(species2019, species2021,
  method = "sorensen",
  nperm = 9999, test.t.perm = TRUE, clock = TRUE
)
res1921$BCD.summary # B = 0.249, C = 0.140, D = 0.390 (63.8% vs. 36.1%)
res1921$t.test_B.C # p.perm = 1e-04
tbi1921 <- as_tibble(res1921$BCD.mat) %>%
  mutate(comparison = "1921")
#### Test plot
plot(res1921, type = "BC")

#### * 2017 vs. 2019 ####
res1719 <- TBI(species2017, species2019,
  method = "sorensen",
  nperm = 9999, test.t.perm = TRUE, clock = TRUE
)
res1719$BCD.summary # B = 0.186, C = 0.212, D = 0.399 (46.7% vs. 53.2%)
res1719$t.test_B.C # p.perm = 0.273
tbi1719 <- as_tibble(res1719$BCD.mat) %>%
  mutate(comparison = "1719")
#### Test plot
plot(res1719, type = "BC")

#### * 2017 vs. 2021 ####
res1721 <- TBI(species2017, species2021,
  method = "sorensen",
  nperm = 9999, test.t.perm = TRUE, clock = TRUE
)
res1721$BCD.summary # B = 0.184, C = 0.450, D = 0.590 (59.0% vs. 40.9%)
res1721$t.test_B.C # p.perm = 0.0021
tbi1721 <- as_tibble(res1721$BCD.mat) %>%
  mutate(comparison = "1721")
#### Test plot
plot(res1721, type = "BC")

#### * Combine datasets ####
data_presence <- bind_rows(tbi1718, tbi1819, tbi1921, tbi1719, tbi1721) %>%
  mutate(presabu = "presence")

### b Calculate TBI Abundance ------------------------------------------

#### * 2017 vs. 2018 ####
res1718 <- TBI(species2017, species2018,
  method = "%diff",
  nperm = 9999, test.t.perm = TRUE, clock = TRUE
)
res1718$BCD.summary # B = 0.213, C = 0.260, D = 0.473 (45.0% vs. 54.9%)
res1718$t.test_B.C # p.perm = 0.1756
tbi1718 <- as_tibble(res1718$BCD.mat) %>%
  mutate(comparison = "1718")
#### Test plot
plot(res1718, type = "BC")

#### * 2018 vs. 2019 ####
res1819 <- TBI(species2018, species2019,
  method = "%diff",
  nperm = 9999, test.t.perm = TRUE, clock = TRUE
)
res1819$BCD.summary # B = 0.167, C = 0.302, D = 0.470 (35.7% vs. 64.2%)
res1819$t.test_B.C # p.perm = 1e-04
tbi1819 <- as_tibble(res1819$BCD.mat) %>%
  mutate(comparison = "1819")
#### Test plot
plot(res1819, type = "BC")

#### * 2019 vs. 2021 ####
res1921 <- TBI(species2019, species2021,
  method = "%diff",
  nperm = 9999, test.t.perm = TRUE, clock = TRUE
)
res1921$BCD.summary # B = 0.331, C = 0.168, D = 0.499 (66.3% vs. 33.6%)
res1921$t.test_B.C # p.perm = 1e-04
tbi1921 <- as_tibble(res1921$BCD.mat) %>%
  mutate(comparison = "1921")
#### Test plot
plot(res1921, type = "BC")

#### * 2017 vs. 2019 ####
res1719 <- TBI(species2017, species2019,
  method = "%diff",
  nperm = 9999, test.t.perm = TRUE, clock = TRUE
)
res1719$BCD.summary # B = 0.210, C = 0.390, D = 0.601 (35.0% vs. 64.9%)
res1719$t.test_B.C # p.perm = 1e-04
tbi1719 <- as_tibble(res1719$BCD.mat) %>%
  mutate(comparison = "1719")
#### Test plot
plot(res1719, type = "BC")

#### * 2017 vs. 2021 ####
res1721 <- TBI(species2017, species2021,
  method = "%diff",
  nperm = 9999, test.t.perm = TRUE, clock = TRUE
)
res1721$BCD.summary # B = 0.301, C = 0.319, D = 0.620 (48.5% vs. 51.4%)
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
### combine abundance and presence data ###
data <- add_row(data_presence, data_abundance) %>%
  mutate(plot = rep(plot, length(data_abundance$comparison) * 2 / 38))
sites_temporal <- sites %>%
  filter(survey_year_factor == "2017") %>%
  left_join(data, by = "plot") %>%
  rename(
    B = "B/(2A+B+C)", C = "C/(2A+B+C)", D = "D=(B+C)/(2A+B+C)",
    change = Change
  ) %>%
  mutate(
    change = C - B,
    plot = as_factor(plot)
  ) %>%
  select(
    plot, block,
    location, locationAbb, locationYear, longitude, latitude,
    riverkm, distanceRiver,
    construction_year,
    exposition, orientation,
    PC1soil, PC2soil, PC3soil, PC1construction_year, PC2construction_year,
    PC3construction_year,
    conf.low, conf.high,
    B, C, D, comparison, presabu
  ) %>%
  mutate(across(
    c(
      PC1soil, PC2soil, PC3soil,
      distanceRiver,
      B, C, D
    ),
    ~ round(.x, digits = 4)
  ))

rm(list = setdiff(ls(), c(
  "sites", "species", "traits", "sites_temporal",
  "pca_construction_year", "pca_soil", "pca_survey_year"
)))


## 9 Finalization ######################################################

### a Rounding ---------------------------------------------------------

sites <- sites %>%
  mutate(across(
    c(
      lcbd,
      syn_total, syn_trend, syn_detrend,
      pc1_soil, pc2_soil, pc3_soil, pc1_survey_year, pc2_survey_year,
      pc3_survey_year,
      pc1_construction_year, pc2_construction_year, pc3_construction_year
    ),
    ~ round(.x, digits = 4)
  )) %>%
  mutate(across(
    c(
      n_total_rerc, target_cover_ratio, graminoid_cover_ratio,
      target_richness_ratio, shannon, eveness
    ),
    ~ round(.x, digits = 3)
  )) %>%
  mutate(across(
    c(river_distance, accumulated_cover),
    ~ round(.x, digits = 1)
  ))

### b Final selection of variables -------------------------------------

sites <- sites %>%
  select(
    id, plot, block, botanist,
    # space
    location, location_abb, location_construction_year, latitude, longitude,
    river_km, river_distance,
    mem1_2017, mem2_2018, mem1_2019, mem1_2021, mem2_2021,
    # time
    survey_year, construction_year, plot_age,
    # local site characteristics
    exposition, orientation, pc1_soil, pc2_soil, pc3_soil,
    # historical factors
    pc1_construction_year, pc2_construction_year, pc3_construction_year,
    # response variables
    accumulated_cover, species_richness, eveness, shannon,
    graminoid_cover_ratio, ruderal_cover,
    target_richness, target_richness_ratio, target_cover_ratio
  )

sites_temporal <- sites_temporal %>%
  pivot_wider(names_from = "presabu", values_from = c("B", "C", "D")) %>%
  select(
    plot, block, comparison,
    # space
    location, locationAbb, locationYear, latitude, longitude, riverkm,
    distanceRiver,
    # time
    construction_year,
    # local site characteristics
    exposition, orientation, PC1soil, PC2soil, PC3soil,
    # historical factors
    PC1construction_year, PC2construction_year, PC3construction_year,
    # temporal beta-diversity indices (TBI)
    B_presence, B_abundance, C_presence, C_abundance,
    D_presence, D_abundance,
    # other variables
    conf.low, conf.high
  )

sites_restoration <- sites %>%
  select(
    id, plot, block,
    # space
    location, locationAbb, locationYear, latitude, longitude, riverkm,
    distanceRiver,
    # time
    survey_year, construction_year, plotAge,
    # local site characteristics
    exposition, orientation, PC1soil, PC2soil, PC3soil,
    # response variables
    accumulated_cover, richness, eveness, graminoid_cover_ratio,
    ruderal_cover,
    # conservation
    target_richness, targetRich_ratio, rlg_richness, target_cover_ratio,
    # legal evaluation
    biotopeType, ffh, changeType, baykompv, biotopePoints, min8, min9,
    botanist, conf.low, conf.high
  )

### c Final selection of plots -----------------------------------------

sites_spatial <- sites %>%
  ### Choose only plots which were surveyed in each year
  filter(accumulated_cover > 0) %>%
  add_count(plot) %>%
  filter(n == max(n)) %>%
  select(-n) %>%
  mutate(plot = factor(plot)) # check number of plots

sites_temporal <- sites_temporal %>%
  semi_join(sites_spatial, by = "plot") %>%
  mutate(plot = factor(plot)) # check number of plots



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# C Save processed data ###############################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


### Data ###
write_csv(
  sites_spatial,
  here("data", "processed", "data_processed_sites_spatial.csv")
)
write_csv(
  species,
  here("data", "processed", "data_processed_species.csv")
)
write_csv(
  traits,
  here("data", "processed", "data_processed_traits.csv")
)
write_csv(
  sites_temporal,
  here("data", "processed", "data_processed_sites_temporal.csv")
)

### Tables ###
write_csv(
  pca_soil,
  here("outputs", "statistics", "pca_soil.csv")
)
write_csv(
  pca_survey_year,
  here("outputs", "statistics", "pca_survey_year.csv")
)
write_csv(
  pca_construction_year,
  here("outputs", "statistics", "pca_construction_year.csv")
)
