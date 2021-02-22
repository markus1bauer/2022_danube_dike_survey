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



### 1 Load data #####################################################################################

### a Sites -------------------------------------------------------------------------------------------
sites <- read_csv2("data_raw_sites_2017_2018_2019.csv", col_names = T, na = "na", col_types = 
                     cols(
                       .default = col_double(),
                       id = col_factor(),
                       location = col_factor(),
                       ageCategory = col_factor(),
                       HCl = col_factor(),
                       phosphorousClass = col_factor(),
                       potassiumClass = col_factor(),
                       magnesiumClass = col_character()
                     )        
)
sites <- select(sites, id, location, RW, HW, constructionYear, pH, Corg, organic, nitrogen, cnRatio, CaCO3, sand, silt, clay, phosphorous, potassium, magnesium)

### b Species -------------------------------------------------------------------------------------------
species <- read_csv("data_raw_species_2017_2018_2019.csv", col_names = F, na = "na")
names <- species %>% 
  slice(1) %>%
  pivot_longer(-(X1:X4), names_to = "id", values_to = "name") %>%
  select(-(X1:X4))
names$name <- str_replace(names$name, " ", "_")
any(duplicated(names$name))
any(duplicated(names$id))
species <- species %>%
  select(-(X1:X3)) %>%
  slice(-c(1:2, 199)) %>%
  pivot_longer(-(X4), names_to = "id", values_to = "abu") %>%
  rename(plot = X4) %>%
  pivot_wider(names_from = plot, names_prefix = "X", values_from = abu)
species <- left_join(names, species, by = "id")
species$name <- as_factor(species$name)
species <- species %>%
  select(-id) %>%
  mutate_if(is.character, as.numeric)
rm(names)

### Check that each species occurs at least one time ###
species <- species %>%
  rowwise() %>%
  mutate(total = sum(c_across(X01_m_2017:X70_m_2019))) %>%
  mutate(presence = if_else(total > 0, 1, 0)) %>%
  filter(presence == 1) %>%
  ungroup() %>%
  select(-total, -presence)
  

    

### c Traits -------------------------------------------------------------------------------------------
#traits <- read_csv("data_raw_traits.csv", col_names = T, na = "na", col_types = 
                      cols(
                        .default = col_double(),
                        name = col_factor(),
                        descriptor = col_factor(),
                        nomenclature = col_factor(),
                        subspecies = col_factor(),
                        germanName = col_factor(),
                        abb = col_factor(),
                        family = col_factor(),
                        iucn = col_factor(),
                        responsibility = col_factor(),
                        legal = col_factor(),
                        rlg = col_factor(),
                        rlb = col_factor(),
                        target = col_factor()
                      )
                   )

traits <- left_join(species, traits, by = "name")
traits <- traits %>%
  select(name, abb, family, rlg, rlb, n, target)

n_miss(sites)
n_miss(species)
n_miss(traits)
miss_var_summary(traits)
vis_miss(traits, cluster = F)
gg_miss_var(traits)
gg_miss_case(traits, order_cases = F)
gg_miss_upset(traits)


### 2 Create variables #####################################################################################

### a Dummies for confidence interval -------------------------------------------------------------------------------------------
sites$conf.low <- c(1:70)
sites$conf.high <- c(1:70)

### b Graminoid's cover ratio -------------------------------------------------------------------------------------------
graminoidsCovratio <- species
graminoidsCovratio$family <- traits$family
graminoidsCovratio <- graminoidsCovratio %>%
  gather("ID", "n", X84I1:X18S18) %>%
  mutate(plot = ID) %>%
  separate(plot, c("year","plot"), sep = 3) %>%
  mutate(year = if_else(year == "X84", 1984, if_else(year == "X93", 1993, if_else(year == "X03", 2003, 2018)))) %>%
  group_by(year, ID, family) %>%
  summarise(sum = sum(n)) %>%
  mutate(herbType = if_else(family == "Poaceae" | family=="Cyperaceae", "graminoids", "other")) %>%
  group_by(year, ID, herbType) %>%
  summarise(sum = sum(sum)) %>%
  spread(herbType, sum) %>%
  transmute(graminoidCovratio = graminoids / (graminoids + other)) %>%
  ungroup() %>%
  select(-year)
sites <- right_join(graminoidsCovratio, sites, by = "ID")
rm(graminoidsCovratio)

### c Species richness -------------------------------------------------------------------------------------------
specRich <- species
specRich$target <- traits$target
specRich$rlg <- traits$rlg
specRich$rlb <- traits$rlb
### Total species richness ###
specRichAll <- specRich %>%
  pivot_longer(names_to = "id", values_to = "n", cols = X01_m_2017:X70_m_2019) %>%
  group_by(id) %>%
  mutate(n = if_else(n > 0, 1, 0)) %>%
  summarise(total = sum(n)) %>%
  group_by(year, id) %>%
  summarise(speciesRichnness = sum(total)) %>%
  mutate_if(is.character, as.factor) %>%
  ungroup()
### Species richness of target species ###
specRichTarget <- specRich %>%
  pivot_longer(names_to = "id", values_to = "n", cols = X01_m_2017:X70_m_2019) %>%
  group_by(id, taarget) %>%
  mutate(n = if_else(n > 0, 1, 0)) %>%
  summarise(total = sum(n)) %>%
  mutate(target = if_else(target == "T" | target == "Te", "target", "non-target")) %>%
  filter(target == "target") %>%
  group_by(year, id) %>%
  summarise(target = sum(total)) %>%
  mutate_if(is.character, as.factor) %>%
  ungroup()
### species richness of species of the red list Germany ###
specRichRLG <- specRich %>%
  pivot_longer(names_to = "id", values_to = "n", cols = X01_m_2017:X70_m_2019) %>%
  group_by(id, target) %>%
  mutate(n = if_else(n > 0, 1, 0)) %>%
  summarise(total = sum(n)) %>%
  filter(rlg == "1" | rlg == "2" | rlg == "3" | rlg == "V") %>%
  group_by(year, id) %>%
  summarise(rlg = sum(total)) %>%
  mutate_if(is.character, as.factor) %>%
  ungroup()
### species richness of species of the red list Bavaria ###
specRichRLB <- specRich %>%
  pivot_longer(names_to = "id", values_to = "n", cols = X01_m_2017:X70_m_2019) %>%
  group_by(id, rlb) %>%
  mutate(n = if_else(n > 0, 1, 0)) %>%
  summarise(total = sum(n)) %>%
  filter(rlb == "1" | rlb == "2" | rlb == "3" | rlb == "V") %>%
  group_by(year, id) %>%
  summarise(rlb = sum(total)) %>%
  mutate_if(is.character, as.factor) %>%
  ungroup()
### implement in sites data set ###
sites <- right_join(specRichAll, sites, by = "id")
sites <- right_join(specRichTarget, sites, by = "id")
sites <- right_join(specRichRLG, sites, by = "id")
sites <- right_join(specRichRLB, sites, by = "id")
sites$tarRichratio <- sites$target / sites$all
rm(specRich, specRichAll, specRichTarget, specRichRLG, specRichRLB)


### d CWM of Ellenberg N -------------------------------------------------------------------------------------------
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
#Npresabs <- dbFD(Ntraits, Nspecies, w.abun = F,
#                        calc.FRic = F, calc.FDiv = F, corr = "sqrt")
### implement in sites data set ###
sites$cwmAbuN <- as.numeric(as.character(Nweighted$CWM$n))
#sites$cwmPresN <- as.numeric(as.character(Npresabs$CWM$n))
rm(Ntraits, Nspecies, Nweighted)

### e Import functional plant traits -------------------------------------------------------------------------------------------
#### * make names equal #### 
dataSLA <- read_table2("data_raw_SLA.txt", col_names = T, na = "na", col_types = 
                         cols(
                           .default = col_double(),
                           sbs_name = col_factor()
                           )
                       )
dataSM <- read_table2("data_raw_seedmass.txt", col_names = T, na = "na", col_types = 
                         cols(
                           .default = col_double(),
                           sbs_name = col_factor()
                           )
                      )
dataH <- read_table2("data_raw_canopy_height.txt", col_names = T, na = "na", col_types = 
                         cols(
                           .default = col_double(),
                           sbs_name = col_factor()
                           )
                     )
n_miss(dataSLA); pct_miss(dataSLA)
n_miss(dataSM); pct_miss(dataSM)
n_miss(dataH); pct_miss(dataH)
dataH <- drop_na(dataH)

#traits[which(!(traits$name %in% dataSLA$sbs_name)),]
#traits[which(!(traits$name %in% dataSM$sbs_name)),]
#traits[which(!(traits$name %in% dataH$sbs_name)),]
##### Synonyme ###
dataSLA$sbs_name <- fct_recode(dataSLA$sbs_name, Achillea_seidlii = "Achillea_pannonica")
dataSM$sbs_name <- fct_recode(dataSM$sbs_name, Achillea_seidlii = "Achillea_pannonica")
dataH$sbs_name <- fct_recode(dataH$sbs_name, Achillea_seidlii = "Achillea_pannonica")
dataSLA$sbs_name <- fct_recode(dataSLA$sbs_name, Anemone_pulsatilla = "Pulsatilla_vulgaris")
dataSM$sbs_name <- fct_recode(dataSM$sbs_name, Anemone_pulsatilla = "Pulsatilla_vulgaris")
dataH$sbs_name <- fct_recode(dataH$sbs_name, Anemone_pulsatilla = "Pulsatilla_vulgaris")
dataSLA$sbs_name <- fct_recode(dataSLA$sbs_name, Erica_herbacea = "Erica_cinerea")
dataSM$sbs_name <- fct_recode(dataSM$sbs_name,  Erica_herbacea = "Erica_cinerea")
dataH$sbs_name <- fct_recode(dataH$sbs_name,  Erica_herbacea = "Erica_cinerea")
dataSLA$sbs_name <- fct_recode(dataSLA$sbs_name, Euphorbia_verrucosa = "Euphorbia_brittingeri")
dataSM$sbs_name <- fct_recode(dataSM$sbs_name,  Euphorbia_verrucosa = "Euphorbia_brittingeri")
dataH$sbs_name <- fct_recode(dataH$sbs_name,  Euphorbia_verrucosa = "Euphorbia_brittingeri")
dataSLA$sbs_name <- fct_recode(dataSLA$sbs_name, Helictotrichon_pratense  = "Avenula_pratensis")
dataSM$sbs_name <- fct_recode(dataSM$sbs_name,  Helictotrichon_pratense = "Avenula_pratensis")
dataH$sbs_name <- fct_recode(dataH$sbs_name,  Helictotrichon_pratense = "Avenula_pratensis")
dataSLA$sbs_name <- fct_recode(dataSLA$sbs_name, Helictotrichon_pubescens = "Avenula_pubescens")
dataSM$sbs_name <- fct_recode(dataSM$sbs_name, Helictotrichon_pubescens = "Avenula_pubescens")
dataH$sbs_name <- fct_recode(dataH$sbs_name, Helictotrichon_pubescens = "Avenula_pubescens")
dataSLA$sbs_name <- fct_recode(dataSLA$sbs_name, Medicago_falcata = "Medicago_sativa_s._falcata")
dataSM$sbs_name <- fct_recode(dataSM$sbs_name, Medicago_falcata = "Medicago_sativa_s._falcata")
dataH$sbs_name <- fct_recode(dataH$sbs_name, Medicago_falcata = "Medicago_sativa_s._falcata")
dataSLA$sbs_name <- fct_recode(dataSLA$sbs_name, Pilosella_officinarum = "Hieracium_pilosella")
dataSM$sbs_name <- fct_recode(dataSM$sbs_name, Pilosella_officinarum = "Hieracium_pilosella")
dataH$sbs_name <- fct_recode(dataH$sbs_name, Pilosella_officinarum = "Hieracium_pilosella")
dataSLA$sbs_name <- fct_recode(dataSLA$sbs_name, Securigera_varia = "Coronilla_varia")
dataSM$sbs_name <- fct_recode(dataSM$sbs_name, Securigera_varia = "Coronilla_varia")
dataH$sbs_name <- fct_recode(dataH$sbs_name, Securigera_varia = "Coronilla_varia")
dataSLA$sbs_name <- fct_recode(dataSLA$sbs_name, Taraxacum_campylodes = "Taraxacum_officinale")
dataSM$sbs_name <- fct_recode(dataSM$sbs_name, Taraxacum_campylodes = "Taraxacum_officinale")
dataH$sbs_name <- fct_recode(dataH$sbs_name, Taraxacum_campylodes = "Taraxacum_officinale")
dataSLA$sbs_name <- fct_recode(dataSLA$sbs_name, Taraxacum_campylodes = "Taraxacum_Sec._Ruderalia")
dataSM$sbs_name <- fct_recode(dataSM$sbs_name, Taraxacum_campylodes = "Taraxacum_Sec._Ruderalia")
dataH$sbs_name <- fct_recode(dataH$sbs_name, Taraxacum_campylodes = "Taraxacum_Sec._Ruderalia")
dataSLA$sbs_name <- fct_recode(dataSLA$sbs_name, Taraxacum_campylodes = "Taraxacum_species")
dataSM$sbs_name <- fct_recode(dataSM$sbs_name, Taraxacum_campylodes = "Taraxacum_species")
dataH$sbs_name <- fct_recode(dataH$sbs_name, Taraxacum_campylodes = "Taraxacum_species")
##### Take value of similar species ###
Fspecies <- species
traits$name <- fct_recode(traits$name, Carex_montana = "Carex_ericetorum")
Fspecies$name <- fct_recode(Fspecies$name, Carex_montana = "Carex_ericetorum")
traits$name <- fct_recode(traits$name, Centaurea_jacea = "Centaurea_pannonica")
Fspecies$name <- fct_recode(Fspecies$name, Centaurea_jacea = "Centaurea_pannonica")
traits$name <- fct_recode(traits$name, Cytisus_scoparius = "Cytisus_ratisbonensis")
Fspecies$name <- fct_recode(Fspecies$name, Cytisus_scoparius = "Cytisus_ratisbonensis")
traits$name <- fct_recode(traits$name, Hieracium_pilosella = "Pilosella_macrantha")
Fspecies$name <- fct_recode(Fspecies$name, Hieracium_pilosella = "Pilosella_macrantha")
traits$name <- fct_recode(traits$name, Pulsatilla_vulgaris = "Pulsatilla_grandis")
Fspecies$name <- fct_recode(Fspecies$name, Pulsatilla_vulgaris = "Pulsatilla_grandis")
##### Subspecies ###
traits$name <- fct_recode(traits$name, Euphrasia_officinalis = "Euphrasia_picta")
Fspecies$name <- fct_recode(Fspecies$name, Euphrasia_officinalis = "Euphrasia_picta")
traits$name <- fct_recode(traits$name, Helianthemum_nummularium = "Helianthemum_ovatum")
Fspecies$name <- fct_recode(Fspecies$name, Helianthemum_nummularium = "Helianthemum_ovatum")
traits$name <- fct_recode(traits$name, Lotus_corniculatus = "Lotus_corniculatus_ssp_corniculatus")
Fspecies$name <- fct_recode(Fspecies$name, Lotus_corniculatus = "Lotus_corniculatus_ssp_corniculatus")
traits$name <- fct_recode(traits$name, Lotus_corniculatus = "Lotus_corniculatus_ssp_hirsutus")
Fspecies$name <- fct_recode(Fspecies$name, Lotus_corniculatus = "Lotus_corniculatus_ssp_hirsutus")
traits$name <- fct_recode(traits$name, Potentilla_cinerea = "Potentilla_incana")
Fspecies$name <- fct_recode(Fspecies$name, Potentilla_cinerea = "Potentilla_incana")

#### * missing values: SLA #### 
####SLA von Cerabolini et al. (2010)
dataSLA <- dataSLA %>% add_row(sbs_name = 'Biscutella_laevigata', value = 18.5)
dataSLA <- dataSLA %>% add_row(sbs_name = "Dorycnium_pentaphyllum_ssp_germanicum", value = 20.8) #Bei cerabolini als D. pentaphyllum gef?hrt
dataSLA <- dataSLA %>% add_row(sbs_name = "Globularia_cordifolia", value = 6.3)
#####SLA von Pierce et al. (2007)
dataSLA <- dataSLA %>% add_row(sbs_name = "Gentiana_verna", value = 16.3)
#####SLA von Pipenbaher et al. (2007)
dataSLA <- dataSLA %>% add_row(sbs_name = "Biscutella_laevigata", value = 14.5)
dataSLA <- dataSLA %>% add_row(sbs_name = "Dorycnium_pentaphyllum_ssp_germanicum", value = 18.2) #Bei Pipenbaher als D. bavaricum gef?hrt
dataSLA <- dataSLA %>% add_row(sbs_name = "Euphorbia_brittingeri", value = 16.2) #Bei Pipenbaher als E. verrucosa gef?hrt
dataSLA <- dataSLA %>% add_row(sbs_name = "Festuca_rupicola", value = 17.7)
dataSLA <- dataSLA %>% add_row(sbs_name = "Globularia_cordifolia", value = 7.0)
dataSLA <- dataSLA %>% add_row(sbs_name = "Helianthemum_nummularium", value = 7.0) #Bei Pipenbaher als H. ovatum gef?hrt
dataSLA <- dataSLA %>% add_row(sbs_name = "Rhinanthus_aristatus", value = 15.9)

#### * missing values: seed mass #### 
#which(!(Fspecies849318$name %in% dataSM$sbs_name));Fspecies849318$name[which(!(Fspecies849318$name %in% dataSM$sbs_name))]
#which(!(Fspecies0318$name %in% dataSM$sbs_name));Fspecies0318$name[which(!(Fspecies0318$name %in% dataSM$sbs_name))]
######Seedmass-Daten: Eigenannahme von Conradi & Kollmann (2016)
dataSM <- dataSM %>% add_row(sbs_name = "Platanthera_bifolia", value = 0.001)
#####Seedmass-Daten: Hintze et al. (2013)
dataSM <- dataSM %>% add_row(sbs_name = "Cytisus_ratisbonensis", value = 5.9) # Bei Hintze als Pulsatilla patens gef?hrt
dataSM <- dataSM %>% add_row(sbs_name = "Dorycnium_pentaphyllum_ssp_germanicum", value = 2.45) # Bei Hintze als D. germanicum gef?hrt
dataSM <- dataSM %>% add_row(sbs_name = "Euphorbia_brittingeri", value = 2.635) # Bei Hintze als E. verrucosa gef?hrt
dataSM <- dataSM %>% add_row(sbs_name = "Potentilla_cinerea", value = 0.131) # Bei Hintze als P. incana gef?hrt
dataSM <- dataSM %>% add_row(sbs_name = "Seseli_annuum", value = 0.5487)
dataSM <- dataSM %>% add_row(sbs_name = "Viola_rupestris", value = 1.12)
dataSM <- dataSM %>% add_row(sbs_name = "Peucedanum_oreoselinum", value = 27.91) # mean(c(36.3,26,21.426))=27.91; niedrige Werte von P. odoratum in LEDA-Datenbank

#### * missing values: canopy height #### 
#which(!(Fspecies849318$name %in% dataH$sbs_name));Fspecies849318$name[which(!(Fspecies849318$name %in% dataH$sbs_name))]
#which(!(Fspecies0318$name %in% dataH$sbs_name));Fspecies0318$name[which(!(Fspecies0318$name %in% dataH$sbs_name))]
#####Height von Pipenbaher et al. (2007)
dataH <- dataH %>% add_row(sbs_name = 'Brachypodium_pinnatum', value = 0.880)
dataH <- dataH %>% add_row(sbs_name = "Dorycnium_pentaphyllum_ssp_germanicum", value = 0.211) #Bei Pipenbaher als D. bavaricum gef?hrt
dataH <- dataH %>% add_row(sbs_name = "Rhinanthus_aristatus", value = 0.190)
#####Height von J?ger (2011)
dataH <- dataH %>% add_row(sbs_name = "Anemone_patens", value = 0.21) # Bei J?ger als Pulsatilla patens gef?hrt
dataH <- dataH %>% add_row(sbs_name = "Cytisus_ratisbonensis", value = 0.35)
dataH <- dataH %>% add_row(sbs_name = "Seseli_annuum", value = 0.5)
dataH <- dataH %>% add_row(sbs_name = "Taraxacum_campylodes", value = 0.28)

#### * put values to traits data frame #### 
dataSLA$sbs_name <- factor(dataSLA$sbs_name)
dataSM$sbs_name <- factor(dataSM$sbs_name)
dataH$sbs_name <- factor(dataH$sbs_name)
dataSLA <- dataSLA %>%
  rename(name = sbs_name, sla = value) %>%
  drop_na() %>%
  group_by(name) %>%
  summarise(sla = median(sla))
dataSM <- dataSM %>%
  rename(name = sbs_name, seedmass = value) %>%
  drop_na() %>%
  group_by(name) %>%
  summarise(seedmass = median(seedmass))
dataH <- dataH %>%
  rename(name = sbs_name, height = value) %>%
  drop_na() %>%
  group_by(name) %>%
  summarise(height = median(height))
traits <- left_join(traits, dataSLA, by = "name")
traits <- left_join(traits, dataSM, by = "name")
traits <- left_join(traits, dataH, by = "name")
### * Check completeness ####
n_miss(traits$sla)
n_miss(traits$seedmass)
n_miss(traits$height)
pct_complete(traits$sla)
pct_complete(traits$seedmass)
pct_complete(traits$height)
vis_miss(traits, cluster = F, sort_miss = T)
gg_miss_var(traits)
gg_miss_case(traits, order_cases = F)

### * Prepare data frames ####
Fspecies <- Fspecies %>%
  group_by(name) %>%
  mutate(sum = sum(c_across(X84I1:X18III16))) %>%
  mutate(presence = if_else(sum > 0, 1, 0)) %>%
  filter(presence == 1) %>%
  ungroup() %>%
  select(-sum, -presence) %>%
  group_by(name) %>%
  summarise(across(starts_with("X"), list(median = median)))
Ftraits <- traits %>%
  select(name, sla, height, seedmass) %>%
  filter(sla > 0 | height > 0 | seedmass > 0) %>%
  group_by(name) %>%
  summarise(sla = median(sla), 
            seedmass = median(seedmass),
            height = median(height)
  )
FtraitsAll <- Ftraits %>%
  select(name, sla, seedmass, height) %>%
  drop_na()
FtraitsSla <- Ftraits %>%
  select(name, sla) %>%
  drop_na()  
FtraitsSeedmass <- Ftraits %>%
  select(name, seedmass) %>%
  drop_na()
FtraitsHeight <- Ftraits %>%
  select(name, height) %>%
  drop_na()
rm(dataSLA,dataSM,dataH,Ftraits)

### f CWM and FDis of functional plant traits -------------------------------------------------------------------------------------------
### * All ####
Tspecies <- semi_join(Fspecies, FtraitsAll, by = "name")
Ttraits <- semi_join(FtraitsAll, Tspecies, by = "name")
Tspecies <- Tspecies %>%
  pivot_longer(-name, "site", "value") %>%
  pivot_wider(site, name) %>%
  column_to_rownames("site")
Ttraits <- column_to_rownames(Ttraits, "name")
log_Ttraits <- log(Ttraits + 1)
TdiversityAbu <- dbFD(log_Ttraits, Tspecies, w.abun = T,
               calc.FRic = F, calc.FDiv = F, corr = "cailliez")
sites$fdisAbuAll <- TdiversityAbu$FDis
### * SLA ####
Tspecies <- semi_join(Fspecies, FtraitsSla, by = "name")
Ttraits <- semi_join(FtraitsSla, Tspecies, by = "name")
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
### * Seed mass ####
Tspecies <- semi_join(Fspecies, FtraitsSeedmass, by = "name")
Ttraits <- semi_join(FtraitsSeedmass, Tspecies, by = "name")
Tspecies <- Tspecies %>%
  pivot_longer(-name, "site", "value") %>%
  pivot_wider(site, name) %>%
  column_to_rownames("site")
Ttraits <- column_to_rownames(Ttraits, "name")
log_Ttraits <- log(Ttraits + 1)
TdiversityAbu <- dbFD(log_Ttraits, Tspecies, w.abun = T, 
                      calc.FRic = F, calc.FDiv = F, corr = "sqrt");
sites$fdisAbuSeedmass <- TdiversityAbu$FDis
sites$cwmAbuSeedmass <- exp(as.numeric(as.character(TdiversityAbu$CWM$seedmass)))
### * Height ####
Tspecies <- semi_join(Fspecies, FtraitsHeight, by = "name")
Ttraits <- semi_join(FtraitsHeight, Tspecies, by = "name")
Tspecies <- Tspecies %>%
  pivot_longer(-name, "site", "value") %>%
  pivot_wider(site, name) %>%
  column_to_rownames("site")
Ttraits <- column_to_rownames(Ttraits, "name")
log_Ttraits <- log(Ttraits + 1)
TdiversityAbu <- dbFD(log_Ttraits, Tspecies, w.abun = T, 
                      calc.FRic = F, calc.FDiv = F, corr = "sqrt");
sites$fdisAbuHeight <- TdiversityAbu$FDis
sites$cwmAbuHeight <- exp(as.numeric(as.character(TdiversityAbu$CWM$height)))
rm(TdiversityAbu,TdiversityPres,Tspecies,Ttraits,Fspecies,FtraitsAll,FtraitsSla,FtraitsSeedmass,FtraitsHeight)


### 3 Separate data sets into 849318 and 0318 #####################################################################################

sites849318 <- filter(sites, dataset == "transects")
sites0318 <- filter(sites, dataset == "blocks")
species849318 <- species %>%
  select(name, starts_with("X84"), starts_with("X93"), starts_with("X18I")) %>%
  group_by(name) %>%
  mutate(sum = sum(c_across(X84I1:X18III16))) %>%
  mutate(presence = if_else(sum > 0, 1, 0)) %>%
  filter(presence == 1) %>%
  ungroup() %>%
  select(-sum, -presence)
species0318 <- species %>%
  select(name, starts_with("X03"), starts_with("X18N"), starts_with("X18M"), starts_with("X18S")) %>%
  group_by(name) %>%
  mutate(sum = sum(c_across(X03N01:X18S18))) %>%
  mutate(presence = if_else(sum > 0, 1, 0)) %>%
  filter(presence == 1) %>%
  ungroup() %>%
  select(-sum, -presence)


### 4 Save processed data #####################################################################################

setwd("Z:/Documents/0_Uni/2018_Projekt_9_Masterthesis/3_Aufnahmen_und_Ergebnisse/2020_monitoring_Garchinger_Heide/data/processed")
write_csv2(sites849318, "data_processed_sites849318.csv")
write_csv2(sites0318, "data_processed_sites0318.csv")
write_csv2(species849318, "data_processed_species849318.csv")
write_csv2(species0318, "data_processed_species0318.csv")
write_csv2(traits, "data_processed_traits.csv")

