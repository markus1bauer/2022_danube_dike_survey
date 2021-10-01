
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
rm(list = setdiff(ls(), c("sites", "species", "traits")))


### 6 Load functional plant traits #####################################################################################

#### * read LEDA data #### 
dataSLA <- data.table::fread("data_raw_traitbase_leda_20210223_sla.txt", 
                             sep = ";",
                             dec = ".",
                             skip = 3,
                             header = T,
                             select = c("SBS name", "single value [mm^2/mg]")
) %>%
  as_tibble() %>%
  rename(name = "SBS name") %>%
  rename(sla = "single value [mm^2/mg]") %>%
  mutate(name = as_factor(str_replace_all(name, " ", "_")))

dataSM <- data.table::fread("data_raw_traitbase_leda_20210223_seedmass.txt", 
                            sep = ";",
                            dec = ".",
                            skip = 3,
                            header = T,
                            select = c("SBS name", "single value [mg]")
) %>%
  as_tibble() %>%
  rename(name = "SBS name") %>%
  rename(seedmass = "single value [mg]") %>%
  mutate(name = as_factor(str_replace_all(name, " ", "_")))

dataH <- data.table::fread("data_raw_traitbase_leda_20210223_canopy_height.txt", 
                           sep = ";",
                           dec = ".",
                           skip = 3,
                           header = T,
                           select = c("SBS name", "single value [m]")
) %>%
  as_tibble() %>%
  rename(name = "SBS name",
         height = "single value [m]") %>%
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
  mutate(name = fct_recode(name, Vicia_sativa_ssp_nigra = "Vicia_sativa_s._nigra")) %>%  
  group_by(name) %>%
  summarise(across(where(is.double), ~median(.x, na.rm = T))) %>%
  right_join(traits, by = "name")
### check completeness of LEDA ###
#(test2 <- traits %>%
#select(name, sla, seedmass, height) %>%
#filter(!complete.cases(.)))

### * read TRY data ####
data <- data.table::fread("data_raw_traitbase_try_20210306_13996.txt", 
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
  mutate(sla = coalesce(sla.x, sla.y),
         seedmass = coalesce(seedmass.x, seedmass.y),
         height = coalesce(height.x, height.y),
         .keep = "unused")
### check completeness of LEDA + TRY ###
#(test2 <- traits %>%
#select(name, sla, seedmass, height) %>%
#filter(!complete.cases(.)))

### * read GrooT data ####
data <- read.csv("data_raw_traitbase_groot_root_traits.csv", header = T, na.strings = c("", "NA")) %>%
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

### * read BiolFlor data ####
data <- full_join(
  read.csv2("data_raw_traitbase_biolflor_flowering_timing.csv", header = T, na.strings = c("", "NA")),
  read.csv2("data_raw_traitbase_biolflor_species_names.csv", header = T, na.strings = c("", "NA")),
  by = "ARTID"
) %>%
  rename_all(.funs = tolower) %>%
  rename(name = species, flowerBegin = beginn, flowerEnd = ende, flowerDuration = dauer) %>%
  select(name, flowerBegin, flowerEnd, flowerDuration) %>%
  mutate(name = str_replace_all(name, " ", "_"))
##### Find Synonyms ###
#traits$name[which(!(traits$name %in% data$name))]
#data %>%
#filter(str_detect(name, "Chamomilla"))
traits <- data %>%
  mutate(name = fct_recode(name, Carex_praecox_ssp_curvata = "Carex_curvata")) %>%
  mutate(name = fct_recode(name, Carex_praecox_ssp_praecox = "Carex_praecox")) %>%
  mutate(name = fct_recode(name, Cerastium_fontanum_ssp_vulgare = "Cerastium_fontanum")) %>%
  mutate(name = fct_recode(name, Cota_tinctoria = "Anthemis_tinctoria")) %>%
  mutate(name = fct_recode(name, Cyanus_segetum = "Centaurea_cyanus")) %>%
  mutate(name = fct_recode(name, Jacobaea_vulgaris = "Senecio_jacobaea")) %>%
  mutate(name = fct_recode(name, Ononis_spinosa_ssp_procurrens = "Ononis_repens")) %>%
  mutate(name = fct_recode(name, Plantago_major_ssp_major = "Plantago_major")) %>%
  mutate(name = fct_recode(name, Plantago_major_ssp_intermedia = "Plantago_intermedia")) %>%
  mutate(name = fct_recode(name, Pilosella_caespitosa = "Hieracium_caespitosum")) %>%
  mutate(name = fct_recode(name, Pilosella_officinarum = "Hieracium_pilosella")) %>%
  mutate(name = fct_recode(name, Pilosella_piloselloides = "Hieracium_piloselloides")) %>%
  mutate(name = fct_recode(name, Ranunculus_serpens_ssp_nemorosus = "Ranunculus_nemorosus")) %>%
  mutate(name = fct_recode(name, "Silene_flos-cuculi" = "Lychnis_flos-cuculi")) %>%
  mutate(name = fct_recode(name, Silene_latifolia_ssp_alba = "Silene_latifolia")) %>%
  mutate(name = fct_recode(name, Stachys_officinalis = "Betonica_officinalis")) %>%
  mutate(name = fct_recode(name, Taraxacum_campylodes = "Taraxacum_sect._Ruderalia_")) %>%
  mutate(name = fct_recode(name, Vicia_sativa_ssp_nigra = "Vicia_sativa")) %>%
  right_join(traits, by = "name")
### check completeness of Roots ###
(test2 <- traits %>%
    select(name, flowerBegin, flowerEnd, flowerDuration) %>%
    filter(!complete.cases(flowerBegin, flowerEnd, flowerDuration)))

### * read CLOPLA data ####
data <- data.table::fread("data_raw_traitbase_clopla _20210608_clonal_traits.csv", 
                          header = T, 
                          sep = ";", 
                          dec = ".", 
                          quote = "") %>%
  as_tibble() #%>%
data %>%
  select(Species, CGO_Type_1, OffspringPerParent_1, LateralSpread_1) %>%#, BudBank_no_>10cm, BudBank_no_>0_to_10cm, BudBank_no_0cm, BudBank_no_0_to_-10cm, BudBank_no_<-10cm) %>%
  mutate(LateralSpread_1 = str_replace_all(LateralSpread_1, "<0.01", "1")) %>%
  mutate(LateralSpread_1 = str_replace_all(LateralSpread_1, "0.01-0.25", "2")) %>%
  mutate(LateralSpread_1 = str_replace_all(LateralSpread_1, ">0.25", "3")) %>%
  mutate(LateralSpread_1 = str_replace_all(LateralSpread_1, "free", "4")) %>%
  mutate(OffspringPerParent_1 = str_replace_all(OffspringPerParent_1, "<1", "1")) %>%
  mutate(OffspringPerParent_1 = str_replace_all(OffspringPerParent_1, "2-10", "2")) %>%
  mutate(OffspringPerParent_1 = str_replace_all(OffspringPerParent_1, ">10", "3")) %>%
  mutate(OffspringPerParent_1 = as.numeric(OffspringPerParent_1)) %>%
  mutate(LateralSpread_1 = as.numeric(LateralSpread_1)) %>%
  mutate(clonalIndex = OffspringPerParent_1 + LateralSpread_1, .keep = "unused") %>%
  mutate(Species = str_replace(Species, " ", "_")) %>%
  mutate(Species = str_extract(Species, "[:alpha:]+_[:alpha:]+")) %>%
  rename(name = "Species") %>%
  group_by(name) %>%
  summarise(across(clonalIndex, ~median(.x, na.rm = T)))

### * check completeness of traits ####
test <- traits %>%
  select(t, n, f, sla, seedmass, height, rmf, rld, srl, lateral, flowerBegin, flowerEnd, flowerDuration)
vis_miss(test, cluster = F, sort_miss = T)
#gg_miss_var(test)
gg_miss_case(test, order_cases = F)
(test2 <- traits %>%
    select(name, sla, seedmass, height, rmf, rld, srl, lateral, flowerBegin, flowerEnd, flowerDuration) %>%
    filter(!complete.cases(sla, seedmass, height, rmf, rld, srl, lateral, flowerBegin, flowerEnd, flowerDuration)))
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
traitsFBE <- traits %>%
  filter(group != "tree" & group != "shrub") %>%
  select(name, flowerBegin) %>%
  drop_na()
traitsFEN <- traits %>%
  filter(group != "tree" & group != "shrub") %>%
  select(name, flowerEnd) %>%
  drop_na()
traitsFDU <- traits %>%
  filter(group != "tree" & group != "shrub") %>%
  select(name, flowerDuration) %>%
  drop_na()
traitsAll <- traits %>%
  filter(group != "tree" & group != "shrub") %>%
  select(name, sla, seedmass, height, srl, flowerBegin) %>%
  drop_na()


### 7 CWM and FDis of functional plant traits #####################################################################################

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
sites$cwmAbuSla <- TdiversityAbu$CWM$sla %>%
  as.character() %>%
  as.numeric() %>%
  exp()

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
sites$cwmAbuSeedmass <- TdiversityAbu$CWM$seedmass %>%
  as.character() %>%
  as.numeric() %>%
  exp()

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
sites$cwmAbuHeight <- TdiversityAbu$CWM$height %>%
  as.character() %>%
  as.numeric() %>%
  exp()
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
  as.character() %>% 
  as.numeric() %>% 
  exp()
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
  as.character() %>% 
  as.numeric() %>% 
  exp()
length(traitsRMF$name) / (herbCount - undefinedSpeciesCount)

### g Flower begin -------------------------------------------------------------------------------------------
Tspecies <- semi_join(species, traitsFBE, by = "name")
Ttraits <- semi_join(traitsFBE, Tspecies, by = "name")
Tspecies <- Tspecies %>%
  pivot_longer(-name, "site", "value") %>%
  pivot_wider(site, name) %>%
  column_to_rownames("site")
Ttraits <- column_to_rownames(Ttraits, "name")
log_Ttraits <- log(Ttraits)
TdiversityAbu <- dbFD(log_Ttraits, Tspecies, w.abun = T, 
                      calc.FRic = F, calc.FDiv = F, corr = "sqrt");
sites$fdisAbuFbe <- TdiversityAbu$FDis
sites$cwmAbuFbe <- TdiversityAbu$CWM$flowerrBegin %>%
  as.character() %>% 
  as.numeric() %>% 
  exp()
length(traitsFBE$name) / (herbCount - undefinedSpeciesCount)

### h All -------------------------------------------------------------------------------------------
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
rm(list=setdiff(ls(), c("sites", "species", "traits")))
