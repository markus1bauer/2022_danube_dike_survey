
### 6 Load functional plant traits #####################################################################################

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
