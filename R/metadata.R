# Beta diversity on dike grasslands
# Prepare meta data ####
# Markus Bauer
# 2022-01-11
# Citation: 
## Bauer M, Huber J, Kollmann J (submitted) 
## Balanced turnover is a main aspect of biodiversity on restored dike grasslands: not only deterministic environmental effects, but also non-directional year and site effects drive spatial and temporal beta diversity.
## Unpublished data.


### Packages ###
library(here)
library(tidyverse)
library(EML)
library(emld)
#remotes::install_github("EDIorg/EMLassemblyline")
library(EMLassemblyline)

### Start ###
rm(list = ls())
setwd(here("data/raw"))



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Collect metadata ##############################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


### 1 Methods and units #####################################################################################

methods_file <- here("data/text/methods.odt")
methods <- set_methods(methods_file)

EMLassemblyline::view_unit_dictionary() # List of standard units, which should be used in metadata file

custom_units <- bind_rows(
  data.frame(id = "milligramPerDezigram", 
             unitType = "massPerMass", 
             parentSI = "gramPerGram", 
             multiplierToSI = 0.00001, 
             description = "milligram of element per 100 gram soil"),
  data.frame(id = "millimeterSquaredPerMilligram", 
             unitType = "specificArea", 
             parentSI = "meterPerGram", 
             multiplierToSI = 1, 
             description = "square millimeters per milligram")
  )

unitList <- set_unitList(custom_units)


### 2 Raw data #####################################################################################

### a data_raw_species  -------------------------------------------------------------------------------------------
attributes <- read_csv("data_raw_species_metadata.csv")

position <- c(
  m = "middle part of the slope",
  o = "upper part of the slope",
  u = "lower part of the slope"
)

factors <- bind_rows(
  data.frame(
    attributeName = "position",
    code = names(position),
    definition = unname(position)
  )
)

attributeList_raw_species <- set_attributes(attributes, 
                                          factors, 
                                          col_classes = c("character", 
                                                          "Date", 
                                                          "factor", 
                                                          "character",
                                                          "numeric", 
                                                          "numeric")
                                          )
                                          
physical_raw_species <- set_physical("data_raw_species.csv")

### b data_raw_traits  -------------------------------------------------------------------------------------------
setwd(here("data/raw"))
attributes <- read_csv("data_raw_traits_metadata.csv") %>%
  select(-type, -factor)

group <- c(
  top_grass = "taller grass species",
  sub_grass = "smaller grass species",
  legume = "species of the family Fabaceae",
  herb = "herbal species"
)
legal <- c(
  BG = "besonders geschuetzt",
  SG = "streng geschuetzt"
)
rlg <- c(
  V = "warning list",
  "3" = "vulnerable",
  "2" = "endangered",
  "1" = "critically endangered"
)
rlb <- c(
  V = "warning list",
  "3" = "vulnerable",
  "2" = "endangered",
  "1" = "critically endangered"
)
neophyte <- c(
  yes = "species is",
  no = "species is not"
)
fchange <- c(
  "1" = "changing moist conditions",
  "2" = "strong chaning moist conditions"
)
targetHerb <- c(
  yes = "species is",
  no = "species is not"
)
targetGrass <- c(
  yes = "species is",
  no = "species is not"
)
targetArrhenatherion <- c(
  yes = "species is",
  no = "species is not"
)
ffh6510 <- c(
  yes = "species is",
  no = "species is not"
)
ffh6210 <- c(
  yes = "species is",
  no = "species is not"
)
table30 <- c(
  yes = "species is",
  no = "species is not"
)
table33 <- c(
  yes = "species is",
  no = "species is not"
)
table34 <- c(
  yes = "species is",
  no = "species is not"
)
nitrogenIndicator <- c(
  yes = "species is",
  no = "species is not"
)
grazingIndicator <- c(
  yes = "species is",
  no = "species is not"
)
ruderalIndicator <- c(
  yes = "species is",
  no = "species is not"
)
others <- c(
  yes = "species is",
  no = "species is not"
)

factors <- bind_rows(
  data.frame(
    attributeName = "group",
    code = names(group),
    definition = unname(group)
  ),
  data.frame(
    attributeName = "legal",
    code = names(legal),
    definition = unname(legal)
  ),
  data.frame(
    attributeName = "rlg",
    code = names(rlg),
    definition = unname(rlg)
  ),
  data.frame(
    attributeName = "rlb",
    code = names(rlb),
    definition = unname(rlb)
  ),
  data.frame(
    attributeName = "neophyte",
    code = names(neophyte),
    definition = unname(neophyte)
  ),
  data.frame(
    attributeName = "fchange",
    code = names(fchange),
    definition = unname(fchange)
  ),
  data.frame(
    attributeName = "targetHerb",
    code = names(targetHerb),
    definition = unname(targetHerb)
  ),
  data.frame(
    attributeName = "targetGrass",
    code = names(targetGrass),
    definition = unname(targetGrass)
  ),
  data.frame(
    attributeName = "targetArrhenatherion",
    code = names(targetArrhenatherion),
    definition = unname(targetArrhenatherion)
  ),
  data.frame(
    attributeName = "ff6510",
    code = names(ffh6510),
    definition = unname(ffh6510)
  ),
  data.frame(
    attributeName = "ffh6210",
    code = names(ffh6210),
    definition = unname(ffh6210)
  ),
  data.frame(
    attributeName = "table30",
    code = names(table30),
    definition = unname(table30)
  ),
  data.frame(
    attributeName = "table33",
    code = names(table33),
    definition = unname(table33)
  ),
  data.frame(
    attributeName = "table34",
    code = names(table34),
    definition = unname(table34)
  ),
  data.frame(
    attributeName = "nitrogenIndicator",
    code = names(nitrogenIndicator),
    definition = unname(nitrogenIndicator)
  ),
  data.frame(
    attributeName = "grazingIndicator",
    code = names(grazingIndicator),
    definition = unname(grazingIndicator)
  ),
  data.frame(
    attributeName = "ruderalIndicator",
    code = names(ruderalIndicator),
    definition = unname(ruderalIndicator)
  ),
  data.frame(
    attributeName = "others",
    code = names(others),
    definition = unname(others)
  )
)

attributeList_raw_traits <- set_attributes(attributes, 
                                        factors, 
                                        col_classes = c("character",
                                                        "character",
                                                        "character",
                                                        "character",
                                                        "factor",
                                                        "factor",
                                                        "factor",
                                                        "factor",
                                                        "factor",
                                                        "numeric",
                                                        "numeric",
                                                        "numeric",
                                                        "numeric",
                                                        "factor",
                                                        "numeric",
                                                        "numeric",
                                                        "factor",
                                                        "factor",
                                                        "factor",
                                                        "factor",
                                                        "factor",
                                                        "factor",
                                                        "factor",
                                                        "factor",
                                                        "factor",
                                                        "factor",
                                                        "factor",
                                                        "factor"
                                                        )
                                        )
physical_raw_traits <- set_physical("data_raw_traits.csv")

### c data_raw_sites  -------------------------------------------------------------------------------------------
setwd(here("data/raw"))
attributes <- read_csv("data_raw_sites_metadata.csv") %>%
  select(-type, -factor)

side <- c(
  water = "water side of the dike",
  land = "land side of the dike"
)
exposition <- c(
  north = "north exposition",
  south = "south exposition",
  east = "east exposition",
  west = "west exposition"
)
ageCategory<- c(
  new = "dike construction in 2016 an younger",
  old = "dike construction in 2015 an older"
)
HCl<- c(
  C1 = "C1",
  C2 = "C2",
  C3 = "C3",
  C4 = "C4"
)
humusLevel <- c(
  h1 = "h1",
  h2 = "h2",
  h3 = "h3",
  h4 = "h4",
  h5 = "h5"
)
cnLevel <- c(
  hq4 = "hq4",
  hq5 = "hq5"
)
phosphorusClass <- c(
  A = "very low",
  B = "low",
  C = "optimum",
  D = "high",
  E = "very high"
)
potassiumClass <- c(
  A = "very low",
  B = "low",
  C = "optimum",
  D = "high",
  E = "very high"
)
magnesiumClass <- c(
  A = "very low",
  B = "low",
  C = "optimum",
  D = "high",
  E = "very high"
)
sceletonLevel <- c(
  "1" = "stones and gravel < 2 vol%",
  "2" = "stones and gravel 2-10 vol%",
  "3" = "stones and gravel 10-25 vol%",
  "4" = "stones and gravel 25-50 vol%",
  "5" = "stones and gravel 50-75 vol%",
  "6" = "stones and gravel >75 vol%"
)

factors <- bind_rows(
  data.frame(
    attributeName = "side",
    code = names(side),
    definition = unname(side)
  ),
  data.frame(
    attributeName = "exposition",
    code = names(exposition),
    definition = unname(exposition)
  ),
  data.frame(
    attributeName = "ageCategory",
    code = names(ageCategory),
    definition = unname(ageCategory)
  ),
  data.frame(
    attributeName = "HCl",
    code = names(HCl),
    definition = unname(HCl)
  ),
  data.frame(
    attributeName = "humusLevel",
    code = names(humusLevel),
    definition = unname(humusLevel)
  ),
  data.frame(
    attributeName = "cnLevel",
    code = names(cnLevel),
    definition = unname(cnLevel)
  ),
  data.frame(
    attributeName = "phosphorusClass",
    code = names(phosphorusClass),
    definition = unname(phosphorusClass)
  ),
  data.frame(
    attributeName = "potassiumClass",
    code = names(potassiumClass),
    definition = unname(potassiumClass)
  ),
  data.frame(
    attributeName = "magnesiumClass",
    code = names(magnesiumClass),
    definition = unname(magnesiumClass)
  ),
  data.frame(
    attributeName = "sceletonLevel",
    code = names(sceletonLevel),
    definition = unname(sceletonLevel)
  )
)

attributeList_raw_sites <- set_attributes(attributes, 
                                            factors, 
                                            col_classes = c("character",
                                                            "character",
                                                            "character",
                                                            "numeric",
                                                            "numeric",
                                                            "numeric",
                                                            "factor",
                                                            "factor",
                                                            "Date",
                                                            "factor",
                                                            "Date",
                                                            "Date",
                                                            "Date",
                                                            "character",
                                                            "character",
                                                            "character",
                                                            "numeric",
                                                            "numeric",
                                                            "numeric",
                                                            "numeric",
                                                            "numeric",
                                                            "numeric",
                                                            "numeric",
                                                            "numeric",
                                                            "numeric",
                                                            "numeric",
                                                            "numeric",
                                                            "numeric",
                                                            "numeric",
                                                            "numeric",
                                                            "numeric",
                                                            "factor",
                                                            "numeric",
                                                            "numeric",
                                                            "numeric",
                                                            "numeric",
                                                            "factor",
                                                            "numeric",
                                                            "numeric",
                                                            "factor",
                                                            "numeric",
                                                            "numeric",
                                                            "numeric",
                                                            "numeric",
                                                            "numeric",
                                                            "factor",
                                                            "numeric",
                                                            "factor",
                                                            "numeric",
                                                            "factor",
                                                            "numeric",
                                                            "numeric",
                                                            "numeric",
                                                            "numeric",
                                                            "numeric",
                                                            "numeric",
                                                            "numeric",
                                                            "factor",
                                                            "numeric",
                                                            "numeric"
                                            )
)


physical_raw_sites <- set_physical("data_raw_sites.csv")

### 3 Processed data #####################################################################################

### a data_processed_species  -------------------------------------------------------------------------------------------
setwd(here("data/processed"))
attributes <- read_csv("data_processed_species_metadata.csv") %>%
  select(-type, -factor)

attributeList_processed_species <- set_attributes(attributes, 
                                                  col_classes = c("character", 
                                                                  "numeric", 
                                                                  "numeric")
                                                  )

physical_processed_species <- set_physical("data_processed_species.csv")


### b data_processed_traits  -------------------------------------------------------------------------------------------
setwd(here("data/processed"))
attributes <- read_csv("data_processed_traits_metadata.csv") %>%
  select(-type, -factor)

group <- c(
  top_grass = "taller grass species",
  sub_grass = "smaller grass species",
  legume = "species of the family Fabaceae",
  herb = "herbal species"
)
legal <- c(
  BG = "besonders geschuetzt",
  SG = "streng geschuetzt"
)
rlg <- c(
  V = "warning list",
  "3" = "vulnerable",
  "2" = "endangered",
  "1" = "critically endangered"
)
rlb <- c(
  V = "warning list",
  "3" = "vulnerable",
  "2" = "endangered",
  "1" = "critically endangered"
)
neophyte <- c(
  yes = "species is",
  no = "species is not"
)
fchange <- c(
  "1" = "changing moist conditions",
  "2" = "strong chaning moist conditions"
)
targetHerb <- c(
  yes = "species is",
  no = "species is not"
)
targetGrass <- c(
  yes = "species is",
  no = "species is not"
)
targetArrhenatherion <- c(
  yes = "species is",
  no = "species is not"
)
ffh6510 <- c(
  yes = "species is",
  no = "species is not"
)
ffh6210 <- c(
  yes = "species is",
  no = "species is not"
)
table30 <- c(
  yes = "species is",
  no = "species is not"
)
table33 <- c(
  yes = "species is",
  no = "species is not"
)
table34 <- c(
  yes = "species is",
  no = "species is not"
)
nitrogenIndicator <- c(
  yes = "species is",
  no = "species is not"
)
grazingIndicator <- c(
  yes = "species is",
  no = "species is not"
)
ruderalIndicator <- c(
  yes = "species is",
  no = "species is not"
)
others <- c(
  yes = "species is",
  no = "species is not"
)

factors <- bind_rows(
  data.frame(
    attributeName = "group",
    code = names(group),
    definition = unname(group)
  ),
  data.frame(
    attributeName = "legal",
    code = names(legal),
    definition = unname(legal)
  ),
  data.frame(
    attributeName = "rlg",
    code = names(rlg),
    definition = unname(rlg)
  ),
  data.frame(
    attributeName = "rlb",
    code = names(rlb),
    definition = unname(rlb)
  ),
  data.frame(
    attributeName = "neophyte",
    code = names(neophyte),
    definition = unname(neophyte)
  ),
  data.frame(
    attributeName = "fchange",
    code = names(fchange),
    definition = unname(fchange)
  ),
  data.frame(
    attributeName = "targetHerb",
    code = names(targetHerb),
    definition = unname(targetHerb)
  ),
  data.frame(
    attributeName = "targetGrass",
    code = names(targetGrass),
    definition = unname(targetGrass)
  ),
  data.frame(
    attributeName = "targetArrhenatherion",
    code = names(targetArrhenatherion),
    definition = unname(targetArrhenatherion)
  ),
  data.frame(
    attributeName = "ff6510",
    code = names(ffh6510),
    definition = unname(ffh6510)
  ),
  data.frame(
    attributeName = "ffh6210",
    code = names(ffh6210),
    definition = unname(ffh6210)
  ),
  data.frame(
    attributeName = "table30",
    code = names(table30),
    definition = unname(table30)
  ),
  data.frame(
    attributeName = "table33",
    code = names(table33),
    definition = unname(table33)
  ),
  data.frame(
    attributeName = "table34",
    code = names(table34),
    definition = unname(table34)
  ),
  data.frame(
    attributeName = "nitrogenIndicator",
    code = names(nitrogenIndicator),
    definition = unname(nitrogenIndicator)
  ),
  data.frame(
    attributeName = "grazingIndicator",
    code = names(grazingIndicator),
    definition = unname(grazingIndicator)
  ),
  data.frame(
    attributeName = "ruderalIndicator",
    code = names(ruderalIndicator),
    definition = unname(ruderalIndicator)
  ),
  data.frame(
    attributeName = "others",
    code = names(others),
    definition = unname(others)
  )
)

attributeList_raw_traits <- set_attributes(attributes, 
                                           factors, 
                                           col_classes = c("character",
                                                           "character",
                                                           "character",
                                                           "character",
                                                           "factor",
                                                           "factor",
                                                           "factor",
                                                           "factor",
                                                           "factor",
                                                           "numeric",
                                                           "numeric",
                                                           "numeric",
                                                           "numeric",
                                                           "factor",
                                                           "numeric",
                                                           "numeric",
                                                           "factor",
                                                           "factor",
                                                           "factor",
                                                           "factor",
                                                           "factor",
                                                           "factor",
                                                           "factor",
                                                           "factor",
                                                           "factor",
                                                           "factor",
                                                           "factor",
                                                           "factor"
                                           )
)

physical_processed_traits <- set_physical("data_raw_traits.csv")


### c data_processed_sites  -------------------------------------------------------------------------------------------
attributes <- read_csv("data_processed_sites_metadata.csv")

position <- c(
  m = "middle part of the slope",
  u = "upper part of the slope",
  l = "lower part of the slope"
)

factors <- bind_rows(
  data.frame(
    attributeName = "position",
    code = names(position),
    definition = unname(position)
  )
)

attributeList_processed_sites <- set_attributes(attributes, 
                                      factors, 
                                      col_classes = c("character", 
                                                      "Date", 
                                                      "factor", 
                                                      "character",
                                                      "numeric", 
                                                      "numeric")
)

physical_processed_sites <- set_physical("data_raw_sites.csv")


### 4 Put data table together #####################################################################################

dataTable <- list(
  list(
    entityName = "data_raw_species.csv",
    entityDescription = "raw species abundances",
    #physical = physical_raw_species,
    attributeList = attributeList_raw_species
  ),
  list(
    entityName = "data_raw_traits.csv",
    entityDescription = "raw plant trait list",
    #physical = physical_raw_traits,
    attributeList = attributeList_raw_traits
  ),
  list(
    entityName = "data_raw_sites.csv",
    entityDescription = "environmental raw data of the sites",
    #physical = physical_raw_sites,
    attributeList = attributeList_raw_sites
  ),
  list(
    entityName = "data_processed_species.csv",
    entityDescription = "processed species abundances",
    #physical = physical_processed_species,
    attributeList = attributeList_processed_species
  )#,
  #list(
    #entityName = "data_processed_traits.csv",
    #entityDescription = "processed plant trait list",
    #physical = physical_processed_traits,
    #attributeList = attributeList_processed_traits
  #),
  #list(
    #entityName = "data_processed_sites.csv",
    #entityDescription = "environmental processed data of the sites",
    #physical = physical_processed_sites,
    #attributeList = attributeList_processed_sites
  #)
)


### 5 Contact #####################################################################################

address <- list(
  deliveryPoint = "Emil-Ramann-Strasse 6",
  city = "Freising",
  administrativeArea = "Bayern",
  postalCode = "85354",
  country = "Germany")

creator <- eml$creator(
  individualName = eml$individualName(
    givenName = "Markus", 
    surName = "Bauer"
  ),
  positionName = "PhD student",
  organizationName = "Technical University of Munich",
  address = address,
  electronicMailAddress = "markusbauer@mailbox.org",
  phone = "0049-152-56391781",
  id = "https://orcid.org/0000-0001-5372-4174"
)

associatedParty <- list(
  eml$associatedParty(
    individualName = eml$individualName(
      givenName = "Jakob", 
      surName = "Huber"
    ),
    role = "Researcher",
    organizationName = "Technical University of Munich",
    electronicMailAddress = "jakob.huber@posteo.de"
  ),
  eml$associatedParty(
    individualName = eml$individualName(
      givenName = "Johannes", 
      surName = "Kollmann"
    ),
    role = "Professor",
    organizationName = "Technical University of Munich",
    address = address,
    electronicMailAddress = "jkollmann@wzw.tum.de",
    phone = "0049-8161-714144",
    id = "https://orcid.org/0000-0002-4990-3636"
  )
)

contact <- 
  list(
    individualName = creator$individualName,
    electronicMailAddress = creator$electronicMailAddress,
    address = address,
    organizationName = "Technical University of Munich",
    onlineUrl = "DOI address to the database"
  )


### 6 Temporal and spatial coverage #####################################################################################

geographicDescription <- "Danube dikes near Deggendorf"

coverage <- set_coverage(
  begin = "2017-06-01", end = "2021-07-31",
  sci_names = list(list(
    Subdivision = "Spermatophytina"
  )),
  geographicDescription = geographicDescription,
  west = 12.58996, east = 13.1162,
  north = 48.90389, south = 48.67502,
  altitudeMin = 309, altitudeMaximum = 315,
  altitudeUnits = "meter"
)


### 7 Description #####################################################################################

pubDate = "2022"

title = "Danube old dikes"

abstract <- "Not written yet"

keywordSet <- list(
  list(
    keywordThesaurus = "LTER controlled vocabulary",
    keyword = list("rivers",
                   "vegetation dynamics",
                   "restoration")
  ),
  list(
    keywordThesaurus = "own vocabulary",
    keyword = list("beta diversity",
                   "temperate grassland",
                   "dike")
  )
)

intellectualRights <- "CC-BY-4.0: https://creativecommons.org/licenses/by/4.0/deed.en"



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B finalize EML ##############################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


dataset <- list(
  title = title,
  pubDate = pubDate,
  creator = creator,
  associatedParty = associatedParty,
  intellectualRights = intellectualRights,
  abstract = abstract,
  keywordSet = keywordSet,
  coverage = coverage,
  contact = contact,
  methods = methods,
  dataTable = dataTable,
  additonalMetadata = list(metadata = list(
    unitList = unitList
  ))
  )

eml <- list(
  packageId = uuid::UUIDgenerate(),
  system = "uuid", # type of identifier
  dataset = dataset
  )

write_eml(eml, here("METADATA.xml"))
eml_validate(here("METADATA.xml"))

emldown::render_eml(here("METADATA.xml"), open = T, outfile = here("METADATA.html"), publish_mode = F)

