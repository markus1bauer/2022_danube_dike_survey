# Beta diversity on dike grasslands
# Prepare meta data ####
# Markus Bauer
# 2022-02-11


### Packages ###
library(here)
library(tidyverse)
library(EML)
library(emld)
#remotes::install_github("EDIorg/EMLassemblyline")

### Start ###
rm(list = ls())
setwd(here("data", "raw"))



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Collect metadata ##########################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


### 1 Methods and units #######################################################

methods_file <- here("data/text/methods.odt")
methods <- set_methods(methods_file)

### List of standard units, which should be used in metadata file ###
EMLassemblyline::view_unit_dictionary()

custom_units <- bind_rows(
  data.frame(
    id = "milligramPerDezigram",
    unitType = "massPerMass",
    parentSI = "gramPerGram",
    multiplierToSI = 0.00001,
    description = "milligram of element per 100 gram soil"
  ),
  data.frame(
    id = "millimeterSquaredPerMilligram",
    unitType = "specificArea",
    parentSI = "meterPerGram",
    multiplierToSI = 1,
    description = "square millimeters per milligram"
  )
)

unitList <- set_unitList(custom_units)


### 2 Raw data ################################################################

### a data_raw_species  -------------------------------------------------------
attributes <- read_csv("data_raw_species_metadata.csv")
physical_raw_species <- set_physical("data_raw_species.csv")

### b data_raw_traits  --------------------------------------------------------
attributes <- read_csv("data_raw_traits_metadata.csv") %>%
  select(-type, -factor)
physical_raw_traits <- set_physical("data_raw_traits.csv")

### c data_raw_sites  ---------------------------------------------------------
attributes <- read_csv("data_raw_sites_metadata.csv") %>%
  select(-type, -factor)
physical_raw_sites <- set_physical("data_raw_sites.csv")


### 3 Processed data ##########################################################

### a data_processed_species  -------------------------------------------------

attributes <- read_csv(
  here("data/processed/data_processed_species_metadata.csv")
  ) %>%
  select(-type, -factor)
physical_processed_species <- set_physical("data_processed_species.csv")


### b data_processed_traits  --------------------------------------------------
attributes <- read_csv(
  here("data/processed/data_processed_traits_metadata.csv")
  ) %>%
  select(-type, -factor)
physical_processed_traits <- set_physical("data_raw_traits.csv")


### c data_processed_sites  ---------------------------------------------------
attributes <- read_csv(here("data/processed/data_processed_sites_metadata.csv"))
physical_processed_sites <- set_physical("data_raw_sites.csv")


### 4 Put data table together #################################################

dataTable <- list(
  list(
    entityName = "data_raw_species.csv",
    # entityDescription = "raw species abundances",
    physical = physical_raw_species,
    # attributeList = attributeList_raw_species
  ),
  list(
    entityName = "data_raw_traits.csv",
    # entityDescription = "raw plant trait list",
    physical = physical_raw_traits,
    # attributeList = attributeList_raw_traits
  ),
  list(
    entityName = "data_raw_sites.csv",
    # entityDescription = "environmental raw data of the sites",
    # physical = physical_raw_sites,
    # attributeList = attributeList_raw_sites
  ),
  list(
    entityName = "data_processed_species.csv",
    # entityDescription = "processed species abundances",
    physical = physical_processed_species,
    # attributeList = attributeList_processed_species
  ),
  list(
    entityName = "data_processed_traits.csv",
    # entityDescription = "processed plant trait list",
    physical = physical_processed_traits,
    # attributeList = attributeList_processed_traits
  ),
  list(
    entityName = "data_processed_sites.csv",
    # entityDescription = "environmental processed data of the sites",
    physical = physical_processed_sites,
    # attributeList = attributeList_processed_sites
  )
)


### 5 Contact #################################################################

address <- list(
  deliveryPoint = "Emil-Ramann-Strasse 6",
  city = "Freising",
  administrativeArea = "Bayern",
  postalCode = "85354",
  country = "Germany"
)

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
    electronicMailAddress = "johannes.kollmann@tum.de",
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
    onlineUrl = "https://doi.org/10.5281/zenodo.6334100"
  )


### 6 Temporal and spatial coverage ###########################################

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


### 7 Description #############################################################

pubDate <- "2022"

title <- "Danube old dikes"

abstract <- "Not written yet"

keywordSet <- list(
  list(
    keywordThesaurus = "LTER controlled vocabulary",
    keyword = list(
      "rivers",
      "vegetation dynamics",
      "restoration"
    )
  ),
  list(
    keywordThesaurus = "own vocabulary",
    keyword = list(
      "beta diversity",
      "temperate grassland",
      "dike"
    )
  )
)

intellectualRights <- "CC-BY-4.0: https://creativecommons.org/licenses/by/4.0/deed.en"



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B finalize EML ##############################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


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
  # dataTable = dataTable,
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

# emldown::render_eml(here("METADATA.xml"), open = TRUE, outfile = here("METADATA.html"), publish_mode = FALSE)
