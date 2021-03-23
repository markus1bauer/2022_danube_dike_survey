# Prepare Metadata ####
# Markus Bauer


### Packages ###
library(here)
library(EML)

### Start ###
rm(list = ls())



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Create infos about persons ##############################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


creator <- eml$creator(
  individualName = eml$individualName(
    givenName = "Markus", 
    surName = "Bauer"),
  electronicMailAddress = "markusbauer@mailbox.org"
)

associatedParty <- list(
  individualName = eml$individualName(
    givenName = "Jakob", 
    surName = "Huber")
)

address <- list(
  deliveryPoint = "Emil-Ramann-Strasse 6",
  city = "Freising",
  administrativeArea = "Bayern",
  postalCode = "85354",
  country = "Germany")

contact <- 
  list(
    individualName = creator$individualName,
    electronicMailAddress = creator$electronicMailAddress,
    address = address,
    organizationName = "Technical University of Munich",
    phone = "0049-152-56391781"
    )



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Create site and experiment infos ##############################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


abstract <- "Not written yet"

keywordSet <- list(
    keywordThesaurus = "LTER controlled vocabulary",
    keyword = list("rivers",
                   "vegetation dynamics",
                   "restoration")
    )

geographicDescription <- "Danube dikes near Deggendorf"

coverage <- set_coverage(
  begin = "2017-06-01", end = "2021-07-31",
  sci_names = list(list(
    Kingdom = "Plantae",
    Division = "Tracheophyta",
    Subdivision = "Spermatophytina"
    )),
  geographicDescription = geographicDescription,
  west = 12.58996, east = 13.1162,
  north = 48.90389, south = 48.67502,
  altitudeMin = 309, altitudeMaximum = 315,
  altitudeUnits = "meter"
  )



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# C finalize EML ##############################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


dataset <- list(
    title = "Danube old dikes",
    creator = creator,
    pubDate = "2022",
    language = "English",
    intellectualRights = "CC BY 4.0",
    abstract = abstract,
    keywordSet = keywordSet,
    coverage = coverage,
    contact = contact
    )

eml <- list(
  packageId = uuid::UUIDgenerate(),
  system = "uuid", # type of identifier
  dataset = dataset
  )

setwd(here())
write_eml(eml, "METADATA.xml")
eml_validate("METADATA.xml")
