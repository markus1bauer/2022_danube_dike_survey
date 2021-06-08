# Prepare Metadata ####
# Markus Bauer


### Packages ###
library(here)
library(EML)
library(emld)

### Start ###
rm(list = ls())
setwd(here())



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Metadata ##############################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


### 1 Contact #####################################################################################

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
    organizationName = "Technical University of Munich"
    )


### 2 Coverage #####################################################################################

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


### 3 Description #####################################################################################

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

intellectualRights <- "CC-BY-4.0"



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

render_eml("METADATA.xml", open = F, outfile = "METADATA.html", publish_mode = F)

