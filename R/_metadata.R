# Beta diversity on dike grasslands
# Prepare meta data ####
# Markus Bauer
# 2022-09-15


### Packages ###
library(here)
library(tidyverse)
library(EML)
library(emld)
#remotes::install_github("EDIorg/EMLassemblyline")
library(EMLassemblyline)

### Start ###
rm(list = ls())



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Collect metadata ##########################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



### 1 Methods and units #######################################################


### List of standard units, which should be used in metadata file ###
#EMLassemblyline::view_unit_dictionary()

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

unit_list <- set_unitList(custom_units)



### 2 Contact #################################################################


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

associated_party <- list(
  eml$associatedParty(
    individualName = eml$individualName(
      givenName = "Jakob K.",
      surName = "Huber"
      ),
    positionName = "Researcher",
    organizationName = "Technical University of Munich",
    electronicMailAddress = "jakob.huber@posteo.de"
  ),
  eml$associatedParty(
    individualName = eml$individualName(
      givenName = "Johannes",
      surName = "Kollmann"
    ),
    positionName = "Professor",
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
    onlineUrl = "https://www3.ls.tum.de/roek/mitarbeiter-in/prof-dr-johannes-kollmann/"
  )



### 3 Temporal and spatial coverage ###########################################


geographic_description <- "Danube dikes near Deggendorf"

coverage <- set_coverage(
  begin = "2017-06-01",
  end = "2021-07-31",
  sci_names = list(list(
    Subdivision = "Spermatophytina"
  )),
  geographicDescription = geographic_description,
  west = 12.58996,
  east = 13.1162,
  north = 48.90389,
  south = 48.67502,
  altitudeMin = 309,
  altitudeMaximum = 315,
  altitudeUnits = "meter"
)



### 4 Description #############################################################


pub_date <- "2022"

title <- "Beta diversity of restored river dike grasslands is strongly influenced by uncontrolled  spatio-temporal variability"

abstract <- "1.	Spatio-temporal dynamics of biodiversity are a key measure when monitoring restoration success. Balanced species turnover is aimed for because it increases overall biodiversity and improves ecosystem stability and multifunctionality, while nestedness should be avoided, as richness gradients indicate low biodiversity of certain restored sites. For predictive restoration, it is important to analyse beta diversity and to identify its control using site characteristics, spatial effects, historical factors and non-directional year effects.\\n
2.	We studied dike grasslands 4–19 years after restoration at River Danube in SE Germany over five years (2017–2021, 41 plots in 12 sites). We calculated beta-diversity indices to describe spatial variation and temporal turnover, including their additive components ‘replacement’ and ‘nestedness’, or ‘gains’ and ‘losses’.\\n
3.	The analysis of the spatial variation of the restored dike grasslands did not reveal homogenisation despite a significant temporal turnover, and was largely dominated by replacement-driven dissimilarity. The drivers of replacement changed over time, although replacement was mainly affected by exposition and spatial factors. Historical factors were inconsistent over time, and no statistically clear drivers were found for nestedness.\\n
4.	The dike grasslands exhibited on average 37 ± 11% year-to-year turnover in species composition, with some spatio-temporal variation. Gains and losses were balanced over time, although prevalences changed over time and were mostly pronounced on south-exposed slopes.\\n
5.	In conclusion, the restored grasslands exhibited spatio-temporal turnover controlled by climate and soil variability, varying management regimes, local stochastic (biotic) dynamics and landscape context. Thus, restoration targets defined as a single state should be extended by a desired variation in species composition. Furthermore, the dominance of replacement move the focus from searching the perfect fit for certain restoration targets to a variation of the approaches to increase the beta diversity within a restoration project.\\n
"

keyword_set <- list(
  list(
    keywordThesaurus = "LTER controlled vocabulary",
    keyword = list(
      "grasslands",
      "meadows",
      "monitoring",
      "permanent plots",
      "plant communities",
      "restoration",
      "rivers",
      "soil samples",
      "spatial properties",
      "species composition",
      "temporal properties",
      "vegetation dynamics"
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

intellectual_rights <- "CC-BY-4.0: https://creativecommons.org/licenses/by/4.0/deed.en"

alternate_identifier <- "https://doi.org/10.5281/zenodo.6107806"

short_name <- "Survey of old danube dikes"

language <- "English"

reference_publication <- "Bauer et al. (2023) EcoEvoRxiv DOI XXX"



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B finalize EML ##############################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



dataset <- list(
  title = title,
  shortName = short_name,
  pubDate = pub_date,
  creator = creator,
  associated_party = associated_party,
  intellectualRights = intellectual_rights,
  alternateIdentifier = alternate_identifier,
  abstract = abstract,
  keywordSet = keyword_set,
  coverage = coverage,
  referencePublication = reference_publication,
  language = language,
  contact = contact,
  additonalMetadata = list(metadata = list(
    unit_list = unit_list
  ))
)

eml <- list(
  packageId = uuid::UUIDgenerate(),
  system = "uuid", # type of identifier
  dataset = dataset
)

write_eml(eml, here("METADATA.xml"))
eml_validate(here("METADATA.xml"))

#emldown::render_eml(here("METADATA.xml"), open = TRUE, outfile = here("METADATA.html"), publish_mode = FALSE)
