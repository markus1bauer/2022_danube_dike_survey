# Prepare Metadata ####
# Markus Bauer


### Packages ###
library(here)
library(EML)
library(emld)

### Start ###
rm(list = ls())
setwd(here())

### Load data ###
metadata <- read_eml("METADATA.xml")
metadata <- metadata$dataset



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Metadata ##############################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

attributes <-
  tibble::tribble(
    ~attributeName, ~attributeDefinition,                                                 ~formatString, ~definition,        ~unit,   ~numberType,
    "run.num",    "which run number (=block). Range: 1 - 6. (integer)",                 NA,            "which run number", NA,       NA,
    "year",       "year, 2012",                                                         "YYYY",        NA,                 NA,       NA,
    "day",        "Julian day. Range: 170 - 209.",                                      "DDD",         NA,                 NA,       NA,
    "hour.min",   "hour and minute of observation. Range 1 - 2400 (integer)",           "hhmm",        NA,                 NA,       NA,
    "i.flag",     "is variable Real, Interpolated or Bad (character/factor)",           NA,            NA,                 NA,       NA,
    "variable",   "what variable being measured in what treatment (character/factor).", NA,            NA,                 NA,       NA,
    "value.i",    "value of measured variable for run.num on year/day/hour.min.",       NA,            NA,                 NA,       NA,
    "length",    "length of the species in meters (dummy example of numeric data)",     NA,            NA,                 "meter",  "real")

i.flag <- c(R = "real",
            I = "interpolated",
            B = "bad")
variable <- c(
  control  = "no prey added",
  low      = "0.125 mg prey added ml-1 d-1",
  med.low  = "0,25 mg prey added ml-1 d-1",
  med.high = "0.5 mg prey added ml-1 d-1",
  high     = "1.0 mg prey added ml-1 d-1",
  air.temp = "air temperature measured just above all plants (1 thermocouple)",
  water.temp = "water temperature measured within each pitcher",
  par       = "photosynthetic active radiation (PAR) measured just above all plants (1 sensor)"
)

value.i <- c(
  control  = "% dissolved oxygen",
  low      = "% dissolved oxygen",
  med.low  = "% dissolved oxygen",
  med.high = "% dissolved oxygen",
  high     = "% dissolved oxygen",
  air.temp = "degrees C",
  water.temp = "degrees C",
  par      = "micromoles m-1 s-1"
)

## Write these into the data.frame format
factors <- rbind(
  data.frame(
    attributeName = "i.flag",
    code = names(i.flag),
    definition = unname(i.flag)
  ),
  data.frame(
    attributeName = "variable",
    code = names(variable),
    definition = unname(variable)
  ),
  data.frame(
    attributeName = "value.i",
    code = names(value.i),
    definition = unname(value.i)
  )
)

attributeList <- set_attributes(attributes, factors, col_classes = c("character", 
                                                                     "Date", 
                                                                     "Date", 
                                                                     "Date", 
                                                                     "factor", 
                                                                     "factor", 
                                                                     "factor", 
                                                                     "numeric")
                                )

physical <- set_physical("hf205-01-TPexp1.csv")

dataTable <- list(
  entityName = "hf205-01-TPexp1.csv",
  entityDescription = "tipping point experiment 1",
  physical = physical,
  attributeList = attributeList)



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# C finalize EML ##############################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


dataset <- list(
  title = metadata$title,
  pubDate = metadata$pubDate,
  creator = metadata$creator,
  associatedParty = metadata$associatedParty,
  intellectualRights = metadata$intellectualRights,
  abstract = metadata$abstract,
  keywordSet = metadata$keywordSet,
  coverage = metadata$coverage,
  contact = metadata$contact
  )

eml <- list(
  packageId = uuid::UUIDgenerate(),
  system = "uuid", # type of identifier
  dataset = dataset
  )

setwd(here())
write_eml(eml, "METADATA.xml")
eml_validate("METADATA.xml")
