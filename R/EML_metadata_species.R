# Prepare Metadata ####
# Markus Bauer


### Packages ###
library(here)
library(EML)
library(emld)
#remotes::install_github("ropenscilabs/emldown", build = F)
library("emldown")

### Start ###
rm(list = ls())
setwd(here())

### Load data ###
metadata <- read_eml("METADATA.xml")
metadata <- metadata$dataset
setwd(here("data/raw"))


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Metadata ##############################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


### 1 Species #####################################################################################

attributes <- read_csv("data_raw_species_metadata.csv")

position <- c(
  m = "middle part of the slope",
  u = "upper part of the slope",
  l = "lower part of the slope"
  )

factors <- rbind(
  data.frame(
    attributeName = "position",
    code = names(position),
    definition = unname(position)
  )
)

attributeList_species <- set_attributes(attributes, 
                                factors, 
                                col_classes = c("character", 
                                                "Date", 
                                                "factor", 
                                                "character",
                                                "numeric", 
                                                "numeric")
                                )

physical_species <- set_physical("data_raw_species.csv")


### 2 Traits #####################################################################################

attributes <- read_csv("data_raw_traits_metadata.csv")

position <- c(
  m = "middle part of the slope",
  u = "upper part of the slope",
  l = "lower part of the slope"
)

factors <- rbind(
  data.frame(
    attributeName = "position",
    code = names(position),
    definition = unname(position)
  )
)

attributeList_traits <- set_attributes(attributes, 
                                        factors, 
                                        col_classes = c("character", 
                                                        "Date", 
                                                        "factor", 
                                                        "character",
                                                        "numeric", 
                                                        "numeric")
)

physical_traits <- set_physical("data_raw_traits.csv")


### 3 Sites #####################################################################################

attributes <- read_csv("data_raw_sites_metadata.csv")

position <- c(
  m = "middle part of the slope",
  u = "upper part of the slope",
  l = "lower part of the slope"
)

factors <- rbind(
  data.frame(
    attributeName = "position",
    code = names(position),
    definition = unname(position)
  )
)

attributeList_sites <- set_attributes(attributes, 
                                       factors, 
                                       col_classes = c("character", 
                                                       "Date", 
                                                       "factor", 
                                                       "character",
                                                       "numeric", 
                                                       "numeric")
)

physical_sites <- set_physical("data_raw_sites.csv")


### 4 Put data table together #####################################################################################

dataTable <- list(
  list(
    entityName = "data_raw_species.csv",
    entityDescription = "species abundances",
    #physical = physical,
    attributeList = attributeList
    ),
  list(
    entityName = "data_raw_traits.csv",
    entityDescription = "plant trait list",
    #physical = physical_traits,
    attributeList = attributeList_traits
  ),
  list(
    entityName = "data_raw_sites.csv",
    entityDescription = "environmental data of the sites",
    #physical = physical_sites,
    attributeList = attributeList_sites
    )
  )



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
  contact = metadata$contact,
  dataTable = dataTable
  )

eml <- list(
  packageId = uuid::UUIDgenerate(),
  system = "uuid", # type of identifier
  dataset = dataset
  )

setwd(here("data/raw"))
write_eml(eml, "data_raw_species_metadata.xml")
eml_validate("data_raw_species_metadata.xml")

render_eml("data_raw_species_metadata.xml", open = T, outfile = "data_raw_species_metadata.html", publish_mode = F)

