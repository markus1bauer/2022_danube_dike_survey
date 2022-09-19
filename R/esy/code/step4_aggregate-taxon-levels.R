################### #
## Taxonomy      ####
################### #

library(tidyverse)

AGG <- parsing.result$species.aggs %>%
  stack() %>%
  mutate(ind = as.character(ind)) %>%
  filter(values != ind & values != "")

obs <- obs %>%
  left_join(AGG, by = c("TaxonName" = "values")) %>%
  mutate(TaxonName = if_else(is.na(ind), TaxonName, ind))
