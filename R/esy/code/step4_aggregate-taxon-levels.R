################### #
## Taxonomy      ####
################### #

aggs <- parsing.result$species.aggs
AGG <- stack(aggs)
AGG$ind <- as.character(AGG$ind)
AGG <- AGG[AGG$values != AGG$ind,]
AGG <- AGG[AGG$values != '',]
index1 <- match(obs$TaxonName, AGG$values)
obs$TaxonName[!is.na(index1)] <- AGG$ind[index1[!is.na(index1)]]
