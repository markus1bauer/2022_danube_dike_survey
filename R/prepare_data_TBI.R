# Prepare data for TBI ####
# Markus Bauer
# Citation: Markus Bauer, Jonathan Kiefer & Harald Albrecht (2020) Berichte der Bayerischen Botanischen Gesellschaft 90, 43-66



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Preparation ################################################################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Explanation: #
# edataT = 1984-vs.-2018 dataset,  edataB = 2002-vs.-2018 dataset


### Packages ###
library(tidyverse)

### Start ###
rm(list = ls())
setwd("Z:/Documents/0_Uni/2018_Projekt_9_Masterthesis/3_Aufnahmen_und_Ergebnisse/2020_BBBG/data/raw")

### Load data ###
edata <- read_csv("data_raw_environment.csv", col_names = T, na = "na", col_types = 
                    cols(
                      .default = col_double(),
                      ID = col_factor(),
                      plot = col_factor(),
                      plot_old = col_factor(),
                      block = col_factor(),
                      year = col_factor(),
                      dataset = col_factor(),
                      botanist = col_factor(),
                      photoName_June = col_factor(),
                      photoDate_June = col_date(),
                      photoName_April = col_factor(),
                      photoDate_April = col_date()
                    )        
)
edataT <- filter(edata, dataset == "transects")
edataB <- filter(edata, dataset == "blocks")

vdata <- read_csv("data_raw_species.csv", col_names = T, na = "na", col_types = 
                      cols(
                        .default = col_double(),
                        name = col_factor(),
                        abb = col_factor()
                      )
)

vdata <- vdata %>%
  filter(name != "Cuscuta_epithymum", name != "Acer_platanoides_juv") %>%
  select(!(starts_with("X73"))) %>%
  group_by(name) %>%
  mutate(sum = sum(c_across(X84I1:X18III16))) %>%
  mutate(presence = if_else(sum > 0, 1, 0)) %>%
  filter(presence == 1) %>%
  ungroup() %>%
  select(-sum, -presence)

vdataT <- select(vdata, "abb", matches("X84"), matches("X18I"))
vdataB <- select(vdata, "abb", matches("X03"), matches("X18N"), matches("X18M"),matches("X18S"))

## 1 Abundance data ####
vdataTabu <- column_to_rownames(vdataT, var = "abb")
vdataBabu <- column_to_rownames(vdataB, var = "abb")
vdataTabu <- as.data.frame(t(vdataTabu))
vdataBabu <- as.data.frame(t(vdataBabu))
vdataTabu <- vdataTabu[,colSums(vdataTabu) > 0]
vdataBabu <- vdataBabu[,colSums(vdataBabu) > 0]
#### 2 Presence-absence data####
vdataTpa <- vdataTabu
vdataTpa[vdataTpa > 0] = 1
vdataBpa <- vdataBabu
vdataBpa[vdataBpa > 0] = 1
####Site names
sitesT <- edataT %>%
  filter(year == 1984) %>% 
  select(plot)
sitesB <- edataB %>% 
  filter(year == 2003) %>% 
  select(plot)
rm(edataT,edataB)
####3 Ellenberg N=2 vs. N=3####
vdataTabu2 <- vdata %>% 
  filter(n == "2") %>% 
  select("abb", matches("X84"), matches("X18I"))
vdataBabu2 <- vdata %>% 
  filter(n == "2") %>% 
  select("abb", matches("X03"), matches("X18N"), matches("X18M"), matches("X18S"))
vdataTabu3 <- vdata %>% 
  filter(n == "3") %>% 
  select("abb", matches("X84"), matches("X18I"))
vdataBabu3 <- vdata %>% 
  filter(n == "3") %>%
  select("abb", matches("X03"), matches("X18N"), matches("X18M"), matches("X18S"))
vdataTabu2 <- column_to_rownames(vdataTabu2, var = "abb")
vdataBabu2 <- column_to_rownames(vdataBabu2, var = "abb")
vdataTabu3 <- column_to_rownames(vdataTabu3, var = "abb")
vdataBabu3 <- column_to_rownames(vdataBabu3, var = "abb")
vdataTabu2 <- as.data.frame(t(vdataTabu2))
vdataBabu2 <- as.data.frame(t(vdataBabu2))
vdataTabu3 <- as.data.frame(t(vdataTabu3))
vdataBabu3 <- as.data.frame(t(vdataBabu3))
vdataTabu2 <- vdataTabu2[,colSums(vdataTabu2) > 0]
vdataBabu2 <- vdataBabu2[,colSums(vdataBabu2) > 0]
vdataTabu3 <- vdataTabu3[,colSums(vdataTabu3) > 0]
vdataBabu3 <- vdataBabu3[,colSums(vdataBabu3) > 0]
vdataTpa2 <- vdataTabu2
vdataTpa2[vdataTpa2>0]=1
vdataBpa2 <- vdataBabu2
vdataBpa2[vdataBpa2>0]=1
vdataTpa3 <- vdataTabu3
vdataTpa3[vdataTpa3>0]=1
vdataBpa3 <- vdataBabu3
vdataBpa3[vdataBpa3>0]=1
#### 4 Ellenberg T=5 vs. T=6####
vdataTabu5 <- vdata %>% 
  filter(t == "5") %>%
  select("abb", matches("X84"), matches("X18I"))
vdataBabu5 <- vdata %>% 
  filter(t == "5") %>% 
  select("abb", matches("X03"), matches("X18N"), matches("X18M"), matches("X18S"))
vdataTabu6 <- vdata %>% 
  filter(t == "6") %>% 
  select("abb", matches("X84"), matches("X18I"))
vdataBabu6 <- vdata %>% 
  filter(t == "6") %>% 
  select("abb", matches("X03"), matches("X18N"), matches("X18M"), matches("X18S"))
vdataTabu5 <- column_to_rownames(vdataTabu5, var = "abb")
vdataBabu5 <- column_to_rownames(vdataBabu5, var = "abb")
vdataTabu6 <- column_to_rownames(vdataTabu6, var = "abb")
vdataBabu6 <- column_to_rownames(vdataBabu6, var = "abb")
vdataTabu5 <- as.data.frame(t(vdataTabu5))
vdataBabu5 <- as.data.frame(t(vdataBabu5))
vdataTabu6 <- as.data.frame(t(vdataTabu6))
vdataBabu6 <- as.data.frame(t(vdataBabu6))
vdataTabu5 <- vdataTabu5[,colSums(vdataTabu5) > 0]
vdataBabu5 <- vdataBabu5[,colSums(vdataBabu5) > 0]
vdataTabu6 <- vdataTabu6[,colSums(vdataTabu6) > 0]
vdataBabu6 <- vdataBabu6[,colSums(vdataBabu6) > 0]
vdataTpa5 <- vdataTabu5
vdataTpa5[vdataTpa5 > 0] = 1
vdataBpa5 <- vdataBabu5
vdataBpa5[vdataBpa5 > 0] = 1
vdataTpa6 <- vdataTabu6
vdataTpa6[vdataTpa6 > 0] = 1
vdataBpa6 <- vdataBabu6
vdataBpa6[vdataBpa6 > 0] = 1
#### 5 Graminoids-Herbs-Legumes####
vdataTabuG <- vdata %>% 
  filter(family == "Poaceae" | family == "Cyperaceae" | family == "Juncaceae") %>%
  select("abb", matches("X84"), matches("X18I"))
vdataTabuH <- vdata %>% 
  filter(family != "Poaceae" | family != "Cyperaceae" | family != "Juncaceae" | family != "Fabaceae") %>% 
  select("abb", matches("X84"), matches("X18I"))
vdataTabuL <- vdata %>% 
  filter(family == "Fabaceae") %>%
  select("abb", matches("X84"), matches("X18I"))
vdataBabuG <- vdata %>%
  filter(family == "Poaceae" | family == "Cyperaceae" | family == "Juncaceae") %>%
  select("abb", matches("X03"), matches("X18N"), matches("X18M"), matches("X18S"))
vdataBabuH <- vdata %>% 
  filter(family != "Poaceae" | family != "Cyperaceae" | family != "Juncaceae" | family != "Fabaceae") %>%
  select("abb", matches("X03"), matches("X18N"), matches("X18M"), matches("X18S"))
vdataBabuL <- vdata %>% 
  filter(family == "Fabaceae") %>%
  select("abb", matches("X03"), matches("X18N"), matches("X18M"), matches("X18S"))
vdataTabuG <- column_to_rownames(vdataTabuG, var = "abb")
vdataTabuH <- column_to_rownames(vdataTabuH, var = "abb")
vdataTabuL <- column_to_rownames(vdataTabuL, var = "abb")
vdataBabuG <- column_to_rownames(vdataBabuG, var = "abb")
vdataBabuH <- column_to_rownames(vdataBabuH, var = "abb")
vdataBabuL <- column_to_rownames(vdataBabuL, var = "abb")
vdataTabuG <- as.data.frame(t(vdataTabuG))
vdataBabuG <- as.data.frame(t(vdataBabuG))
vdataTabuH <- as.data.frame(t(vdataTabuH))
vdataBabuH <- as.data.frame(t(vdataBabuH))
vdataTabuL <- as.data.frame(t(vdataTabuL))
vdataBabuL <- as.data.frame(t(vdataBabuL))
vdataTabuG <- vdataTabuG[,colSums(vdataTabuG) > 0]
vdataBabuG <- vdataBabuG[,colSums(vdataBabuG) > 0]
vdataTabuH <- vdataTabuH[,colSums(vdataTabuH) > 0]
vdataBabuH <- vdataBabuH[,colSums(vdataBabuH) > 0]
vdataTabuL <- vdataTabuL[,colSums(vdataTabuL) > 0]
vdataBabuL <- vdataBabuL[,colSums(vdataBabuL) > 0]
vdataTpaG <- vdataTabuG
vdataTpaG[vdataTpaG > 0] = 1
vdataBpaG <- vdataBabuG
vdataBpaG[vdataBpaG > 0] = 1
vdataTpaH <- vdataTabuH
vdataTpaH[vdataTpaH > 0] = 1
vdataBpaH <- vdataBabuH
vdataBpaH[vdataBpaH > 0] = 1
vdataTpaL <- vdataTabuL
vdataTpaL[vdataTpaL > 0] = 1
vdataBpaL <- vdataBabuL
vdataBpaL[vdataBpaL > 0] = 1
