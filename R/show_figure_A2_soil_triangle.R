# Show map of the Danube dike experiment ####
# Markus Bauer
# Citation: Markus Bauer, Jakob Huber, Katharina Beck, Johannes Kollmann (xxxx)



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Preparation ################################################################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

### Packages ###
library(here)
library(tidyverse)
library(soiltexture)

### Start ###
rm(list = ls())


### Load data ###
selection <- read_csv(here("data/processed/data_processed_sites.csv"), col_names = T, na = "na", col_types = 
                    cols(
                      .default = "?"
                    )) %>%
  filter(accumulatedCov > 0) %>%
  add_count(plot) %>%
  filter(surveyYear == 2017 & n == max(n)) %>%
  select(plot)
sites <- read_csv(here("data/raw/data_raw_sites.csv"), col_names = T, na = "na", col_types = 
                     cols(
                       .default = "?",
                       vegetationCov_2017 = "d"
                     )) %>%
  select(id, sandPerc, siltPerc, clayPerc, NtotalPerc) %>%
  mutate(id = str_extract(id, "\\d\\d")) %>%
  rename(SAND = sandPerc, SILT = siltPerc, CLAY = clayPerc, plot = id) %>%
  right_join(selection, by = "plot") %>%
  as.data.frame()



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Plot ##############################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


tiff(here("outputs/figures/figure_a2_(800dpi_22x22cm).tiff"),
     res = 72, width = 22, height = 22, units = "cm", compression = "none")
TT.plot(
  class.sys = "DE.BK94.TT", 
  tri.data = sites, 
  z.name = "NtotalPerc"
  )
dev.off()
