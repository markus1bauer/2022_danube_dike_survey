# Beta diversity on dike grasslands
# Figure A1 ####
# Markus Bauer
# 2022-01-11



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Preparation ##########################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


### Packages ###
library(here)
library(tidyverse)
library(soiltexture)

### Start ###
rm(list = ls())


### Load data ###
selection <- read_csv(here("data", "processed",
                           "data_processed_sites_spatial.csv"),
  col_names = TRUE, na = "na", col_types =
    cols(
      .default = "?"
    )
) %>%
  filter(accumulatedCov > 0) %>%
  add_count(plot) %>%
  filter(surveyYear == 2017 & n == max(n)) %>%
  select(plot)
sites <- read_csv(here("data", "raw", "data_raw_sites.csv"),
  col_names = TRUE,
  na = "na", col_types =
    cols(
      .default = "?",
      vegetationCov_2017 = "d"
    )
) %>%
  select(id, sandPerc, siltPerc, clayPerc, NtotalPerc) %>%
  mutate(id = str_extract(id, "\\d\\d")) %>%
  rename(SAND = sandPerc, SILT = siltPerc, CLAY = clayPerc, plot = id) %>%
  right_join(selection, by = "plot") %>%
  as.data.frame()



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Plot #################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


tiff(here("outputs", "figures", "figure_a1_800dpi_22x22cm.tiff"),
  res = 72, width = 22, height = 22, units = "cm", compression = "none"
)
TT.plot(
  class.sys = "DE.BK94.TT",
  tri.data = sites,
  # z.name = "NtotalPerc",
  col = "black"
)
dev.off()
