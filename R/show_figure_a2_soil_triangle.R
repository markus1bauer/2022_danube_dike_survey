# Beta diversity on dike grasslands
# Figure A1 ####
# Markus Bauer
# 2022-01-11



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Preparation ###############################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



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
    )) %>%
  filter(accumulated_cover > 0) %>%
  add_count(plot) %>%
  filter(survey_year == 2017 & n == max(n)) %>%
  select(plot)
sites <- read_csv(here("data", "raw", "data_raw_sites.csv"),
  col_names = TRUE, na = "na", col_types =
    cols(
      .default = "?",
      vegetation_cover.2017 = "d"
    )) %>%
  select(id, sand, silt, clay, n_total_ratio) %>%
  mutate(id = str_extract(id, "\\d\\d")) %>%
  rename(SAND = sand, SILT = silt, CLAY = clay, plot = id) %>%
  right_join(selection, by = "plot") %>%
  as.data.frame()



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Plot ######################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



tiff(here("outputs", "figures", "figure_a2_800dpi_22x22cm.tiff"),
  res = 72, width = 22, height = 22, units = "cm", compression = "none")
TT.plot(
  class.sys = "DE.BK94.TT",
  tri.data = sites,
  #z.name = "n_total_ratio",
  col = "black"
  )
dev.off()
