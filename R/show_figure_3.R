# Beta diversity on dike grasslands
# Plot Fig 3 ####
# Markus Bauer
# 2022-01-11



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Preparation ##########################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


### Packages ###
library(here)
library(tidyverse)
library(patchwork)

### Start ###
rm(list = setdiff(ls(), c("graph_a", "graph_b", "graph_c", "graph_d")))



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Plot #################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


((graph_a + theme(legend.position = "none")) | graph_b) /
  (plot_spacer()  | graph_c) +
  plot_layout(heights = c(5, 4), guides = 'keep') +
  plot_annotation(tag_levels = "A", tag_prefix = "", tag_suffix = "") &
  theme(plot.tag = element_text(size = 10, face = "bold"))

### Save ###
ggsave(here("outputs", "figures", "figure_3_800dpi_17x17cm.tiff"),
       dpi = 800, width = 17, height = 17, units = "cm")
