# Beta diversity on dike grasslands
# Plot Fig A9 completely ####
# Markus Bauer
# 2022-09-08



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Preparation ###############################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



### Packages ###
library(here)
library(tidyverse)
library(patchwork)

### Start ###
rm(list = setdiff(ls(), c("graph_a", "graph_b", "graph_c", "graph_d")))



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Plot #######################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



((graph_a + theme(legend.position = "none")) |
   graph_b) /
  ((graph_c + theme(legend.position = c(.19, .88))) |
     (graph_d + theme(legend.position = "none"))) +
  plot_layout(guides = "keep") +
  plot_annotation(tag_levels = "A", tag_prefix = "", tag_suffix = "") &
  theme(plot.tag = element_text(size = 10, face = "bold"))

### Save ###
ggsave(
  here("outputs", "figures", "figure_a9_800dpi_17x17cm.tiff"),
  dpi = 800, width = 17, height = 17, units = "cm"
)
