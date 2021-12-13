# Multifunctionality of dike grasslands 
# Figure 2A and 2B to figure 2
# Michaela Moosner
# 2021-11-12
# Citation: 
## Teixeira LH, Bauer M, Moosner M, Kollmann J (submitted) 
## Multifunctionality of dike grasslands: Trade-offs between flood protection, biodiversity, recreation and management. 
## unpublished data.



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Preparation ################################################################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

### Packages ###
library(here)
library(patchwork)

### Start ###
rm(list = setdiff(ls(), c("graph_a", "graph_b")))



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Plot ################################################################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


(graph_a | graph_b) +
  plot_layout(guides = 'collect', ncol = 2) +
  plot_annotation(tag_levels = "A", tag_prefix = "", tag_suffix = "") &
  theme(plot.tag = element_text(size = 10, face = "bold"),
        legend.position = "bottom")

### Save ###
ggsave(here("outputs/figures/figure_2_(800dpi_17x10cm).tiff"),
       dpi = 800, width = 17, height = 10, units = "cm")