# Beta diversity on dike grasslands
# Figure A4 ####
# Markus Bauer
# 2022-01-11
# Citation: 
## Bauer M, Huber J, Kollmann J (submitted) 
## Balanced turnover is a main aspect of biodiversity on restored dike grasslands: not only deterministic environmental effects, but also non-directional year and site effects drive spatial and temporal beta diversity.
## Unpublished data.



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Preparation ################################################################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


### Packages ###
library(here)
library(tidyverse)

### Start ###
rm(list = ls())
setwd(here("data/processed"))

### Load data ###
sites <- read_csv("data_processed_sites_temporal.csv", col_names = T, na = c("na", "NA"), col_types = 
                  cols(
                    .default = "?",
                    block = "f",
                    plot = "f",
                    locationYear = "f",
                    exposition = "f",
                    side = "f",
                    comparison = "f"
                  )) %>%
  select(plot, D_presence, D_abundance)

### * Functions ####
themeMB <- function(){
  theme(
    panel.background = element_rect(fill = "white"),
    strip.text = element_text(size = 10),
    axis.text = element_text(angle = 0, hjust = 0.5, size = 9, color = "black"),
    axis.title = element_text(angle = 0, hjust = 0.5, size = 9, color = "black"),
    axis.line = element_line(),
    legend.key = element_rect(fill = "white"),
    legend.position = "bottom",
    legend.margin = margin(0, 0, 0, 0, "cm"),
    plot.margin = margin(0, 0, 0, 0, "cm")
  )
}



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Plot ##############################################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


ggplot(data = sites, aes(x = D_abundance, y = D_presence)) +
  stat_density_2d(aes(fill = after_stat(level)),
                  geom = "polygon",
                  contour = T,
                  bins = 8,
                  contour_var = "ndensity",
                  show.legend = F) +
  geom_point(size = 1) +
  coord_fixed(ratio = 1, xlim = c(0, 1), ylim = c(0, 1)) +
  scale_x_continuous(breaks = seq(-100, 100, 0.2)) +
  scale_y_continuous(breaks = seq(-100, 100, 0.2)) +
  scale_fill_distiller(palette = "Blues", direction = 1) +
  labs(x = expression(Temporal~"beta"~diversity~"["*italic('D')[bc]*"]"), y = expression(Temporal~"beta"~diversity~"["*italic('D')[sor]*"]")) +
  themeMB()

### Save ###
ggsave(here("outputs/figures/figure_a4_(800dpi_8x8cm).tiff"), 
       dpi = 800, width = 8, height = 8, units = "cm")

