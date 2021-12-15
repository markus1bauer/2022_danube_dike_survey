# Show figure 2 ####
# Markus Bauer
# Citation: Markus Bauer 



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
tbi <- read_csv("data_processed_tbi.csv", col_names = T, na = c("na", "NA"), col_types = 
                  cols(
                    .default = "?",
                    id = "f",
                    locationAbb = "f",
                    block = "f",
                    plot = "f",
                    locationYear = "f",
                    exposition = "f",
                    side = "f",
                    comparison = "f"
                  )) %>%
  filter(comparison %in% c("1718", "1819", "1921")) %>%
  mutate(comparison = factor(comparison)) %>%
  pivot_wider(names_from = "presabu", values_from = "D") %>%
  group_by(plot, comparison, exposition, side, locationYear) %>%
  summarise(across(c(presence, abundance, PC1soil, PC2soil, PC3soil, distanceRiver), ~ max(.x, na.rm = T))) %>%
  ungroup()

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


ggplot(data = tbi, aes(x = presence, y = abundance)) +
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
  labs(x = expression(Dissimilarity~"["*TBI[sor]*"]"), y = expression(Dissimilarity~"["*TBI[bc]*"]")) +
  themeMB()

### Save ###
ggsave(here("outputs/figures/figure_tbi_d_d_(800dpi_8x8cm).tiff"), 
       dpi = 800, width = 8, height = 8, units = "cm")

