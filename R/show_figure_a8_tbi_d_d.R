# Beta diversity on dike grasslands
# Figure A7 ####
# Markus Bauer
# 2022-01-11



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Preparation ###############################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


### Packages ###
library(here)
library(tidyverse)

### Start ###
rm(list = ls())
setwd(here("data", "processed"))

### Functions ###
theme_mb <- function() {
  theme(
    panel.background = element_rect(fill = "white"),
    strip.text = element_text(size = 10),
    axis.text = element_text(angle = 0, hjust = 0.5, size = 9,
                             color = "black"),
    axis.title = element_text(angle = 0, hjust = 0.5, size = 9,
                              color = "black"),
    axis.line = element_line(),
    legend.key = element_rect(fill = "white"),
    legend.position = "bottom",
    legend.margin = margin(0, 0, 0, 0, "cm"),
    plot.margin = margin(0, 0, 0, 0, "cm")
  )
}

### Load data ###
sites <- read_csv("data_processed_sites_temporal.csv",
  col_names = TRUE,
  na = c("na", "NA"), col_types =
    cols(
      .default = "?",
      block = "f",
      plot = "f",
      location_construction_year = "f",
      exposition = "f",
      orientation = "f",
      comparison = "f"
    )) %>%
  select(plot, d, comparison, presabu, pool) %>%
  filter(
    (comparison == "1718" | comparison == "1819" | comparison == "1921") &
      pool == "all") %>%
  mutate(comparison = factor(comparison)) %>%
  pivot_wider(names_from = "presabu", values_from = "d")




#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Plot ######################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


ggplot(
  data = sites,
  aes(x = abundance, y = presence)) +
  stat_density_2d(aes(fill = after_stat(level)),
    geom = "polygon",
    contour = TRUE,
    bins = 8,
    contour_var = "ndensity",
    show.legend = FALSE
  ) +
  geom_point(size = 1) +
  coord_fixed(ratio = 1, xlim = c(0, 1), ylim = c(0, 1)) +
  scale_x_continuous(breaks = seq(-100, 100, 0.2)) +
  scale_y_continuous(breaks = seq(-100, 100, 0.2)) +
  scale_fill_distiller(palette = "Blues", direction = 1) +
  labs(
    x = expression(Temporal ~ "beta" ~ diversity ~ "[" * italic("D")[bc] * "]"),
    y = expression(Temporal ~ "beta" ~ diversity ~ "[" * italic("D")[sor] * "]")
  ) +
  theme_mb()

### Save ###
ggsave(here("outputs", "figures", "figure_a8_800dpi_8x8cm.tiff"),
  dpi = 800, width = 8, height = 8, units = "cm")
