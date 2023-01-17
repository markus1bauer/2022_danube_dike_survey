# Beta diversity on dike grasslands
# Plot Fig 4B ####
# Markus Bauer
# 2023-01-17



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Preparation ###############################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



### Packages ###
library(here)
library(tidyverse)
library(blme)
library(ggeffects)
library(ggbeeswarm)

### Start ###
rm(list = setdiff(ls(), c("graph_a", "graph_b", "graph_c", "graph_d")))
setwd(here("data", "processed"))

###  Functions ###
theme_mb <- function() {
  theme(
    panel.background = element_rect(fill = "white"),
    text  = element_text(size = 9, color = "black"),
    strip.text = element_text(size = 10),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    legend.key = element_rect(fill = "white"),
    legend.position = "none",
    legend.margin = margin(0, 0, 0, 0, "cm"),
    plot.margin = margin(0, 0, 0, 0, "cm")
  )
}

### Load data ###
sites <- read_csv("data_processed_sites_temporal.csv",
                  col_names = TRUE,
                  na = c("", "na", "NA"), col_types =
                    cols(
                      .default = "?",
                      plot = "f",
                      block = "f",
                      comparison = "f",
                      exposition = "f",
                      orientation = "f",
                      location_construction_year = "f"
                    )) %>%
  filter(
    (comparison == "1718" | comparison == "1819" | comparison == "1921") &
      pool == "all" & presabu == "presence") %>%
  mutate(
    y = d,
    comparison = factor(comparison),
    location_construction_year = fct_relevel(
      location_construction_year, "HOF-2012", after = Inf
    ),
    x = location_construction_year,
    river_km_scaled = scale(river_km),
    river_distance_scaled = scale(river_distance),
    biotope_distance_scaled = scale(biotope_distance),
    biotope_area_scaled = scale(biotope_area)
  )



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Plot ######################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

(graph_b <- ggplot() +
    geom_hline(
      yintercept = c(mean(sites$y),
                     mean(sites$y) + 0.5 * sd(sites$y),
                     mean(sites$y) - 0.5 * sd(sites$y)),
      linetype = c(1, 2, 2),
      color = "grey70"
    ) +
    geom_quasirandom(
      data = sites,
      aes(x = x, y),
      dodge.width = .6, size = 1, shape = 16,
      color = "grey70"
    ) +
    geom_boxplot(
      data = sites,
      aes(x = x, y = y),
      fill = "transparent"
    ) +
    scale_y_continuous(limits = c(0, .83), breaks = seq(0, 400, .1)) +
    scale_shape_manual(values = c("circle", "circle open")) +
    labs(x = "", y = expression(Temporal ~ "beta" ~ diversity ~
                                  "[" * italic("D")[sor] * "]")) +
    theme_mb())

### Save ###
ggsave(here("outputs", "figures", "figure_4b_800dpi_8x8cm.tiff"),
       dpi = 800, width = 8, height = 8, units = "cm")
