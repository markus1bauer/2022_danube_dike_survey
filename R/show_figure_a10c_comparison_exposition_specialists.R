# Beta diversity on dike grasslands
# Plot Fig A10C ####

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

### Functions ###
theme_mb <- function() {
  theme(
    panel.background = element_rect(fill = "white"),
    text = element_text(size = 9, color = "black"),
    strip.text = element_text(size = 10),
    axis.text.y = element_text(angle = 0, hjust = 0.5, size = 9,
                               color = "black"),
    axis.text.x = element_text(angle = 90, hjust = 0.5, size = 9,
                               color = "black"),
    axis.title = element_text(angle = 0, hjust = 0.5, size = 9,
                              color = "black"),
    axis.line = element_line(),
    legend.key = element_rect(fill = "white"),
    legend.position = "none",
    legend.margin = margin(0, 0, 0, 0, "cm"),
    plot.margin = margin(0, 0, 0, 0, "cm")
  )
}

### Load data ###
sites <- read_csv(
  here("data", "processed", "data_processed_sites_temporal.csv"),
  col_names = TRUE, na = c("", "na", "NA"),
  col_types =
    cols(
      .default = "?",
      plot = "f",
      block = "f",
      comparison = "f",
      location = "f",
      location_construction_year = "f",
      exposition = col_factor(levels = c("south", "north")),
      orientation = col_factor(levels = c("land", "water"))
    )
) %>%
  filter(
    (comparison == "1718" | comparison == "1819" | comparison == "1921") &
      pool == "target" & presabu == "presence") %>%
  mutate(
    y = c - b,
    comparison = factor(comparison),
    location_construction_year = fct_relevel(
      location_construction_year, "HOF-2012", after = Inf
    ),
    river_km_scaled = scale(river_km),
    river_distance_scaled = scale(river_distance),
    biotope_distance_scaled = scale(biotope_distance),
    biotope_area_scaled = scale(biotope_area)
  )

### * Model ####
load(file = here("outputs", "models", "model_tbi_bc_specialist_2.Rdata"))
m <- m2
m@call



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Plot ######################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



data_model <- ggeffects::ggeffect(
  m, type = "emm", c("comparison", "exposition"), back.transform = TRUE
) %>%
  mutate(
    cross = if_else(
      x %in% c("1718", "1921") & group == "north",
      "open", "filled"
    ),
    x = fct_recode(
      x,
      "2017 vs 2018" = "1718",
      "2018 vs 2019" = "1819",
      "2019 vs 2021" = "1921"
    ),
    group = fct_recode(group, "North" = "north", "South" = "south")
  )


data <- sites %>%
  rename(predicted = y, x = comparison, group = exposition) %>%
  mutate(
    x = fct_recode(
      x,
      "2017 vs 2018" = "1718",
      "2018 vs 2019" = "1819",
      "2019 vs 2021" = "1921"
    ),
    group = fct_recode(group, "North" = "north", "South" = "south")
  )


(graph_c <- ggplot() +
    geom_quasirandom(
      data = data,
      aes(x = x, y = predicted),
      dodge.width = .6, size = 1, shape = 16, color = "grey70"
    ) +
    geom_errorbar(
      data = data_model,
      aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high),
      width = 0.0, size = 0.4
    ) +
    geom_point(
      data = data_model,
      aes(x = x, y = predicted, shape = cross),
      size = 2
    ) +
    geom_hline(
      yintercept = 0, linetype = 2,  color = "grey70"
    ) +
    facet_wrap(~group) +
    scale_y_continuous(limits = c(-.6, .55), breaks = seq(-1, 400, .2)) +
    scale_shape_manual(values = c("circle", "circle open")) +
    guides(shape = "none") +
    labs(x = "", shape = "", color = "", group = "",
         y = expression(Gains ~ -~Losses ~
                          "[" * italic("C")[sor] - italic("B")[sor] * "]")) +
    theme_mb())

### Save ###
# ggsave(
#   here("outputs", "figures",
#        "figure_a10c_comparison_exposition_800dpi_8x8cm.tiff"),
#        dpi = 800, width = 8, height = 8, units = "cm"
#   )
