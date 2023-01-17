# Beta diversity on dike grasslands
# Plot Fig A8B ####
# Markus Bauer
# 2023-01-17



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Preparation ###############################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



### Packages ###
library(here)
library(tidyverse)
library(blme)
library(dotwhisker)

### Start ###
rm(list = setdiff(ls(), c("graph_a", "graph_b", "graph_c", "graph_d")))
setwd(here("data", "processed"))

### Functions ###
theme_mb <- function() {
  theme(
    panel.background = element_rect(fill = "white"),
    text = element_text(size = 9, color = "black"),
    strip.text = element_text(size = 10),
    axis.text.y = element_text(angle = 0, hjust = 1, size = 9,
                               color = "black"),
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 9,
                               color = "black"),
    axis.title.x = element_text(angle = 0, hjust = 0.5, size = 9,
                                color = "black"),
    axis.title.y = element_blank(),
    axis.line = element_line(),
    legend.key = element_rect(fill = "white"),
    legend.position = "none",
    legend.margin = margin(0, 0, 0, 0, "cm"),
    plot.margin = margin(0, 0, 0, 0, "cm")
  )
}

### Load data ###
sites <- read_csv("data_processed_sites_temporal.csv",
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
                    )) %>%
  filter(
    (comparison == "1718" | comparison == "1819" | comparison == "1921") &
      pool == "all" & presabu == "presence") %>%
  mutate(
    y = c - b,
    comparison = factor(comparison),
    river_km_scaled = scale(river_km),
    river_distance_scaled = scale(river_distance),
    biotope_distance_scaled = scale(biotope_distance),
    biotope_area_scaled = scale(biotope_area)
    )

### * Model ####
load(file = here("outputs", "models", "model_tbi_bc_all_2.Rdata"))
m <- m2
m@call



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Plot #################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



(graph_b <- m %>%
  broom.mixed::tidy(conf.int = TRUE, conf.level = .95) %>%
  filter(
    !str_detect(term, "location*") &
      !str_detect(term, "comparison") &
      !str_detect(term, "sd_*") &
      !str_detect(term, "(Intercept)")
  ) %>%
  mutate(
    term = fct_relevel(term, c(
      "river_km_scaled", "pc3_soil", "pc2_soil", "pc1_soil",
      "river_distance_scaled", "biotope_distance_scaled",
      "orientationwater", "expositionnorth"
    ))
  ) %>%
  ggplot(aes(x = estimate, y = term, xmin = conf.low, xmax = conf.high)) +
  geom_vline(xintercept = 0, linetype = 2, color = "black") +
  geom_point(size = 2, shape = "circle open") +
  geom_linerange() +
  labs(x = expression(Estimate ~ "[" * italic("C")[sor] -
                        italic("B")[sor] * "]")) +
  theme_mb() +
   theme(axis.line.y = element_blank(),
         axis.text.y = element_blank(),
         axis.ticks.y = element_blank()))

### Save ###
ggsave(here("outputs", "figures", "figure_a8b_800dpi_8x8cm.tiff"),
  dpi = 800, width = 8, height = 8, units = "cm")
