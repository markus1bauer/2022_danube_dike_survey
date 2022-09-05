# Beta diversity on dike grasslands
# Plot Fig 4D ####
# Markus Bauer
# 2022-09-05



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Preparation ################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



### Packages ###
library(here)
library(tidyverse)
library(blme)
library(ggeffects)
library(ggbeeswarm)

### Start ###
rm(list = setdiff(ls(), c("graph_a", "graph_b", "graph_c", "graph_d")))
setwd(here("data", "processed"))

### Functions ####
theme_mb <- function() {
  theme(
    panel.background = element_rect(fill = "white"),
    text = element_text(size = 9, color = "black"),
    strip.text = element_text(size = 10),
    axis.text.x = element_text(
      angle = 90, vjust = 0.5, hjust = 0,
      size = 9, color = "black"
    ),
    axis.line.x = element_line(),
    axis.line.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
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
  mutate(y = c - b,
         comparison = factor(comparison)) %>%
  mutate(across(c("river_km", "river_distance"), scale))

### * Model ####
m3 <- blmer(
  y ~ comparison * exposition + pc1_soil + pc2_soil + pc3_soil +
    orientation + river_distance + location_construction_year +
    (1 | plot),
  REML = TRUE,
  control = lmerControl(optimizer = "Nelder_Mead"),
  cov.prior = wishart,
  data = sites
)




#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Plot #######################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



data_model <- ggeffect(m3, type = "emm", c("location_construction_year"),
                       back.transform = TRUE)


data <- sites %>%
  rename(predicted = y, x = locationYear)


(graph_b <- ggplot() +
  geom_quasirandom(
    data = data,
    aes(x = x, predicted),
    dodge.width = .6, size = 1, shape = 16, color = "grey70"
  ) +
  geom_errorbar(
    data = data_model,
    aes(x, predicted, ymin = conf.low, ymax = conf.high),
    width = 0.0, size = 0.4
  ) +
  geom_point(
    data = data_model,
    aes(x, predicted),
    size = 2,
    shape = 1
  ) +
  geom_hline(yintercept = 0, linetype = 2) +
  scale_y_continuous(limits = c(-.6, .5), breaks = seq(-1, 400, .1)) +
  labs(x = "", y = expression(Gains ~ -~Losses ~ "[" * TBI[sor] * "]")) +
  theme_mb())

### Save ###
ggsave(
  here("outputs", "figures", "figure_3b_location_800dpi_8x8cm.tiff"),
  dpi = 800, width = 8, height = 8, units = "cm"
  )
