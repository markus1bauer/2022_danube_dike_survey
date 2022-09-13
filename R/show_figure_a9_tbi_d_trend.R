# Beta diversity on dike grasslands
# Plot Fig A7 ####
# Markus Bauer
# 2022-09-05



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
rm(list = ls())
setwd(here("data", "processed"))

### Functions ###
theme_mb <- function() {
  theme(
    panel.background = element_rect(fill = "white"),
    text = element_text(size = 9, color = "black"),
    strip.text = element_text(size = 10),
    axis.text = element_text(angle = 0, hjust = 0.5, size = 9, color = "black"),
    axis.title = element_text(angle = 0,
                              hjust = 0.5, size = 9, color = "black"),
    axis.line = element_line(),
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
    (comparison == "1718" | comparison == "1719" | comparison == "1721") &
      pool == "all" & presabu == "presence") %>%
  mutate(y = d,
         comparison = factor(comparison),
         across(where(is.numeric) & !y, scale),
         comparison = fct_recode(comparison, ""))



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Plot ######################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


data <- sites %>%
  mutate(comparison = fct_recode(comparison,
                        "2017 vs 2018" = "1718",
                        "2017 vs 2019" = "1719",
                        "2017 vs 2021" = "1721"))

(graph_a <- ggplot(
  data = data,
  aes(x = comparison, y = y)
  ) +
  geom_quasirandom(
    dodge.width = .6, size = 1, shape = 16, color = "grey70"
  ) +
  geom_hline(
    yintercept = c(
      mean(sites$y),
      mean(sites$y) + 0.5 * sd(sites$y),
      mean(sites$y) - 0.5 * sd(sites$y)
    ),
    linetype = c(1, 2, 2),
    color = "grey70"
  ) +
  geom_boxplot(fill = "transparent") +
  scale_y_continuous(limits = c(0, .92), breaks = seq(-100, 400, .1)) +
  scale_shape_manual(values = c("circle", "circle open")) +
  labs(x = "",
       y = expression(Temporal ~ "beta" ~ diversity ~ "[" * italic("D")[sor] * "]")) +
  theme_mb())

### Save ###
ggsave(here("outputs", "figures", "figure_a9_800dpi_8x8cm.tiff"),
  dpi = 800, width = 8, height = 8, units = "cm")
