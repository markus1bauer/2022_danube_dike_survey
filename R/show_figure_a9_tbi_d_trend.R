# Beta diversity on dike grasslands
# Plot Fig A9 ####

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
rm(list = ls())

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
sites <- read_csv(
  here("data", "processed", "data_processed_sites_temporal.csv"),
  col_names = TRUE, na = c("", "na", "NA"), col_types = cols(
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
  mutate(
    y = d,
    comparison = factor(comparison),
    river_km_scaled = scale(river_km),
    river_distance_scaled = scale(river_distance),
    biotope_distance_scaled = scale(biotope_distance),
    biotope_area_scaled = scale(biotope_area)
    )

### Model ###
m <- blmer(
  log(y) ~ comparison + exposition + pc1_soil + pc2_soil + pc3_soil +
    orientation + river_distance_scaled + biotope_distance_scaled +
    river_km_scaled +
    (1 | plot),
  REML = FALSE,
  control = lmerControl(optimizer = "Nelder_Mead"),
  cov.prior = wishart,
  data = sites
)
DHARMa::simulateResiduals(m, plot = TRUE)
MuMIn::r.squaredGLMM(m)



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Plot ######################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



data_model <- ggeffect(m, type = "emm", c("comparison"),
                       back.transform = TRUE) %>%
  mutate(
    predicted = exp(predicted),
    conf.low = exp(conf.low),
    conf.high = exp(conf.high),
    x = fct_recode(
      x,
      "2017 vs 2018" = "1718",
      "2017 vs 2019" = "1719",
      "2017 vs 2021" = "1721"
    )
  )

data <- sites %>%
  rename(predicted = y, x = comparison) %>%
  mutate(
    x = fct_recode(
      x,
      "2017 vs 2018" = "1718",
      "2017 vs 2019" = "1719",
      "2017 vs 2021" = "1721"
    )
  )

(graph_a <- ggplot(
  data = data,
  aes(x = x, y = predicted)
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
    geom_errorbar(
      data = data_model,
      aes(x, predicted, ymin = conf.low, ymax = conf.high),
      width = 0.0, size = 0.4
    ) +
    geom_point(
      data = data_model,
      aes(x, predicted),
      size = 2
    ) +
    scale_y_continuous(limits = c(0, .92), breaks = seq(-100, 400, .1)) +
    scale_shape_manual(values = c("circle", "circle open")) +
    labs(x = "",
         y = expression(Temporal ~ "beta" ~ diversity ~
                          "[" * italic("D")[sor] * "]")) +
    theme_mb())

### Save ###
# ggsave(
#   here("outputs", "figures", "figure_a9_800dpi_8x8cm.tiff"),
#   dpi = 800, width = 8, height = 8, units = "cm"
#   )
