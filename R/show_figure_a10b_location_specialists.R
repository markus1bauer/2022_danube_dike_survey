# Beta diversity on dike grasslands
# Plot Fig A10B ####
# Markus Bauer
# 2022-09-14



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
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
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
      pool == "target" & presabu == "presence") %>%
  mutate(y = d,
         comparison = factor(comparison),
         location_construction_year = fct_reorder(
           location_construction_year, construction_year
         ),
         across(c("river_km", "river_distance"), scale))

### * Model ####
m5 <- blmer(
  log(y) ~ comparison + exposition + pc1_soil + pc2_soil + pc3_soil +
    orientation + river_distance + location_construction_year +
    (1 | plot),
  REML = FALSE,
  control = lmerControl(optimizer = "Nelder_Mead"),
  cov.prior = wishart,
  data = sites
)



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Plot ######################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



data_model <- ggeffect(m5, type = "emm", c("location_construction_year"),
                       back.transform = TRUE) %>%
  mutate(predicted = exp(predicted),
         conf.low = exp(conf.low),
         conf.high = exp(conf.high),
         cross = if_else(
           x %in% c("STO-2002", "IRL-2003", "OEB-2006", "PFE-2008"),
           "filled", "open"
         ))

data <- sites %>%
  rename(predicted = y, x = location_construction_year)

(graph_b <- ggplot() +
    geom_quasirandom(
      data = data,
      aes(x = x, predicted),
      dodge.width = .6, size = 1, shape = 16,
      color = "grey70"
    ) +
    geom_hline(
      yintercept = c(mean(sites$y),
                     mean(sites$y) + 0.5 * sd(sites$y),
                     mean(sites$y) - 0.5 * sd(sites$y)),
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
      aes(x, predicted, shape = cross),
      size = 2
    ) +
    scale_y_continuous(limits = c(0, .84), breaks = seq(0, 400, .1)) +
    scale_shape_manual(values = c("circle", "circle open")) +
    labs(x = "", y = expression(Temporal ~ "beta" ~ diversity ~
                                  "[" * italic("D")[sor] * "]")) +
    theme_mb())

### Save ###
ggsave(here("outputs", "figures", "figure_a10b_800dpi_8x8cm.tiff"),
       dpi = 800, width = 8, height = 8, units = "cm")
