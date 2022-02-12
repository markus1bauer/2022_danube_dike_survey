# Beta diversity on dike grasslands
# Plot Fig A5 ####
# Markus Bauer
# 2022-01-11



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Preparation #########################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


### Packages ###
library(here)
library(tidyverse)
library(blme)
library(ggeffects)
library(ggbeeswarm)

### Start ###
rm(list = setdiff(ls(), c("graph_a", "graph_b", "graph_c", "graph_d")))
setwd(here("data", "processed"))

### Load data ###
sites <- read_csv("data_processed_sites_temporal.csv", col_names = TRUE,
                  na = c("", "na", "NA"), col_types =
                    cols(
                      .default = "?",
                      plot = "f",
                      block = "f",
                      comparison = "f",
                      exposition = "f",
                      side = "f",
                      locationYear = "f"
                    )) %>%
  rename(y = D_presence) %>%
  filter(comparison == "1718" | comparison == "1719" | comparison == "1721") %>%
  mutate(across(where(is.numeric) & !y, scale))

### * Model ####
m2 <- blmer(log(y) ~ comparison + exposition * PC1soil + PC2soil + PC3soil + 
              side + distanceRiver + locationYear + abundance +
              (1|plot), 
            REML = T,
            control = lmerControl(optimizer = "Nelder_Mead"),
            cov.prior = wishart,
            data = sites)

### * Functions ####
themeMB <- function() {
  theme(
    panel.background = element_rect(fill = "white"),
    text  = element_text(size = 9, color = "black"),
    strip.text = element_text(size = 10),
    axis.text = element_text(angle = 0, hjust = 0.5, size = 9, color = "black"),
    axis.title = element_text(angle = 0, hjust = 0.5, size = 9, color = "black"),
    axis.line = element_line(),
    legend.key = element_rect(fill = "white"),
    legend.position = "none",
    legend.margin = margin(0, 0, 0, 0, "cm"),
    plot.margin = margin(0, 0, 0, 0, "cm")
  )
}


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Plot ################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

data_model <- ggeffect(m2, type = "emm", c("comparison"), back.transform = TRUE) %>%
  mutate(predicted = exp(predicted),
         conf.low = exp(conf.low),
         conf.high = exp(conf.high),
         cross = if_else(x %in% c("1819"), "filled", "open"),
         x = fct_recode(x, "2017 vs 2018" = "1718", "2017 vs 2019" = "1719", "2017 vs 2021" = "1721"))

data <- sites %>%
  rename(predicted = y, x = comparison) %>%
  mutate(x = fct_recode(x, "2017 vs 2018" = "1718", "2017 vs 2019" = "1719", "2017 vs 2021" = "1721"))

(graph_a <- ggplot() +
    geom_quasirandom(data = data, 
                     aes(x = x, predicted),
                     dodge.width = .6, size = 1, shape = 16, color = "grey70") + 
    geom_hline(yintercept = c(mean(sites$y),
                              mean(sites$y) + 0.5 * sd(sites$y),
                              mean(sites$y) - 0.5 * sd(sites$y)), 
               linetype = c(1, 2, 2),
               color = "grey70") +
    geom_errorbar(data = data_model, 
                  aes(x, predicted, ymin = conf.low, ymax = conf.high), 
                  width = 0.0, size = 0.4) +
    geom_point(data = data_model,
               aes(x, predicted, shape = cross),
               size = 2) +
    scale_y_continuous(limits = c(0, .92), breaks = seq(-100, 400, .1)) +
    scale_shape_manual(values = c("circle", "circle open")) +
    labs(x = "", y = expression(Temporal~"beta"~diversity~"["*italic('D')[sor]*"]")) +
    themeMB())

### Save ###
ggsave(here("outputs", "figures", "figure_a5_800dpi_8x8cm.tiff"),
       dpi = 800, width = 8, height = 8, units = "cm")
