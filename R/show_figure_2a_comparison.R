# Show figure 2 ####
# Markus Bauer
# Citation: Markus Bauer 



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Preparation ################################################################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


### Packages ###
library(here)
library(tidyverse)
library(lme4)
library(ggeffects)
library(ggbeeswarm)

### Start ###
rm(list = setdiff(ls(), c("graph_a", "graph_b", "graph_c", "graph_d")))
setwd(here("data/processed"))


### Load data ###
tbi <- read_csv("data_processed_tbi.csv", col_names = T, na = c("", "na", "NA"), col_types = 
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
  filter(comparison %in% c("1718", "1819", "1921") & presabu == "presence") %>%
  mutate(comparison = factor(comparison),
         locationYear = factor(locationYear),
         plot = factor(plot),
         exposition = factor(exposition)
  ) %>%
  rename(y = D) %>%
  mutate(across(where(is.numeric) & !y, scale))

### * Model ####
m5 <- lmer(log(y) ~ comparison + exposition * PC2soil + PC1soil + PC3soil + side + distanceRiver + locationYear + 
             (1|plot), 
           REML = T,
           data = tbi)

### * Functions ####
themeMB <- function(){
  theme(
    panel.background = element_rect(fill = "white"),
    text  = element_text(size = 9, color = "black"),
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

data_model <- ggeffect(m5, type = "emm", c("comparison"), back.transform = T) %>%
  mutate(predicted = exp(predicted),
         conf.low = exp(conf.low),
         conf.high = exp(conf.high),
         x = fct_recode(x, "2017 vs 2018" = "1718", "2018 vs 2019" = "1819", "2019 vs 2021" = "1921"))

data <- tbi %>%
  rename(predicted = y, x = comparison) %>%
  mutate(x = fct_recode(x, "2017 vs 2018" = "1718", "2018 vs 2019" = "1819", "2019 vs 2021" = "1921"))

(graph_a <- ggplot() +
    geom_quasirandom(data = data, 
                     aes(x = x, predicted),
                     dodge.width = .6, size = 1, shape = 16, color = "grey70") + 
    geom_hline(yintercept = c(data_model$predicted[1], data_model$conf.low[1], data_model$conf.high[1]), 
               linetype = c(1, 2, 2),
               color = "grey70") +
    geom_errorbar(data = data_model, 
                  aes(x, predicted, ymin = conf.low, ymax = conf.high), 
                  width = 0.0, size = 0.4) +
    geom_point(data = data_model,
               aes(x, predicted),
               size = 2) +
    annotate("text", 
             label = c("ab", "a", "b"), 
             x = c(1, 2, 3), 
             y = c(.8, .8, .8)) +
    scale_y_continuous(limits = c(0, .8), breaks = seq(-100, 400, .1)) +
    labs(x = "", y = expression(Temporal~"beta"~diversity~"["*italic('D')[sor]*"]")) +
    themeMB())

### Save ###
ggsave(here("outputs/figures/figure_2a_(800dpi_8x8cm).tiff"),
       dpi = 800, width = 8, height = 8, units = "cm")
