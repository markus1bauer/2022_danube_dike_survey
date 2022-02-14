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
rm(list = ls())
setwd(here("data/processed"))


### Load data ###
tbi <- read_csv("data_processed_tbi.csv", col_names = T, na = c("na", "NA"), col_types = 
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
  select(-matches("PC.constructionYear"), -conf.low, -conf.high, -B, -C) %>%
  filter(comparison %in% c("1718", "1819", "1921") & presabu == "presence") %>%
  mutate(comparison = factor(comparison),
         locationYear = factor(locationYear),
         block = factor(block),
         plot = factor(plot),
         exposition = factor(exposition)
  ) %>%
  rename(y = D)

### * Model ####
m5 <- lmer(log(y) ~ comparison + exposition * (PC2soil) + PC1soil + PC3soil + side + log(distanceRiver) + locationYear + 
             (1|plot), 
           REML = T,
           data = tbi)

### * Functions ####
themeMB <- function(){
  theme(
    panel.background = element_rect(fill = "white"),
    text  = element_text(size = 9, color = "black"),
    strip.text = element_text(size = 10),
    axis.text.y = element_text(angle = 0, hjust = 0.5, size = 9, color = "black"),
    axis.text.x = element_text(angle = 90, hjust = 0.5, size = 9, color = "black"),
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

data_model <- ggeffect(m5, type = "emm", c("locationYear"), back.transform = T) %>%
  mutate(predicted = exp(predicted),
         conf.low = exp(conf.low),
         conf.high = exp(conf.high))

data <- tbi %>%
  rename(predicted = y, x = locationYear) 

(graph_a <- ggplot() +
    #geom_quasirandom(data = data, 
    #                 aes(x = x, predicted),
    #                 dodge.width = .6, size = 1, shape = 16, color = "grey75") + 
    geom_hline(yintercept = c(mean(tbi$y), mean(tbi$y) + 0.5 * sd(tbi$y), mean(tbi$y) - 0.5 * sd(tbi$y)), 
               linetype = c(1, 2, 2),
               color = "grey50") +
    geom_errorbar(data = data_model, 
                  aes(x, predicted, ymin = conf.low, ymax = conf.high), 
                  width = 0.0, size = 0.4) +
    geom_point(data = data_model,
               aes(x, predicted),
               size = 2) +
    annotate("text", 
             label = expression(italic(Chi)^2*"(11) ="~13.4*","~italic(p)~"="~2.7%*%10^-1), 
             x = 9, 
             y = 1,
             size = 3) +
    scale_y_continuous(limits = c(-.15, 1), breaks = seq(0, 400, .1)) +
    labs(x = "", y = expression(Dissimilarity~"["*TBI[sor]*"]")) +
    themeMB())

### Save ###
ggsave(here("outputs/figures/figure_tbi_d_presence_location_(800dpi_8x8cm).tiff"),
       dpi = 800, width = 8, height = 8, units = "cm")
