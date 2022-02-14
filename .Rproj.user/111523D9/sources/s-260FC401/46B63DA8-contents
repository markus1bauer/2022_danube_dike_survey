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

data_model <- ggeffect(m5, type = "emm", c("PC2soil", "exposition"), back.transform = T) %>%
  mutate(predicted = exp(predicted),
         conf.low = exp(conf.low),
         conf.high = exp(conf.high),
         group = fct_recode(group, "North" = "north", "South" = "south"))

data <- tbi %>%
  rename(predicted = y, x = PC2soil, group = exposition) %>%
  mutate(group = fct_recode(group, "North" = "north", "South" = "south"))

(graph_a <- ggplot() +
    geom_ribbon(data = data_model,
                aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, fill = group),
                alpha = .3) +
    geom_point(data = data,
               aes(x = x, y = predicted, color = group), 
               size = 1) +
    geom_line(data = data_model,
              aes(x = x, y = predicted, color = group)) +
    annotate("text", 
             label = expression(italic(p)~"="~8.4%*%10^-2), 
             x = .5, 
             y = .8,
             size = 2.5) +
    scale_y_continuous(limits = c(0, .8), breaks = seq(-100, 400, .1)) +
    scale_x_continuous(breaks = seq(-100, 400, .5)) +
    annotate("text",
             label = c("Phosphorus", "Sandy"),
             x = c(-1.2, 0.5),
             y = c(0, 0),
             size = 2.5) +
    labs(x = "PC2", y = expression(Dissimilarity~"["*TBI[sor]*"]"), color = "Exposition", fill = "Exposition") +
    themeMB())

### Save ###
ggsave(here("outputs/figures/figure_tbi_d_presence_pc2_exposition_(800dpi_8x8cm).tiff"),
       dpi = 800, width = 8, height = 8, units = "cm")

ggpredict(m5, c("PC2soil", "exposition"), alpha = .1) %>%
  plot(add.data = T, limit.range = T) +
  scale_y_continuous(limits = c(0, .95), breaks = seq(-100, 400, .1)) +
  scale_x_continuous(breaks = seq(-100, 400, .5)) +
  annotate("text",
           label = c("Phosphorus", "Sandy"),
           x = c(-1.5, 0.5),
           y = c(0, 0)) +
  labs(x = "PC2", y = expression(Dissimilarity~"["*TBI[sor]*"]"), color = "Exposition") +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 9, color = "black"),
        axis.title = element_text(size = 9, color = "black"),
        axis.line = element_line(color = "black"))