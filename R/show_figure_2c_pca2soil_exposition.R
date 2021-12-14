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
                    exposition = col_factor(levels = c("north", "south")),
                    side = "f",
                    comparison = "f"
                  )) %>%
  filter(comparison %in% c("1718", "1819", "1921") & presabu == "presence") %>%
  mutate(comparison = factor(comparison),
         locationYear = factor(locationYear),
         block = factor(block),
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
    axis.text.y = element_text(angle = 0, hjust = 0.5, size = 9, color = "black"),
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 9, color = "black"),
    axis.title = element_text(angle = 0, hjust = 0.5, size = 9, color = "black"),
    axis.line = element_line(),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    legend.key = element_rect(fill = "white"),
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

(graph_c <- ggplot() +
    geom_ribbon(data = data_model,
                aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, fill = group),
                alpha = .3) +
    geom_point(data = data,
               aes(x = x, y = predicted, shape = group), 
               size = 1, color = "grey70", fill = "grey70") +
    geom_line(data = data_model,
              aes(x = x, y = predicted, group = group)) +
    scale_y_continuous(limits = c(0, .8), breaks = seq(-100, 400, .1)) +
    scale_x_continuous(breaks = seq(-100, 400, .5)) +
    scale_fill_manual(values = c("grey40", "grey70")) +
    scale_shape_manual(values = c("circle filled", "circle open")) +
    annotate("text",
             label = c("Phosphorus", "Sand"),
             x = c(-1.6, 0.7),
             y = c(0, 0),
             size = 2.5) +
    labs(x = "PC2", y = expression(Temporal~"beta"~diversity~"["*italic('D')[sor]*"]"), color = "Exposition", fill = "Exposition", shape = "Exposition") +
    themeMB() +
    theme(legend.position = c(.8, .9)))

### Save ###
ggsave(here("outputs/figures/figure_2c_(800dpi_8x8cm).tiff"),
       dpi = 800, width = 8, height = 8, units = "cm")
