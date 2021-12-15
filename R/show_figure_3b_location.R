# Show figure 2 ####
# Markus Bauer
# Citation: Markus Bauer 



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Preparation ################################################################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


### Packages ###
library(here)
library(tidyverse)
library(blme)
library(ggeffects)
library(ggbeeswarm)

### Start ###
rm(list = setdiff(ls(), c("graph_a", "graph_b", "graph_c", "graph_d")))
setwd(here("data/processed"))


### Load data ###
tbi <- read_csv("data_processed_tbi.csv", col_names = T, na = c("", "na", "NA"), col_types = 
                  cols(
                    .default = "?",
                    exposition = col_factor(levels = c("south", "north")),
                    side = col_factor(levels = c("land", "water"))
                  )) %>%
  filter(comparison %in% c("1718", "1819", "1921") & presabu == "presence") %>%
  mutate(plot = factor(plot),
         block = factor(block),
         comparison = factor(comparison),
         exposition = factor(exposition),
         side = factor(side),
         constructionYear = factor(constructionYear),
         locationYear = factor(locationYear)) %>%
  mutate(across(c("longitude", "latitude", "riverkm", "distanceRiver"), scale)) %>%
  mutate(y = C - B)

### * Model ####
m3 <- blmer(y ~ comparison * exposition + PC1soil + PC2soil + PC3soil + 
              side + distanceRiver + locationYear + 
              (1|plot), 
            REML = F,
            control = lmerControl(optimizer = "Nelder_Mead"),
            cov.prior = wishart,
            data = tbi)

### * Functions ####
themeMB <- function(){
  theme(
    panel.background = element_rect(fill = "white"),
    text  = element_text(size = 9, color = "black"),
    strip.text = element_text(size = 10),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0, size = 9, color = "black"),
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



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Plot ##############################################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

data_model <- ggeffect(m3, type = "emm", c("locationYear"), back.transform = T) 


data <- tbi %>%
  rename(predicted = y, x = locationYear) 


(graph_b <- ggplot() +
    geom_quasirandom(data = data, 
                     aes(x = x, predicted),
                     dodge.width = .6, size = 1, shape = 16, color = "grey70") + 
    geom_errorbar(data = data_model, 
                  aes(x, predicted, ymin = conf.low, ymax = conf.high), 
                  width = 0.0, size = 0.4) +
    geom_point(data = data_model,
               aes(x, predicted),
               size = 2,
               shape = 1) +
    geom_hline(yintercept = 0, linetype = 2) +
    scale_y_continuous(limits = c(-.6, .5), breaks = seq(-1, 400, .1)) +
    labs(x = "", y = expression(Gains~-~Losses~"["*TBI[sor]*"]")) +
    themeMB())

### Save ###
ggsave(here("outputs/figures/figure_3b_location_(800dpi_8x8cm).tiff"),
       dpi = 800, width = 8, height = 8, units = "cm")
