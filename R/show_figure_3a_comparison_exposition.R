# Beta diversity on dike grasslands
# Plot Fig 3A ####
# Markus Bauer
# 2022-01-11
# Citation: 
## Bauer M, Huber J, Kollmann J (submitted) 
## Balanced turnover is a main aspect of biodiversity on restored dike grasslands: not only deterministic environmental effects, but also non-directional year and site effects drive spatial and temporal beta diversity.
## Unpublished data.



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
sites <- read_csv("data_processed_sites_temporal.csv", col_names = T, na = c("", "na", "NA"), col_types = 
                  cols(
                    .default = "?",
                    plot = "f",
                    block = "f",
                    comparison = "f",
                    locationYear = "f",
                    exposition = col_factor(levels = c("south", "north")),
                    side = col_factor(levels = c("land", "water"))
                  )) %>%
  mutate(across(c("longitude", "latitude", "riverkm", "distanceRiver"), scale)) %>%
  mutate(y = C_presence - B_presence)

### * Model ####
m3 <- blmer(y ~ comparison * exposition + PC1soil + PC2soil + PC3soil + 
              side + distanceRiver + locationYear + 
              (1|plot), 
            REML = T,
            control = lmerControl(optimizer = "Nelder_Mead"),
            cov.prior = wishart,
            data = sites)

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
    legend.position = "none",
    legend.margin = margin(0, 0, 0, 0, "cm"),
    plot.margin = margin(0, 0, 0, 0, "cm")
  )
}



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Plot ##############################################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

data_model <- ggeffect(m3, type = "emm", c("comparison", "exposition"), back.transform = T) %>%
  mutate(cross = if_else(x %in% c("1718", "1921") & group == "north", "open", "filled"),
         x = fct_recode(x, "2017 vs 2018" = "1718", "2018 vs 2019" = "1819", "2019 vs 2021" = "1921"),
         group = fct_recode(group, "North" = "north", "South" = "south"))


data <- sites %>%
  rename(predicted = y, x = comparison, group = exposition) %>%
  mutate(x = fct_recode(x, "2017 vs 2018" = "1718", "2018 vs 2019" = "1819", "2019 vs 2021" = "1921"),
         group = fct_recode(group, "North" = "north", "South" = "south"))


(graph_a <- ggplot() +
    geom_quasirandom(data = data, 
                     aes(x = x, predicted),
                     dodge.width = .6, size = 1, shape = 16, color = "grey70") + 
    geom_errorbar(data = data_model, 
                  aes(x, predicted, ymin = conf.low, ymax = conf.high), 
                  width = 0.0, size = 0.4) +
    geom_point(data = data_model,
               aes(x, predicted, shape = cross),
               size = 2) +
    facet_wrap(~group) +
    geom_hline(yintercept = 0, linetype = 2) +
    scale_y_continuous(limits = c(-.6, .5), breaks = seq(-1, 400, .2)) +
    scale_shape_manual(values = c("circle", "circle open")) +
    labs(x = "", y = expression(Gains~-~Losses~"["*italic('C')[sor]-italic('B')[sor]*"]")) +
    themeMB())

### Save ###
ggsave(here("outputs/figures/figure_3a_comparison_exposition_(800dpi_8x8cm).tiff"),
       dpi = 800, width = 8, height = 8, units = "cm")
