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
                    plot = "c",
                    locationYear = "f",
                    exposition = "f",
                    side = "c",
                    comparison = "f"
                  )) %>%
  select(-matches("PC.constructionYear"), -conf.low, -conf.high) %>%
  filter(comparison %in% c("1718", "1819", "1921") & presabu == "presence") %>%
  mutate(side = if_else(side == "water_creek", "water", side),
         ageGroup = if_else(constructionYear %in% c(2002, 2003), "0203", if_else(
           constructionYear %in% c(2006, 2007), "0607", if_else(
             constructionYear == 2008, "2008", if_else(
               constructionYear %in% c(2010, 2011), "1011", if_else(
                 constructionYear %in% c(2012, 2013), "1213", "other"
               ))))),
         ageGroup = fct_relevel(ageGroup, "2008", after = 2),
         plot = factor(plot),
         comparison = factor(comparison),
         exposition = factor(exposition),
         side = factor(side),
         ageGroup = factor(ageGroup),
         constructionYearF = factor(constructionYear)) %>%
  mutate(y = C - B)

### * Model ####
m8 <- lmer(y ~ comparison + (exposition + PC1soil + (PC2soil) + PC3soil + side) + log(distanceRiver) + locationYear + 
             comparison:exposition +
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

data_model <- ggeffect(m8, type = "emm", c("locationYear"), back.transform = T) 


data <- tbi %>%
  rename(predicted = y, x = locationYear) 


(graph_a <- ggplot() +
    #geom_quasirandom(data = data, 
    #                 aes(x = x, predicted),
    #                 dodge.width = .6, size = 1, shape = 16, color = "grey75") + 
    geom_errorbar(data = data_model, 
                  aes(x, predicted, ymin = conf.low, ymax = conf.high), 
                  width = 0.0, size = 0.4) +
    geom_point(data = data_model,
               aes(x, predicted),
               size = 2) +
    geom_hline(yintercept = 0, linetype = 2) +
    scale_y_continuous(limits = c(-.5, .5), breaks = seq(-1, 400, .1)) +
    labs(x = "", y = expression(Gains~-~Losses~"["*TBI[sor]*"]")) +
    themeMB())

### Save ###
ggsave(here("outputs/figures/figure_tbi_bc_presence_location_(800dpi_8x8cm).tiff"),
       dpi = 800, width = 8, height = 8, units = "cm")
