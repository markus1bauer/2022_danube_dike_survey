# Show Figure targetRichness ####
# Markus Bauer



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Preparation ################################################################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


### Packages ###
library(here)
library(tidyverse)
library(ggbeeswarm)
library(lme4)
library(emmeans)
library(ggeffects)

### Start ###
rm(list = ls())
setwd(here("data/processed"))

### Load data ###
sites <- read_csv2("data_processed_sites.csv", col_names = T, na = c("na", "NA"), col_types = 
                     cols(
                       .default = col_guess(),
                       id = "f",
                       location = "f",
                       block = "f",
                       plot = "f",
                       exposition = "f",
                       side = "f",
                       ffh = "f",
                       changeType = col_factor(c("FFH6510", "any-FFH", "better", "change", "worse", "non-FFH"))
                     )) %>%
  select(targetRichness, id, surveyYear, constructionYear, plotAge, location, block, plot, side, exposition, PC1, PC2, PC3, conf.low, conf.high) %>%
  mutate(n = targetRichness) %>%
  mutate(location = factor(location, levels = unique(location[order(constructionYear)]))) %>%
  mutate(plotAge = scale(plotAge, scale = T, center = F)) %>%
  mutate(surveyYear = scale(surveyYear, scale = T, center = F)) %>%
  mutate(surveyYearF = as_factor(surveyYear)) %>%
  mutate(constructionYear = scale(constructionYear, scale = T, center = F)) %>%
  mutate(constructionYearF = as_factor(constructionYear)) %>%
  filter(exposition != "east",
         exposition != "west") %>%
  group_by(plot) %>%
  mutate(count = n()) %>%
  ungroup() %>%
  slice_max(count, n = 1) %>%
  mutate(exposition = factor(exposition)) %>%
  mutate(location = factor(location)) %>%
  mutate(plot = factor(plot)) %>%
  mutate(block = factor(block)) %>%
  select(-count)

### Choosen model ###
m3 <- lmer((n) ~ (exposition + side + PC1 + PC2 + PC3 + surveyYearF + location) + 
             (1|plot), sites, REML = F);

### Functions ###
themeMB <- function(){
  theme(
    panel.background = element_rect(fill = "white"),
    text  = element_text(size = 8, color = "black"),
    strip.text = element_text(size = 10),
    axis.text.y = element_text(angle = 0, hjust = 0.5),
    axis.line = element_line(),
    legend.key = element_rect(fill = "white"),
    legend.position = "none",
    legend.margin = margin(0, 0, 0, 0, "cm"),
    plot.margin = margin(0, 0, 0, 0, "cm")
  )
}



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Plot ################################################################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


pdata <- ggemmeans(m3, terms = c("PC2"), type = "fe")
sites2 <- rename(sites, predicted = targetRichness, x = PC2, group = surveyYearF)
predicted <- c(29.1, 16.6)
conf.low <- c(24.3, 11.1)
conf.high <- c(33.9, 22)
group <- c("ffh6510", "non-ffh")
meandata <- tibble(predicted, conf.low, conf.high, group)
rm(predicted, conf.low, conf.high, group)
pd <- position_dodge(.6)
(graph <- ggplot(pdata, 
                 aes(x, predicted, ymin = conf.low, ymax = conf.high)) +
    geom_point(data = sites2, aes(x, predicted), 
               color = "grey70", size = 0.7) +
    geom_hline(aes(yintercept = predicted, color = group), meandata, 
               size = .25) +
    geom_hline(aes(yintercept = conf.low, color = group), meandata, 
               linetype = "dashed", size = .25) +
    geom_hline(aes(yintercept = conf.high, color = group), meandata, 
               linetype = "dashed", size = .25) +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .3) +
    geom_line() +
    annotate("text", label = expression(italic(p)==3.2%*%10^-5), x = .5, y = 40, size = 2) +
    annotate("text", label = "mean FFH6510 plots", x = -1.2, y = 30, hjust = T, size = 2) +
    annotate("text", label = "mean non-FFH plots", x = -1.2, y = 17.5, hjust = T, size = 2) +
    annotate("text", label = "Phosphorous", x = -1.2, y = 1, hjust = T, size = 2) +
    annotate("text", label = "Calcium carbonate", x = .7, y = 1, hjust = T, size = 2) +
    scale_y_continuous(limits = c(0, 40), breaks = seq(-100, 100, 10)) +
    scale_x_continuous(limits = c(-2, 1), breaks = seq(-100, 100, .5)) +
    #scale_shape_manual(values = c(1,16)) +
    labs(x = "PC2", y = expression(paste("Target species richness [#]"))) +
    scale_color_manual(values = c("grey70", "grey90")) +
    themeMB()
)

### Save ###
setwd(here("outputs/figures"))
ggsave("figure_targetRichness_pc2_(800dpi_8x7cm).tiff",
       dpi = 800, width = 8, height = 7, units = "cm")
