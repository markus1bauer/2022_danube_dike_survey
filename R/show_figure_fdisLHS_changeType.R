# Show Figure targetRichness ####
# Markus Bauer



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Preparation ################################################################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


### Packages ###
library(tidyverse)
library(ggbeeswarm)
library(lme4)
library(emmeans)
library(ggeffects)

### Start ###
rm(list = ls())
setwd(here("data/processed"))

### Load data ###
sites <- read_csv2("data_processed_sites.csv", col_names = T, na = "na", col_types = 
                     cols(
                       .default = col_guess(),
                       id = col_factor(),
                       location = col_factor(),
                       block = col_factor(),
                       plot = col_factor(),
                       exposition = col_factor(),
                       ffh = col_factor(),
                       changeType = col_factor(c("FFH6510", "better", "change", "worse", "any-FFH", "non-FFH"))
                     )) %>%
  select(fdisAbuLHS, surveyYear, constructionYear, id, block, location, plot, exposition, plotAge, PC1, PC2, PC3, changeType, ffh, conf.low, conf.high) %>%
  mutate(n = fdisAbuLHS) %>%
  mutate(plotAge = scale(plotAge, scale = T, center = T)) %>%
  mutate(surveyYear = scale(surveyYear, scale = T, center = T)) %>%
  mutate(surveyYearF = as_factor(surveyYear)) %>%
  mutate(constructionYearF = as_factor(constructionYear)) %>%
  filter(ffh != "6210", 
         changeType != "any-FFH", 
         !is.na(changeType),
         exposition != "east",
         exposition != "west") %>%
  mutate(changeType = fct_collapse(changeType, "change/worse" = c("change", "worse")))

### Choosen model ###
m2 <- lmer(log(n) ~ (exposition + PC1 + PC2 + PC3 + ffh + changeType) + 
             (surveyYear|location/plot), sites, REML = F)

### Functions ###
themeMB <- function(){
  theme(
    panel.background = element_rect(fill = "white"),
    text  = element_text(size = 8, color = "black"),
    strip.text = element_text(size = 10),
    axis.text.y = element_text(angle = 0, hjust = 0.5),
    axis.line.y = element_line(),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.key = element_rect(fill = "white"),
    legend.position = "right",
    legend.margin = margin(0, 0, 0, 0, "cm"),
    plot.margin = margin(0, 0, 0, 0, "cm")
  )
}



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Plot ################################################################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


pdata <- ggemmeans(m2, terms = c("changeType"), type = "fe")
data2 <- rename(sites, predicted = fdisAbuLHS, x = changeType, group = surveyYearF)
meandata <- filter(pdata, x == "FFH6510")
pd <- position_dodge(.6)
(graph <- ggplot(pdata, 
                 aes(x, predicted, ymin = conf.low, ymax = conf.high)) +
    geom_quasirandom(data = data2, aes(x, predicted), 
                     color = "grey70", dodge.width = .6, size = 0.7) +
    geom_hline(aes(yintercept = predicted), meandata, 
               color = "grey70", size = .25) +
    geom_hline(aes(yintercept = conf.low), meandata, 
               color = "grey70", linetype = "dashed", size = .25) +
    geom_hline(aes(yintercept = conf.high), meandata, 
               color = "grey70", linetype = "dashed", size = .25) +
    geom_errorbar(position = pd, width = 0.0, size = 0.4) +
    geom_point(position = pd, size = 2.5) +
    #facet_grid(~ group) +
    annotate("text", label = "n.s.", x = 4.2, y = 2) +
    scale_y_continuous(limits = c(0, 2), breaks = seq(-100, 100, .5)) +
    #scale_shape_manual(values = c(1,16)) +
    labs(x = "FFH type", y = expression(paste("FDis of LHS"))) +
    themeMB()
)

### Save ###
setwd(here("outputs/figures"))
ggsave("figure_fdisLHS_(800dpi_8x7cm).tiff",
       dpi = 800, width = 8, height = 7, units = "cm")
