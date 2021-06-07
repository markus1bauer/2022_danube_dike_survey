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
                       .default = "?",
                       id = "f",
                       locationAbb = "f",
                       block = "f",
                       plot = "f",
                       exposition = col_factor(c("north", "south", "west", "east")),
                       side = "f",
                       ffh = "f",
                       changeType = col_factor(c("FFH6510", "any-FFH", "better", "change", "worse", "non-FFH"))
                     )) %>%
  select(cwmAbuHeight, id, surveyYear, constructionYear, plotAge, locationAbb, block, plot, side, exposition, PC1, PC2, PC3, changeType, conf.low, conf.high) %>%
  mutate(n = cwmAbuHeight) %>%
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
  select(-count) %>%
  mutate(exposition = factor(exposition)) %>%
  mutate(locationAbb = factor(locationAbb)) %>%
  mutate(locationAbb = factor(locationAbb, levels = unique(locationAbb[order(constructionYear)]))) %>%
  mutate(plot = factor(plot)) %>%
  mutate(block = factor(block))
sites %>%
  group_by(changeType) %>%
  summarise(mean = mean(n), sd = sd(n)) %>%
  mutate(min = mean - sd) %>%
  mutate(max = mean + sd)


### Choosen model ###
m5 <- lmer((n) ~ (exposition + side + PC1 + PC2 + PC3 + surveyYearF + locationAbb) +
             locationAbb:exposition + locationAbb:surveyYearF +
             (1|plot), sites, REML = F);

### Functions ###
themeMB <- function(){
  theme(
    panel.background = element_rect(fill = "white"),
    text  = element_text(size = 8, color = "black"),
    strip.text = element_text(size = 10),
    axis.text.y = element_text(angle = 0, hjust = 0.5),
    axis.text.x = element_text(angle = 90, hjust = 0.5),
    axis.line = element_line(),
    legend.key = element_rect(fill = "white"),
    legend.position = "bottom",
    legend.margin = margin(0, 0, 0, 0, "cm"),
    plot.margin = margin(0, 0, 0, 0, "cm")
  )
}



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Plot ################################################################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


data2 <- rename(sites, predicted = cwmAbuHeight, x = locationAbb, group = exposition)
predicted <- c(0.321, 0.340)
conf.low <- c(0.225, 0.230)
conf.high <- c(0.418, 0.451)
group <- c("ffh6510", "non-ffh")
meandata <- tibble(predicted, conf.low, conf.high, group)
rm(predicted, conf.low, conf.high, group)
pd <- position_dodge(.6)
(graph <- ggplot(data2, 
                 aes(x, predicted, ymin = conf.low, ymax = conf.high)) +
    geom_hline(aes(yintercept = predicted, color = group), meandata, 
               size = .25, show.legend = F) +
    geom_hline(aes(yintercept = conf.low, color = group), meandata, 
               linetype = "dashed", size = .25, show.legend = F) +
    geom_hline(aes(yintercept = conf.high, color = group), meandata, 
               linetype = "dashed", size = .25, show.legend = F) +
    geom_boxplot(data = data2, aes(x, predicted)) +
    geom_quasirandom(data = data2, aes(x, predicted), 
                     color = "grey70", dodge.width = .6, size = 0.7) +
    #facet_grid(~ group) +
    annotate("text", label = expression(italic(p)==4.6%*%10^-2), x = 14, y = 1.1, size = 2) +
    #annotate("text", label = "mean FFH6510 plots", x = 2.4, y = 30, hjust = T, size = 2) +
    #annotate("text", label = "mean non-FFH plots", x = 2.4, y = 17.5, hjust = T, size = 2) +
    scale_y_continuous(limits = c(0, 1.1), breaks = seq(-100, 100, .5)) +
    scale_color_manual(values = c("grey70", "grey90")) +
    #scale_fill_manual(values = c("cyan3", "lightcoral")) +
    labs(x = "Location and restoration year", fill = "Exposition: ", y = expression(Delta~CWM~canopy~height~"["*m*"]")) +
    themeMB()
)

### Save ###
setwd(here("outputs/figures"))
ggsave("figure_cwmAbuHeight_location_(800dpi_16x7cm).tiff",
       dpi = 800, width = 16, height = 7, units = "cm")
