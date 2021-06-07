# Show Figure cwmSla ####
# Markus Bauer



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Preparation ################################################################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


### Packages ###
library(here)
library(tidyverse)
library(lme4)
library(dotwhisker)


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
                       exposition = "f",
                       side = "f",
                       ffh = "f",
                       changeType = col_factor(c("FFH6510", "any-FFH", "better", "change", "worse", "non-FFH"))
                     )) %>%
  select(cwmAbuSla, id, surveyYear, constructionYear, plotAge, locationAbb, block, plot, side, exposition, PC1, PC2, PC3, changeType, ffh) %>%
  mutate(n = cwmAbuSla) %>%
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

m <- broom::tidy(m5) %>% 
  filter(term != "(Intercept)",
         !grepl("location", term),
         !grepl("surveyYear", term),
         !grepl("sd", term)) %>%
  mutate(group = c("Dike", "Dike", "Soil", "Soil", "Soil"))
dwplot(m,
       dot_args = list(colour = "black"),
       whisker_args = list(colour = "black"),
       vline = geom_vline(xintercept = 0, colour = "black", linetype = 2)) %>%
  relabel_predictors(c(expositionnorth = "Exposition (North vs. south)",
                       sideland = "Side (Land vs. water)")) +
  #scale_x_continuous(limits = c(-2, 4), breaks = seq(-100, 100, 2)) +
  xlab(expression(Delta~Specific~leaf~area~"["*mm^2*mg^-1*"]")) +
  annotate("text", label = expression(paste("R"^2*""[m]*"= 0.02")), x = 4, y = 5.3, hjust = T, size = 2) +
  annotate("text", label = expression(paste("R"^2*""[m2]*"= 0.41")), x = 4, y = 4.9, hjust = T, size = 2) +
  annotate("text", label = expression(paste("R"^2*""[c]*"= 0.43")), x = 4, y = 4.5, hjust = T, size = 2) +
  themeMB()

### Save ###
setwd(here("outputs/figures"))
ggsave("figure_cwmSla_(800dpi_16x5cm).tiff",
       dpi = 800, width = 16, height = 5, units = "cm")
