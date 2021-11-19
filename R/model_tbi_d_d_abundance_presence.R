# Calculate TBI ####
# Markus Bauer
# Citation: Bauer M, Huber J, Kollmann J (202x)



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Preparation ################################################################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


### Packages ###
library(here)
library(tidyverse)
library(ggbeeswarm)
library(lme4)  
library(DHARMa)
library(emmeans)

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
  filter(comparison %in% c("1718", "1819", "1921")) %>%
  mutate(comparison = factor(comparison)) %>%
  pivot_wider(names_from = "presabu", values_from = "D") %>%
  group_by(plot, comparison, exposition, side, locationYear) %>%
  summarise(across(c(presence, abundance, PC1soil, PC2soil, PC3soil, distanceRiver), ~ max(.x, na.rm = T))) %>%
  ungroup()



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Statistics ################################################################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


## 1 Data exploration #####################################################################################

### a Graphs ---------------------------------------------------------------------------------------------
#main
ggplot(tbi, aes(x = presence, y = abundance)) + 
  geom_point() +
  geom_smooth(method = "lm")
#2way
ggplot(tbi, aes(x = presence, y = abundance, color = exposition)) + 
  geom_point() +
  geom_smooth(method = "lm")
ggplot(data = tbi, aes(x = presence, y = PC1soil)) +
  geom_point() +
  geom_density_2d_filled(alpha = .5,
                         contour_var = "ndensity",
                         show.legend = F) +
  #facet_wrap(~comparison) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_continuous(limits = c(0, 1))
ggplot(data = tbi, aes(x = presence, y = abundance)) +
  geom_point() +
  geom_density_2d_filled(alpha = .5,
                         contour_var = "ndensity",
                         show.legend = F) +
  #facet_wrap(~comparison) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_continuous(limits = c(0, 1))
ggsave(here("outputs/figures/figure_tbi_d_d_presabu_(800dpi_11x10cm).tiff"), dpi = 800, width = 11, height = 10, units = "cm")
range(tbi$abundance)
range(tbi$presence)

ggplot(tbi, aes(PC2soil)) + geom_density()
ggplot(tbi, aes(exp(PC2soil))) + geom_density()

ggplot(data = tbi, aes(x = presence, y = abundance, color = PC2soil)) +
  geom_point()
m <- lmer((abundance) ~ presence + (PC2soil) + comparison + exposition + PC1soil + PC3soil + side + log(distanceRiver) + locationYear +
            (1|plot), 
          REML = T,
          data = tbi);simulateResiduals(m, plot = T)
car::Anova(m, type = 2)
sjPlot::plot_model(m, type = "emm", terms = "presence", show.data = T)
ggsave(here("outputs/figures/figure_tbi_d_d_presabu_reg_(800dpi_11x10cm).tiff"), dpi = 800, width = 11, height = 10, units = "cm")
visreg::visreg2d(m, "presence", "PC2soil", trans = log)


### b Outliers, zero-inflation, transformations? -----------------------------------------------------
dotchart((tbi$y), groups = factor(tbi$exposition), main = "Cleveland dotplot")
tbi %>% count(locationYear)
boxplot(tbi$y);#identify(rep(1, length(etbi$rgr13)), etbi$rgr13, labels = c(etbi$n))
plot(table((tbi$y)), type = "h", xlab = "Observed values", ylab = "Frequency")
ggplot(tbi, aes(y)) + geom_density()
ggplot(tbi, aes(sqrt(y))) + geom_density()


## 2 Model building ################################################################################

### a models ----------------------------------------------------------------------------------------
#random structure
m1a <- lmer(y ~ 1 + (1|locationYear), data = tbi, REML = F)
VarCorr(m1a) 
m1b <- lmer(y ~ 1 + (1|locationAbb), data = tbi, REML = F)
VarCorr(m1b)
#fixed effects
m1 <- lm((y) ~ comparison * exposition * (PC1soil + PC2soil + PC3soil) + side + plot, 
           data = tbi)
simulateResiduals(m1, plot = T)
m2 <- lm((y) ~ comparison * exposition + (PC1soil + PC2soil + PC3soil) + side + plot, 
           data = tbi)
simulateResiduals(m2, plot = T)
m3 <- lm((y) ~ comparison + exposition * (PC1soil + PC2soil + PC3soil) + side + plot, 
           data = tbi)
simulateResiduals(m3, plot = T)
m4 <- lm((y) ~ exposition + comparison * (PC1soil + PC2soil + PC3soil) + side + plot, 
           data = tbi)
simulateResiduals(m4, plot = T)
m5 <- lm((y) ~ comparison + exposition + PC1soil + PC2soil + PC3soil + side + plot,
           data = tbi)
simulateResiduals(m5, plot = T)
m6 <- lm((y) ~ comparison + exposition * PC1soil + PC2soil + PC3soil + side + plot,
         data = tbi)
simulateResiduals(m6, plot = T)
m7 <- lm((y) ~ exposition + comparison * PC1soil + PC2soil + PC3soil + side + plot,
         data = tbi)
simulateResiduals(m7, plot = T)

### b comparison -----------------------------------------------------------------------------------------
AICcmodavg::aictab(cand.set = list("m1" = m1, "m2" = m2, "m3" = m3, "m4" = m4, "m5" = m5, "m6" = m6, "m7" = m7))
rm(m1a, m1b, m1c, m2, m3, m4, m5)

### c model check -----------------------------------------------------------------------------------------
simulationOutput <- simulateResiduals(m2, plot = T)
plotResiduals(simulationOutput$scaledResiduals, tbi$locationYear)
plotResiduals(simulationOutput$scaledResiduals, tbi$plot)
plotResiduals(simulationOutput$scaledResiduals, tbi$plotAge)
plotResiduals(simulationOutput$scaledResiduals, tbi$comparison)
plotResiduals(simulationOutput$scaledResiduals, tbi$exposition)
plotResiduals(simulationOutput$scaledResiduals, tbi$side)
plotResiduals(simulationOutput$scaledResiduals, tbi$PC1soil)
plotResiduals(simulationOutput$scaledResiduals, tbi$PC2soil)
plotResiduals(simulationOutput$scaledResiduals, tbi$PC3soil)
plotResiduals(simulationOutput$scaledResiduals, tbi$distanceRiver)


## 3 Chosen model output ################################################################################

m2 <- lmer((y) ~ comparison * exposition + (PC1soil + PC2soil + PC3soil) + side + locationYear +
             (1|plot), 
           REML = T,
           data = tbi)
m2 <- lm((y) ~ comparison * exposition + (PC1soil + PC2soil + PC3soil) + side + locationYear + plot, 
           data = tbi)
isSingular(m2)
anova(m2)
### * Model output ####
MuMIn::r.squaredGLMM(m) #R2m = 0.363, R2c = 0.416
VarCorr(m2)
sjPlot::plot_model(m2, type = "re", show.values = T)
car::Anova(m2, type = 2)
anova(m2)
### * Effect sizes ####
(emm <- emmeans(m4, revpairwise ~ comparison * exposition, type = "response"))
plot(emm, comparison = T)
(emm <- emmeans(m4, revpairwise ~ side, type = "response"))

### * Save ####
table <- broom::tidy(car::Anova(m2, type = 3))
write.csv(table, here("outputs/statistics/table_anova_tbi_bc_presence.csv"))

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# C Plotten ################################################################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

library(ggbeeswarm)
library(ggeffects)

themeMB <- function(){
  theme(
    panel.background = element_rect(fill = "white"),
    text  = element_text(size = 9, color = "black"),
    axis.text  = element_text(size = 9, color = "black"),
    legend.text  = element_text(size = 8, color = "black"),
    strip.text = element_text(size = 9, color = "black"),
    axis.text.y = element_text(angle = 0, hjust = 0.5),
    axis.line.y = element_line(),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.key = element_rect(fill = "transparent"),
    legend.justification = c("right", "top"),
    legend.direction = "vertical",
    legend.position = c(.95, .45),
    legend.background = element_rect(fill = "transparent"),
    legend.margin = margin(0, 0, 0, 0, "cm"),
    plot.margin = margin(0, 0, 0, 0, "cm")
  )
}

pdata <- ggemmeans(m6, terms = c("comparison", "index","exposition"), type = "fe") %>%
  mutate(group = fct_recode(group, Losses = "B", Gains = "C", Total = "y")) %>%
  mutate(x = fct_recode(x, "'17 vs. '18" = "1718", "'18 vs. '19" = "1819", "'19 vs. '21" = "1921")) %>%
  mutate(facet = fct_recode(facet, South = "south", North = "north"))
pd <- position_dodge(.1)
(graph <- ggplot(pdata, 
                 aes(x, predicted, ymin = conf.low, ymax = conf.high, shape = group)) +
    #geom_quasirandom(tbi = tbi, aes(x, predicted), 
    #                 color = "grey70", dodge.width = .6, size = 0.7) +
    geom_errorbar(position = pd, width = 0.0, size = 0.4) +
    geom_point(position = pd, size = 2.5, fill = "white") +
    facet_grid(~ facet) +
    annotate("text", label = "n.s.", x = 2.2, y = 0.6) +
    scale_y_continuous(limits = c(0, 0.6), breaks = seq(-100, 100, 0.1)) +
    scale_shape_manual(values = c(25, 24, 16)) +
    labs(x = "", y = "Temporal beta-diversity", shape = "") +
    themeMB() +
    theme()
)

### Save ###
ggsave(here("outputs/figures/figure_4_tbi_presence_year_(800dpi_9x5cm).tiff"),
       dpi = 800, width = 9, height = 5, units = "cm")

