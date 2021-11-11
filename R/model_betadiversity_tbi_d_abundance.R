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
library(adespatial)
library(lme4)  
library(DHARMa)
library(AICcmodavg)
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
                    plot = "c",
                    locationYear = "f",
                    exposition = "f",
                    side = "c",
                    comparison = "f"
                  )) %>%
  select(-matches("PC.constructionYear"), -conf.low, -conf.high, -locationYear) %>%
  filter(comparison %in% c("1718", "1819", "1921") & presabu == "abundance") %>%
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
  rename(y = D)



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Statistics ################################################################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


## 1 Data exploration #####################################################################################

#### * Graphs #####

#main
ggplot(tbi, aes(x = comparison, y = y)) + 
  geom_boxplot() +
  geom_quasirandom()
ggplot(tbi, aes(x = exposition, y = y)) + 
  geom_boxplot() +
  geom_quasirandom()
ggplot(tbi, aes(x = side, y = y)) + 
  geom_boxplot() +
  geom_quasirandom()
ggplot(tbi, aes(x = riverkm, y = (y))) + 
  geom_point() +
  geom_smooth(method = "loess")
ggplot(tbi, aes(x = log(distanceRiver), y = (y))) + 
  geom_point() +
  geom_smooth(method = "loess")
ggplot(tbi, aes(x = locationYear, y = y)) + 
  geom_boxplot()
ggplot(tbi, aes(x = ageGroup, y = y)) + 
  geom_boxplot()
ggplot(tbi, aes(x = constructionYearF, y = y)) + 
  geom_boxplot()
ggplot(tbi, aes(x = PC1soil, y = (y))) + 
  geom_point() +
  geom_smooth(method = "loess")
ggplot(tbi, aes(x = exp(PC2soil), y = y)) + 
  geom_point() +
  geom_smooth(method = "loess")
#2way
ggplot(tbi, aes(x = exposition, y = y, color = comparison)) + 
  geom_boxplot() +
  geom_quasirandom(dodge.width = .8)
ggplot(tbi, aes(x = PC1soil, y = y, color = comparison)) + 
  geom_point() +
  geom_smooth(method = "loess")
ggplot(tbi, aes(x = PC2soil, y = y, color = comparison)) + 
  geom_point() +
  geom_smooth(method = "loess")
ggplot(tbi, aes(x = PC1soil, y = y, color = exposition)) + 
  geom_point() +
  geom_smooth(method = "loess")
ggplot(tbi, aes(x = exp(PC2soil), y = y, color = exposition)) + 
  geom_point() +
  geom_smooth(method = "lm")
ggplot(tbi, aes(x = PC1soil, y = y, color = ageGroup)) + 
  geom_point() +
  geom_smooth(method = "lm")
#3way
ggplot(tbi, aes(x = constructionYear, y = y, color = comparison)) + 
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~exposition)
ggplot(tbi, aes(x = PC1soil, y = y, color = comparison)) + 
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~exposition)
ggplot(tbi, aes(x = PC2soil, y = y, color = comparison)) + 
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~exposition)

#### * Outliers, zero-inflation, transformations? ####

dotchart((tbi$y), groups = factor(tbi$exposition), main = "Cleveland dotplot")
tbi %>% count(locationYear)
boxplot(tbi$y);#identify(rep(1, length(etbi$rgr13)), etbi$rgr13, labels = c(etbi$n))
plot(table((tbi$y)), type = "h", xlab = "Observed values", ylab = "Frequency")
ggplot(tbi, aes(y)) + geom_density()
ggplot(tbi, aes(log(y))) + geom_density()
ggplot(tbi, aes(PC2soil)) + geom_density()
ggplot(tbi, aes(exp(PC2soil))) + geom_density()

#### * check collinearity ####
data <- tbi %>%
  select(where(is.numeric), -conf.low, -conf.high, -constructionYear, -B, -C, -y)
GGally::ggpairs(data, lower = list(continuous = "smooth_loess"))
#--> xx ~ xx has r > 0.7 (Dormann et al. 2013 Ecography) --> xx has to be excluded


## 2 Model building ################################################################################

### a models ----------------------------------------------------------------------------------------

#random structure
m1a <- lmer(y ~ 1 + (1|locationYear), data = tbi, REML = T)
VarCorr(m1a) 
m1b <- lmer(y ~ 1 + (1|locationYear/plot), data = tbi, REML = T)
VarCorr(m1b)
m1c <- lmer(y ~ 1 + (1|plot), data = tbi, REML = T)
VarCorr(m1c)
#fixed effects
m1 <- lmer(log(y) ~ locationAbb + (comparison + exposition + PC1soil + PC2soil)^2 + PC3soil + side + constructionYearF + log(distanceRiver) +
             (1|plot), 
           REML = F,
           data = tbi)
simulateResiduals(m1, plot = T)
m2 <- lmer(log(y) ~ locationAbb + (comparison + exposition + PC1soil)^2 + (PC2soil) + PC3soil + side + constructionYearF + log(distanceRiver) +
             (1|plot), 
           REML = F, 
           data = tbi)
simulateResiduals(m2, plot = T)
m3 <- lmer(log(y) ~ locationAbb + (comparison + exposition + (PC2soil))^2 + PC1soil + PC3soil + side + constructionYearF + log(distanceRiver) +
             (1|plot), 
           REML = F, 
           data = tbi)
simulateResiduals(m3, plot = T)
m4 <- lmer(log(y) ~ locationAbb + comparison + exposition * PC1soil + (PC2soil) + PC3soil + side + constructionYearF + log(distanceRiver) +
             (1|plot), 
           REML = F, 
           data = tbi)
simulateResiduals(m4, plot = T)
m5 <- lmer(log(y) ~ comparison + exposition * (PC2soil) + PC1soil + PC3soil + side + constructionYearF + log(distanceRiver) +
             (1|plot), 
           REML = F,
           data = tbi)
simulateResiduals(m5, plot = T)
m6 <- lmer(log(y) ~ locationAbb + comparison * exposition + PC1soil + (PC2soil) + PC3soil + side + constructionYearF + log(distanceRiver) +
             (1|plot), 
           REML = F,
           data = tbi)
simulateResiduals(m6, plot = T)
m7 <- lmer(log(y) ~ locationAbb + comparison * PC1soil + exposition + (PC2soil) + PC3soil + side + constructionYearF + log(distanceRiver) +
             (1|plot), 
           REML = F,
           data = tbi)
simulateResiduals(m7, plot = T)
m8 <- lmer(log(y) ~ locationAbb + comparison + exposition + PC1soil + (PC2soil) + PC3soil + side + constructionYearF + log(distanceRiver) +
             (1|plot), 
           REML = F,
           data = tbi)
simulateResiduals(m8, plot = T)

### b comparison -----------------------------------------------------------------------------------------

aictab(cand.set = list("m1a" = m1a, "m1b" = m1b, "m1c" = m1c))
aictab(cand.set = list("m1" = m1, "m2" = m2, "m3" = m3, "m4" = m4, "m5" = m5, "m6" = m6, "m7" = m7, "m8" = m8))
sjPlot::plot_model(m5, type = "emm", terms = c("PC2soil", "exposition"), show.data = T)
ggsave(here("outputs/figures/figure_tbi_d_presence_PC2_exposition_(800dpi_9x5cm).tiff"),
       dpi = 800, width = 9, height = 5, units = "cm")
sjPlot::plot_model(m5, type = "emm", terms = c("side"), show.data = F)
sjPlot::plot_model(m5, type = "emm", terms = c("comparison"), show.data = F)
sjPlot::plot_model(m5, type = "emm", terms = c("exposition"), show.data = F)
sjPlot::plot_model(m5, type = "emm", terms = c("constructionYearF"), show.data = F)
sjPlot::plot_model(m5, type = "emm", terms = c("distanceRiver"), show.data = F)
sjPlot::plot_model(m5)
ggsave(here("outputs/figures/figure_tbi_d_presence_(800dpi_9x5cm).tiff"),
       dpi = 800, width = 9, height = 5, units = "cm")
sjPlot::plot_model(m4, type = "emm", terms = c("PC1soil", "exposition"), show.data = T)
dotwhisker::dwplot(list(m5, m4, m8), show_intercept = T)
rm(m1a, m1b, m1c, m2, m3, m4, m5)
car::Anova(m5, type = 3)

#### c model check -----------------------------------------------------------------------------------------

simulationOutput <- simulateResiduals(m4, plot = T)
plotResiduals(simulationOutput$scaledResiduals, tbi$locationYear)
plotResiduals(simulationOutput$scaledResiduals, tbi$block)
plotResiduals(simulationOutput$scaledResiduals, tbi$plot)
plotResiduals(simulationOutput$scaledResiduals, tbi$plotAge)
plotResiduals(simulationOutput$scaledResiduals, tbi$comparison)
plotResiduals(simulationOutput$scaledResiduals, tbi$exposition)
plotResiduals(simulationOutput$scaledResiduals, tbi$side)
plotResiduals(simulationOutput$scaledResiduals, tbi$PC1soil)
plotResiduals(simulationOutput$scaledResiduals, tbi$PC2soil)
plotResiduals(simulationOutput$scaledResiduals, tbi$PC3soil)
plotResiduals(simulationOutput$scaledResiduals, tbi$distanceRiver)


### 3 Chosen model output ################################################################################

m4 <- lmer((y) ~ comparison + exposition + PC1soil + PC2soil + PC3soil + side + plotAge + locationYear + 
             (1|plot), 
           REML = T,
           data = tbi)


### * Model output ####
MuMIn::r.squaredGLMM(m4) #R2m = 0.252, R2c = 0.376
VarCorr(m4)
sjPlot::plot_model(m4, type = "re", show.values = T)
car::Anova(m4, type = 2)

### * Effect sizes ####
(emm <- emmeans(m4, revpairwise ~ comparison, type = "response"))
plot(emm, comparison = T)
(emm <- emmeans(m4, revpairwise ~ side, type = "response"))

### * Save ####
table <- broom::tidy(car::Anova(m4, type = 2))
write.csv(table, here("outputs/statistics/table_anova_tbi_d_presence.csv"))

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

