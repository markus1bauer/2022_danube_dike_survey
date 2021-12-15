# Calculate TBI ####
# Markus Bauer
# Citation: Bauer M, Huber J, Kollmann J (2022)



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Preparation ######################################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


### Packages ###
library(here)
library(tidyverse)
library(ggbeeswarm)
library(lmerTest)
library(DHARMa)
library(emmeans)

### Start ###
rm(list = ls())
setwd(here("data/processed"))

### Load data ###
tbi <- read_csv("data_processed_tbi.csv", col_names = T, na = c("", "na", "NA"), col_types = 
                  cols(
                    .default = "?"
                  )) %>%
  filter(comparison %in% c("1718", "1819", "1921")) %>%
  pivot_wider(id_cols = c(plot, comparison, exposition, side, block, 
                          location, locationYear, constructionYear, 
                          longitude, latitude, riverkm, distanceRiver, 
                          PC1soil, PC2soil, PC3soil, conf.low, conf.high), 
              names_from = "presabu", 
              values_from = "D") %>%
  mutate(plot = factor(plot),
         block = factor(block),
         comparison = factor(comparison),
         exposition = factor(exposition),
         side = factor(side),
         constructionYear = factor(constructionYear),
         locationYear = factor(locationYear)) %>%
  rename(y = presence)

data_collinearity <- tbi %>%
  select(where(is.numeric))

tbi <- tbi %>%
  mutate(across(c("longitude", "latitude", "riverkm", "distanceRiver"), scale))



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Statistics #############################################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


## 1 Data exploration #####################################################################################

#### * Graphs #####
#main
ggplot(tbi, aes(x = abundance, y = y)) + 
  geom_point() +
  geom_smooth(method = "loess")
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
ggplot(tbi, aes(x = (distanceRiver), y = log(y))) + 
  geom_point() +
  geom_smooth(method = "lm")
ggplot(tbi, aes(x = constructionYear, y = y)) + 
  geom_point() + 
  geom_smooth(method = "loess")
ggplot(tbi, aes(x = PC1soil, y = (y))) + 
  geom_point() +
  geom_smooth(method = "lm")
ggplot(tbi, aes(x = exp(PC2soil), y = y)) + 
  geom_point() +
  geom_smooth(method = "lm")
#2way
ggplot(tbi, aes(x = exposition, y = y, color = comparison)) + 
  geom_boxplot() +
  geom_quasirandom(dodge.width = .8)
ggplot(tbi, aes(x = PC1soil, y = y, color = comparison)) + 
  geom_point() +
  geom_smooth(method = "lm")
ggplot(tbi, aes(x = PC2soil, y = y, color = comparison)) + 
  geom_point() +
  geom_smooth(method = "lm")
ggplot(tbi, aes(x = PC1soil, y = y, color = exposition)) + 
  geom_point() +
  geom_smooth(method = "lm")
ggplot(tbi, aes(x = (PC2soil), y = y, color = exposition)) + 
  geom_point() +
  geom_smooth(method = "lm")

#### * Outliers, zero-inflation, transformations? ####
dotchart(log(tbi$y), groups = factor(tbi$exposition), main = "Cleveland dotplot")
tbi %>% count(locationYear)
tbi %>% count(plot) %>% count(n)
boxplot(tbi$y);#identify(rep(1, length(etbi$rgr13)), etbi$rgr13, labels = c(etbi$n))
plot(table((tbi$y)), type = "h", xlab = "Observed values", ylab = "Frequency")
ggplot(tbi, aes(y)) + geom_density()
ggplot(tbi, aes(log(y))) + geom_density()

#### * check collinearity ####
#GGally::ggpairs(data_collinearity, lower = list(continuous = "smooth_loess"))
#--> riverkm ~ longitude/latitude has r > 0.7 (Dormann et al. 2013 Ecography)
rm(data_collinearity)


## 2 Model building ################################################################################

### a models ----------------------------------------------------------------------------------------

### * random structure ####
m1a <- blmer(log(y) ~ 1 + (1|locationYear), data = tbi, REML = T)
m1b <- blmer(log(y) ~ 1 + (1|locationYear/plot), data = tbi, REML = T)
m1c <- blmer(log(y) ~ 1 + (1|plot), data = tbi, REML = T)
MuMIn::AICc(m1a, m1b, m1c) # m1c most parsimonous

#### * fixed effects ####
m1 <- blmer(log(y) ~ (comparison + exposition + PC1soil)^2 + PC2soil + PC3soil + side + distanceRiver + locationYear + abundance + 
             (1|plot), 
           REML = F,
           control = lmerControl(optimizer = "Nelder_Mead"),
           cov.prior = wishart,
           data = tbi)
simulateResiduals(m1, plot = T)
m2 <- blmer(log(y) ~ comparison + exposition * PC1soil + PC2soil + PC3soil + side + distanceRiver + locationYear + abundance + 
             (1|plot), 
           REML = F,
           control = lmerControl(optimizer = "Nelder_Mead"),
           cov.prior = wishart,
           data = tbi)
simulateResiduals(m2, plot = T)
m3 <- blmer(log(y) ~ comparison * exposition + PC1soil + PC2soil + PC3soil + side + distanceRiver + locationYear + abundance + 
             (1|plot), 
            REML = F,
            control = lmerControl(optimizer = "Nelder_Mead"),
            cov.prior = wishart,
            data = tbi)
simulateResiduals(m3, plot = T)
m4 <- blmer(log(y) ~ comparison * PC1soil + exposition + PC2soil + PC3soil + side + distanceRiver + locationYear + abundance + 
             (1|plot), 
           REML = F,
           control = lmerControl(optimizer = "Nelder_Mead"),
           cov.prior = wishart,
           data = tbi)
simulateResiduals(m4, plot = T)
m5 <- blmer(log(y) ~ comparison + exposition + PC1soil + PC2soil + PC3soil + side + distanceRiver + locationYear + abundance + 
             (1|plot), 
           REML = F,
           control = lmerControl(optimizer = "Nelder_Mead"),
           cov.prior = wishart,
           data = tbi)
simulateResiduals(m5, plot = T)

### b comparison -----------------------------------------------------------------------------------------

MuMIn::AICc(m1, m2, m3, m4, m5) # m2 most parsimonious; Use AICc and not AIC since ratio n/K < 40 (Burnahm & Anderson 2002 p. 66)
dotwhisker::dwplot(list(m2, m5), 
                   show_intercept = F,
                   vline = geom_vline(
                     xintercept = 0,
                     colour = "grey60",
                     linetype = 2)
                   ) +
  theme_classic()
m <- update(m2, REML = T)
rm(list = setdiff(ls(), c("tbi", "m")))

### c model check -----------------------------------------------------------------------------------------

simulationOutput <- simulateResiduals(m, plot = T)
plotResiduals(simulationOutput$scaledResiduals, tbi$locationYear)
plotResiduals(simulationOutput$scaledResiduals, tbi$plot)
plotResiduals(simulationOutput$scaledResiduals, tbi$comparison)
plotResiduals(simulationOutput$scaledResiduals, tbi$exposition)
plotResiduals(simulationOutput$scaledResiduals, tbi$side)
plotResiduals(simulationOutput$scaledResiduals, tbi$PC1soil)
plotResiduals(simulationOutput$scaledResiduals, tbi$PC2soil)
plotResiduals(simulationOutput$scaledResiduals, tbi$PC3soil)
plotResiduals(simulationOutput$scaledResiduals, tbi$distanceRiver)
plotResiduals(simulationOutput$scaledResiduals, tbi$abundance)
plotResiduals(simulationOutput$scaledResiduals, tbi$riverkm)
car::vif(m) # remove riverkm since > 3 oder 10 (Zuur et al. 2010 Methods Ecol Evol)


## 3 Chosen model output ################################################################################

### * Model output ####
MuMIn::r.squaredGLMM(m) #R2m = 0.375, R2c = 0.534
VarCorr(m)
sjPlot::plot_model(m, type = "re", show.values = T)
dotwhisker::dwplot(m, 
                   show_intercept = F,
                   vline = geom_vline(
                     xintercept = 0,
                     colour = "grey60",
                     linetype = 2)
                   ) +
  theme_classic()

### * Effect sizes ####
(emm <- emmeans(m, revpairwise ~ side, type = "response"))
plot(emm, comparison = T)
(emm <- emmeans(m, revpairwise ~ comparison, type = "response"))
sjPlot::plot_model(m, type = "emm", terms = c("PC1soil", "exposition"), show.data = T)
