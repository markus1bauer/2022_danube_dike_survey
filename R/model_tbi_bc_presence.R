# Beta diversity on dike grasslands
# Ratio of gains and losses of TBI (presence-absence data) ####
# Markus Bauer
# 2022-01-11
# Citation: 
## Bauer M, Huber J, Kollmann J (submitted) 
## Balanced turnover is a main aspect of biodiversity on restored dike grasslands: not only deterministic environmental effects, but also non-directional year and site effects drive spatial and temporal beta diversity.
## Unpublished data.



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Preparation ################################################################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


### Packages ###
library(here)
library(tidyverse)
library(ggbeeswarm)
library(blme)
library(DHARMa)
library(emmeans)

### Start ###
rm(list = ls())
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
  mutate(y = C_presence - B_presence)

data_collinearity <- sites %>%
  select(ends_with("ude"), riverkm, distanceRiver, starts_with("PC"))

sites <- sites %>%
  mutate(across(c("longitude", "latitude", "riverkm", "distanceRiver"), scale))



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Statistics ################################################################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


## 1 Data exploration #####################################################################################

#### * Graphs #####
#main
ggplot(sites, aes(x = comparison, y = y)) + 
  geom_boxplot() +
  geom_quasirandom()
ggplot(sites, aes(x = exposition, y = y)) + 
  geom_boxplot() +
  geom_quasirandom()
ggplot(sites, aes(x = side, y = y)) + 
  geom_boxplot() +
  geom_quasirandom()
ggplot(sites, aes(x = riverkm, y = (y))) + 
  geom_point() +
  geom_smooth(method = "loess")
ggplot(sites, aes(x = (distanceRiver), y = (y))) + 
  geom_point() +
  geom_smooth(method = "loess")
ggplot(sites, aes(x = as.double(constructionYear), y = y)) + 
  geom_point() + 
  geom_smooth(method = "loess")
ggplot(sites, aes(x = PC1soil, y = (y))) + 
  geom_point() +
  geom_smooth(method = "lm")
ggplot(sites, aes(x = (PC2soil), y = y)) + 
  geom_point() +
  geom_smooth(method = "lm")
#2way
ggplot(sites, aes(x = exposition, y = y, color = comparison)) + 
  geom_boxplot() +
  geom_quasirandom(dodge.width = .8)
ggplot(sites, aes(x = PC1soil, y = y, color = comparison)) + 
  geom_point() +
  geom_smooth(method = "lm")
ggplot(sites, aes(x = PC2soil, y = y, color = comparison)) + 
  geom_point() +
  geom_smooth(method = "lm")
ggplot(sites, aes(x = PC1soil, y = y, color = exposition)) + 
  geom_point() +
  geom_smooth(method = "lm")
ggplot(sites, aes(x = (PC2soil), y = y, color = exposition)) + 
  geom_point() +
  geom_smooth(method = "lm")

#### * Outliers, zero-inflation, transformations? ####

dotchart((sites$y), groups = factor(sites$exposition), main = "Cleveland dotplot")
sites %>% count(locationYear)
sites %>% count(plot) %>% count(n)
boxplot(sites$y);#identify(rep(1, length(etbi$rgr13)), etbi$rgr13, labels = c(etbi$n))
plot(table((sites$y)), type = "h", xlab = "Observed values", ylab = "Frequency")
ggplot(sites, aes(y)) + geom_density()
ggplot(sites, aes(log(y))) + geom_density()

#### * check collinearity ####
#GGally::ggpairs(data_collinearity, lower = list(continuous = "smooth_loess"))
#--> riverkm ~ longitude/latitude has r > 0.7 (Dormann et al. 2013 Ecography)
rm(data_collinearity)


## 2 Model building ################################################################################

### a models ----------------------------------------------------------------------------------------

### * Random structure ####
m1a <- blmer(y ~ 1 + (1|locationYear), data = sites, REML = T)
m1b <- blmer(y ~ 1 + (1|locationYear/plot), data = sites, REML = T)
m1c <- blmer(y ~ 1 + (1|plot), data = sites, REML = T)
MuMIn::AICc(m1a, m1b, m1c) #m1b most parsimonious

### * Fixed effects ####
m1 <- blmer(y ~ (comparison + exposition + PC1soil)^2 + PC2soil + PC3soil + side + distanceRiver + locationYear + 
              (1|plot), 
            REML = F,
            control = lmerControl(optimizer = "Nelder_Mead"),
            cov.prior = wishart,
            data = sites)
simulateResiduals(m1, plot = T)
m2 <- blmer(y ~ comparison + exposition * PC1soil + PC2soil + PC3soil + side + distanceRiver + locationYear + 
              (1|plot), 
            REML = F,
            control = lmerControl(optimizer = "Nelder_Mead"),
            cov.prior = wishart,
            data = sites)
simulateResiduals(m2, plot = T)
m3 <- blmer(y ~ comparison * exposition + PC1soil + PC2soil + PC3soil + side + distanceRiver + locationYear + 
              (1|plot), 
            REML = F,
            control = lmerControl(optimizer = "Nelder_Mead"),
            cov.prior = wishart,
            data = sites)
simulateResiduals(m3, plot = T)
m4 <- blmer(y ~ comparison * PC1soil + exposition + PC2soil + PC3soil + side + distanceRiver + locationYear + 
              (1|plot), 
            REML = F,
            control = lmerControl(optimizer = "Nelder_Mead"),
            cov.prior = wishart,
            data = sites)
simulateResiduals(m4, plot = T)
m5 <- blmer(y ~ comparison + exposition + PC1soil + PC2soil + PC3soil + side + distanceRiver + locationYear + 
              (1|plot), 
            REML = F,
            control = lmerControl(optimizer = "Nelder_Mead"),
            cov.prior = wishart,
            data = sites)
simulateResiduals(m5, plot = T)

### b comparison -----------------------------------------------------------------------------------------

MuMIn::AICc(m1, m2, m3, m4, m5) # m4 most parsimonious; Use AICc and not AIC since ratio n/K < 40 (Burnahm & Anderson 2002 p. 66)
dotwhisker::dwplot(list(m4, m3), #m4 bad model critique m3 ok 
                   show_intercept = F,
                   vline = geom_vline(
                     xintercept = 0,
                     colour = "grey60",
                     linetype = 2)
                   ) +
  theme_classic()
m <- update(m3, REML = T)
rm(list = setdiff(ls(), c("sites", "m")))

### c model check -----------------------------------------------------------------------------------------

simulationOutput <- simulateResiduals(m, plot = T)
plotResiduals(simulationOutput$scaledResiduals, sites$locationYear)
plotResiduals(simulationOutput$scaledResiduals, sites$plot)
plotResiduals(simulationOutput$scaledResiduals, sites$comparison)
plotResiduals(simulationOutput$scaledResiduals, sites$exposition)
plotResiduals(simulationOutput$scaledResiduals, sites$side)
plotResiduals(simulationOutput$scaledResiduals, sites$PC1soil)
plotResiduals(simulationOutput$scaledResiduals, sites$PC2soil)
plotResiduals(simulationOutput$scaledResiduals, sites$PC3soil)
plotResiduals(simulationOutput$scaledResiduals, sites$distanceRiver)
plotResiduals(simulationOutput$scaledResiduals, sites$riverkm)
car::vif(m) # --> remove riverkm > 3 oder 10 (Zuur et al. 2010 Methods Ecol Evol) 


## 3 Chosen model output ################################################################################

### * Model output ####
MuMIn::r.squaredGLMM(m) #R2m = 0.413, R2c = 0.438
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
(emm <- emmeans(m, revpairwise ~ comparison | exposition, type = "response"))
plot(emm, comparison = T)
sjPlot::plot_model(m, type = "emm", terms = c("comparison", "exposition"), show.data = F)
