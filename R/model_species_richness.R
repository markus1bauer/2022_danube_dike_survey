# Model for graminoid's cover ratio ####
# Markus Bauer
# Citation: Markus Bauer & Harald Albrecht (2020) Basic and Applied Ecology 42, 15-26
# https://doi.org/10.1016/j.baae.2019.11.003



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Preparation ################################################################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

### Packages ###
library(tidyverse)
library(ggbeeswarm)
library(lmerTest)
library(DHARMa)
library(MuMIn)
library(car)
library(emmeans)

### Start ###
rm(list = ls())
setwd("Z:/Documents/0_Uni/Projekt_9_Masterthesis/3_Aufnahmen_und_Ergebnisse/Basic_Appl_Ecol/data/processed")

### Load data ###
sites <- read_csv2("data_processed_sites849318.csv", col_names = T, na = "na", col_types = 
                    cols(
                      .default = col_double(),
                      ID = col_factor(),
                      plot = col_factor(),
                      block = col_factor(),
                      dataset = col_factor(),
                      year = col_factor()
                    )        
)

(sites <- sites %>% 
    select(ID, plot, block, year, all, target, rlg, rlb) %>%
    pivot_longer(c(all, target, rlg, rlb), names_to = "type", values_to = "n")
  )
sites$type <- factor(sites$type)


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Statistics ################################################################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


### 1 Data exploration #####################################################################################

#### a Graphs ---------------------------------------------------------------------------------------------
#simple effects:
par(mfrow = c(2,2))
plot(n ~ year, sites)
plot(n ~ type, sites)
plot(n ~ block, sites)
#2way type:year
ggplot(sites, aes(type, n, color = year)) + geom_boxplot() + geom_quasirandom(dodge.width = .7, groupOnX = T)
#2way with block:
ggplot(sites, aes(block, n, color = year)) + geom_boxplot() + geom_quasirandom(dodge.width = .7, groupOnX = T)
ggplot(sites, aes(block, n, color = type)) + geom_boxplot() + geom_quasirandom(dodge.width = .7, groupOnX = T)

##### b Outliers, zero-inflation, transformations? -----------------------------------------------------
par(mfrow = c(2,2))
dotchart((sites$n), groups = factor(sites$year), main = "Cleveland dotplot")
dotchart((sites$n), groups = factor(sites$type), main = "Cleveland dotplot")
dotchart((sites$n), groups = factor(sites$block), main = "Cleveland dotplot")
par(mfrow=c(1,1));
boxplot(sites$n);#identify(rep(1, length(edata$rgr13)), edata$rgr13, labels = c(edata$n))
plot(table((sites$n)), type = "h", xlab = "Observed values", ylab = "Frequency")
ggplot(sites, aes(n)) + geom_density()
ggplot(sites, aes(sqrt(n))) + geom_density()


## 2 Model building ################################################################################

#### a models ----------------------------------------------------------------------------------------
#random structure
m1 <- lmer(n ~ year * type + (1|block/plot), sites, REML = F)
VarCorr(m1)
#2w-model
m2 <- lmer(sqrt(n) ~ year * type + (1|block/plot), sites, REML = F)
isSingular(m2)
simulateResiduals(m2, plot = T)
#1w-model
m3 <- lmer(sqrt(n) ~ year + type + (1|block/plot), sites, REML = F)
isSingular(m3)
simulateResiduals(m3, plot = T)
#2w glmm model
m4 <- glmer(n ~ year * type + (1|block/plot), sites, family = poisson)
isSingular(m4)
simulateResiduals(m4, plot = T)
#1w glmm model
m5 <- glmer(n ~ year + type + (1|block/plot), sites, family = poisson)
isSingular(m5)
simulateResiduals(m5, plot = T)

#### b comparison -----------------------------------------------------------------------------------------
anova(m2,m3)
anova(m4,m5)
rm(m1,m3,m4,m5)

#### c model check -----------------------------------------------------------------------------------------
simulationOutput <- simulateResiduals(m2, plot = T)
par(mfrow=c(2,2));
plotResiduals(main = "year", simulationOutput$scaledResiduals, sites$year)
plotResiduals(main = "type", simulationOutput$scaledResiduals, sites$type)
plotResiduals(main = "block", simulationOutput$scaledResiduals, sites$block)


## 3 Chosen model output ################################################################################

### Model output ---------------------------------------------------------------------------------------------
m2 <- lmer(sqrt(n) ~ year * type + (1|block/plot), sites, REML = F)
MuMIn::r.squaredGLMM(m2) #R2m = 0.75, R2c = 0.92
VarCorr(m2)
sjPlot::plot_model(m2, type = "re", show.values = T)
car::Anova(m2, type = 3)

### Effect sizes -----------------------------------------------------------------------------------------
(emm <- emmeans(m2, revpairwise ~ year | type, type = "response"))
plot(emm, comparison = T)
contrast(emmeans(m2, ~ year | type, type = "response"), "trt.vs.ctrl", ref = 3)
