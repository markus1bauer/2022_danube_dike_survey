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
                    plot = "c",
                    locationYear = "f",
                    exposition = "f",
                    side = "c",
                    comparison = "f"
                  )) %>%
  select(-matches("PC.constructionYear"), -conf.low, -conf.high) %>%
  filter(comparison %in% c("1718", "1819", "1921") & presabu == "abundance") %>%
  mutate(side = if_else(side == "water_creek", "water", side),
         ageGroup = if_else(constructionYear %in% c(2002, 2003), "0203", if_else(
           constructionYear %in% c(2006, 2007), "0607", if_else(
             constructionYear == 2008, "2008", if_else(
               constructionYear %in% c(2010, 2011), "1011", if_else(
                 constructionYear %in% c(2012, 2013), "1213", "other"
               ))))),
         constructionYearF = if_else(constructionYear %in% c(2003, 2006), 2003.6, constructionYear),
         constructionYearF = factor(constructionYearF),
         ageGroup = fct_relevel(ageGroup, "2008", after = 2),
         plot = factor(plot),
         comparison = factor(comparison),
         exposition = factor(exposition),
         side = factor(side),
         ageGroup = factor(ageGroup),
         locationYear = factor(locationYear)) %>%
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
  geom_smooth(method = "lm")
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
tbi %>% count(constructionYearF)
boxplot(tbi$y);#identify(rep(1, length(etbi$rgr13)), etbi$rgr13, labels = c(etbi$n))
plot(table((tbi$y)), type = "h", xlab = "Observed values", ylab = "Frequency")
ggplot(tbi, aes(y)) + geom_density()
ggplot(tbi, aes(log(y))) + geom_density()
ggplot(tbi, aes(PC2soil)) + geom_density()
ggplot(tbi, aes(exp(PC2soil))) + geom_density()

#### * check collinearity ####
data <- tbi %>%
  select(where(is.numeric), -constructionYear, -B, -C, -y)
GGally::ggpairs(data, lower = list(continuous = "smooth_loess"))
#--> riverkm ~ longitutde/latitude has r > 0.7 (Dormann et al. 2013 Ecography) --> longitude and latitude have to be excluded


## 2 Model building ################################################################################

### a models ----------------------------------------------------------------------------------------

#random structure
m1a <- lmer(y ~ 1 + (1|locationAbb), data = tbi, REML = T)
VarCorr(m1a) 
m1b <- lmer(y ~ 1 + (1|locationAbb/plot), data = tbi, REML = T)
VarCorr(m1b)
m1c <- lmer(y ~ 1 + (1|plot), data = tbi, REML = T)
VarCorr(m1c)
#fixed effects
m1 <- lmer(y ~ scale(longitude) + scale(latitude) + (comparison + exposition + PC1soil + PC2soil)^2 + PC3soil + side + constructionYearF + log(distanceRiver) +
             (1|plot), 
           REML = F,
           data = tbi)
simulateResiduals(m1, plot = T)
m2 <- lmer(y ~ scale(longitude) + scale(latitude) + (comparison + exposition + PC1soil)^2 + (PC2soil) + PC3soil + side + constructionYearF + log(distanceRiver) +
             (1|plot), 
           REML = F, 
           data = tbi)
simulateResiduals(m2, plot = T)
m3 <- lmer(y ~ scale(longitude) + scale(latitude) + (comparison + exposition + (PC2soil))^2 + PC1soil + PC3soil + side + constructionYearF + log(distanceRiver) +
             (1|plot), 
           REML = F, 
           data = tbi)
simulateResiduals(m3, plot = T)
m4 <- lmer(y ~ scale(longitude) + scale(latitude) + comparison + exposition * PC1soil + (PC2soil) + PC3soil + side + constructionYearF + log(distanceRiver) +
             (1|plot), 
           REML = F, 
           data = tbi)
simulateResiduals(m4, plot = T)
m5 <- lmer(y ~ scale(longitude) + scale(latitude) + comparison + exposition * (PC2soil) + PC1soil + PC3soil + side + constructionYearF + log(distanceRiver) +
             (1|plot), 
           REML = F,
           data = tbi)
simulateResiduals(m5, plot = T)
m6 <- lmer(y ~ scale(longitude) + scale(latitude) + comparison * exposition + PC1soil + (PC2soil) + PC3soil + side + constructionYearF + log(distanceRiver) +
             (1|plot), 
           REML = F,
           data = tbi)
simulateResiduals(m6, plot = T)
m7 <- lmer(y ~ scale(longitude) + scale(latitude) + comparison * PC1soil + exposition + (PC2soil) + PC3soil + side + constructionYearF + log(distanceRiver) +
             (1|plot), 
           REML = F,
           data = tbi)
simulateResiduals(m7, plot = T)
m8 <- lmer(y ~ comparison + exposition + PC1soil + (PC2soil) + PC3soil + side + log(distanceRiver) + locationYear + 
             (1|plot), 
           REML = F,
           data = tbi)
isSingular(m8);simulateResiduals(m8, plot = T)

### b comparison -----------------------------------------------------------------------------------------

AICcmodavg::aictab(cand.set = list("m1a" = m1a, "m1b" = m1b, "m1c" = m1c))
AICcmodavg::aictab(cand.set = list("m1" = m1, "m2" = m2, "m3" = m3, "m4" = m4, "m5" = m5, "m6" = m6, "m7" = m7, "m8" = m8))
dotwhisker::dwplot(list(m8, m4, m7), show_intercept = T)
m8 <- update(m8, REML = T)
rm(list = setdiff(ls(), c("tbi", "m8")))

#### c model check -----------------------------------------------------------------------------------------
simulationOutput <- simulateResiduals(m8, plot = T)
plotResiduals(simulationOutput$scaledResiduals, tbi$locationYear)
plotResiduals(simulationOutput$scaledResiduals, tbi$block)
plotResiduals(simulationOutput$scaledResiduals, tbi$plot)
plotResiduals(simulationOutput$scaledResiduals, tbi$comparison)
plotResiduals(simulationOutput$scaledResiduals, tbi$exposition)
plotResiduals(simulationOutput$scaledResiduals, tbi$side)
plotResiduals(simulationOutput$scaledResiduals, tbi$PC1soil)
plotResiduals(simulationOutput$scaledResiduals, tbi$PC2soil)
plotResiduals(simulationOutput$scaledResiduals, tbi$PC3soil)
plotResiduals(simulationOutput$scaledResiduals, tbi$distanceRiver)
car::vif(m8) # all < 3 (Zuur et al. 2010 Methods Ecol Evol)


### 3 Chosen model output ################################################################################

### * Model output ####
MuMIn::r.squaredGLMM(m8) #R2m = 0.271, R2c = 0.367
VarCorr(m8)
sjPlot::plot_model(m8, type = "re", show.values = T)
car::Anova(m8, type = 2)

### * Save ####
table <- broom::tidy(car::Anova(m8, type = 2))
write.csv(table, here("outputs/statistics/table_anova_tbi_d_abundance.csv"))