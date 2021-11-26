# Calculate TBI ####
# Markus Bauer
# Citation: Bauer M, Huber J, Kollmann J (202x)



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Preparation ######################################################################################
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
  select(-matches("PC.constructionYear"), -conf.low, -conf.high, -B, -C) %>%
  filter(comparison %in% c("1718", "1819", "1921") & presabu == "presence") %>%
  mutate(comparison = factor(comparison),
         locationYear = factor(locationYear),
         block = factor(block),
         plot = factor(plot),
         exposition = factor(exposition)
         ) %>%
  rename(y = D)



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Statistics #############################################################################################
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
ggplot(tbi, aes(PC2soil)) + geom_density()
ggplot(tbi, aes(exp(PC2soil))) + geom_density()

#### * check collinearity ####
data <- tbi %>%
  select(where(is.numeric), -constructionYear, -y)
#GGally::ggpairs(data, lower = list(continuous = "smooth_loess"))
#--> riverkm ~ longitude/latitude has r > 0.7 (Dormann et al. 2013 Ecography)


## 2 Model building ################################################################################

### a models ----------------------------------------------------------------------------------------

#random structure
m1a <- lmer(log(y) ~ 1 + (1|locationYear), data = tbi, REML = T)
VarCorr(m1a) 
m1b <- lmer(log(y) ~ 1 + (1|locationYear/plot), data = tbi, REML = T)
VarCorr(m1b)
m1c <- lmer(log(y) ~ 1 + (1|plot), data = tbi, REML = T)
VarCorr(m1c)
#fixed effects
m1 <- lmer(log(y) ~ (comparison + exposition + PC1soil + (PC2soil))^2 + PC3soil + side + log(distanceRiver) + locationYear +
             (1|plot), 
           REML = F,
         data = tbi)
simulateResiduals(m1, plot = T)
m2 <- lmer(log(y) ~ (comparison + exposition + PC1soil)^2 + (PC2soil) + PC3soil + side + log(distanceRiver) + locationYear +
             (1|plot), 
           REML = F, 
         data = tbi)
simulateResiduals(m2, plot = T)
m3 <- lmer(log(y) ~ (comparison + exposition + (PC2soil))^2 + PC1soil + PC3soil + side + log(distanceRiver) + locationYear +
             (1|plot), 
           REML = F, 
         data = tbi)
simulateResiduals(m3, plot = T)
m4 <- lmer(log(y) ~ comparison + exposition * PC1soil + (PC2soil) + PC3soil + side + log(distanceRiver) + locationYear +
             (1|plot), 
           REML = F, 
         data = tbi)
simulateResiduals(m4, plot = T)
m5 <- lmer(log(y) ~ comparison + exposition * (PC2soil) + PC1soil + PC3soil + side + log(distanceRiver) + locationYear + 
             (1|plot), 
           REML = F,
         data = tbi)
simulateResiduals(m5, plot = T)
m6 <- lmer(log(y) ~ comparison * exposition + PC1soil + (PC2soil) + PC3soil + side + log(distanceRiver) + locationYear +
             (1|plot), 
           REML = F,
         data = tbi)
simulateResiduals(m6, plot = T)
m7 <- lmer(log(y) ~ comparison * PC1soil + exposition + (PC2soil) + PC3soil + side + log(distanceRiver) + locationYear +
             (1|plot), 
           REML = F,
         data = tbi)
simulateResiduals(m7, plot = T)
m8 <- lmer(log(y) ~ comparison * PC2soil + exposition + PC1soil + PC3soil + side + log(distanceRiver) + locationYear +
             (1|plot), 
           REML = F,
           data = tbi)
simulateResiduals(m8, plot = T)
m9 <- lmer(log(y) ~ comparison + exposition + PC1soil + (PC2soil) + PC3soil + side + log(distanceRiver) + locationYear +
             (1|plot), 
           REML = F,
           data = tbi)
simulateResiduals(m9, plot = T)

### b comparison -----------------------------------------------------------------------------------------

AICcmodavg::aictab(cand.set = list("m1a" = m1a, "m1b" = m1b, "m1c" = m1c)) # m1c is parsimonous
AICcmodavg::aictab(cand.set = list("m1" = m1, "m2" = m2, "m3" = m3, "m4" = m4, "m5" = m5, "m6" = m6, "m7" = m7, "m8" = m8, "m9" = m9))
dotwhisker::dwplot(list(m5, m4, m9), 
                   show_intercept = F,
                   vline = geom_vline(
                     xintercept = 0,
                     colour = "grey60",
                     linetype = 2)
                   ) +
  theme_classic()
m5 <- update(m5, REML = T)
rm(list = setdiff(ls(), c("tbi", "m5")))

### c model check -----------------------------------------------------------------------------------------

simulationOutput <- simulateResiduals(m5, plot = T)
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
recalculateOutput = recalculateResiduals(simulationOutput , group = tbi$id)
car::vif(m5) # all < 3 (Zuur et al. 2010 Methods Ecol Evol)


## 3 Chosen model output ################################################################################

### * Model output ####
MuMIn::r.squaredGLMM(m5) #R2m = 0.373, R2c = 0.427
VarCorr(m5)
sjPlot::plot_model(m5, type = "re", show.values = T)
car::Anova(m5, type = 3)

### * Effect sizes ####
(emm <- emmeans(m5, revpairwise ~ side, type = "response"))
(emm <- emmeans(m5, revpairwise ~ comparison, type = "response"))
plot(emm, comparison = T)

### * Save ####
table <- broom::tidy(car::Anova(m5, type = 3))
write.csv(table, here("outputs/statistics/table_anova_tbi_d_presence.csv"))
