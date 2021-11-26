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
  filter(comparison %in% c("1718", "1819", "1921") & presabu == "presence") %>%
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
  mutate(y = C - B)



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
ggplot(tbi, aes(x = locationAbb, y = y)) + 
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
tbi %>% count(locationAbb)
tbi %>% count(constructionYearF)
boxplot(tbi$y);#identify(rep(1, length(etbi$rgr13)), etbi$rgr13, labels = c(etbi$n))
plot(table((tbi$y)), type = "h", xlab = "Observed values", ylab = "Frequency")
ggplot(tbi, aes(y)) + geom_density()
ggplot(tbi, aes(log(y))) + geom_density()
ggplot(tbi, aes(PC2soil)) + geom_density()
ggplot(tbi, aes(exp(PC2soil))) + geom_density()

#### * check collinearity ####
data <- tbi %>%
  select(where(is.numeric), -constructionYear, -D, -B, -C, -y)
GGally::ggpairs(data, lower = list(continuous = "smooth_loess"))
#--> xx ~ xx has r > 0.7 (Dormann et al. 2013 Ecography) --> xx has to be excluded
rm(data)

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
m1 <- lmer(y ~ (comparison + exposition + PC1soil + exp(PC2soil))^2 + PC3soil + side + locationAbb + ageGroup + log(distanceRiver) +
             (1|plot), 
           REML = F,
           data = tbi);simulateResiduals(m1, plot = T)
m2 <- lmer(y ~ (comparison + exposition + PC1soil)^2 + (PC2soil) + PC3soil + side + constructionYearF + locationAbb + log(distanceRiver) +
             (1|plot), 
           REML = F, 
           data = tbi)
simulateResiduals(m2, plot = T)
m3 <- lmer(y ~ (comparison + exposition + exp(PC2soil))^2 + PC1soil + PC3soil + side + constructionYearF + locationAbb + log(distanceRiver) +
             (1|plot), 
           REML = F, 
           data = tbi)
simulateResiduals(m3, plot = T)
m4 <- lmer(y ~ comparison + exposition * PC1soil + (PC2soil) + PC3soil + side + constructionYearF + locationAbb + log(distanceRiver) +
             (1|plot), 
           REML = F, 
           data = tbi)
simulateResiduals(m4, plot = T)
m5 <- lmer(y ~ comparison + exposition * (PC2soil) + PC1soil + PC3soil + side + constructionYearF + locationAbb + log(distanceRiver) +
             (1|plot), 
           REML = F,
           data = tbi)
simulateResiduals(m5, plot = T)
m6 <- lmer(y ~ comparison * exposition + PC1soil + (PC2soil) + PC3soil + side + constructionYearF + locationAbb + log(distanceRiver) +
             (1|plot), 
           REML = F,
           data = tbi)
simulateResiduals(m6, plot = T)
m7 <- lmer(y ~ comparison * PC1soil + exposition + (PC2soil) + PC3soil + side + constructionYearF + locationAbb + log(distanceRiver) +
             (1|plot), 
           REML = F,
           data = tbi)
simulateResiduals(m7, plot = T)
m8 <- lmer(y ~ comparison + (exposition + PC1soil + (PC2soil) + PC3soil + side) + log(distanceRiver) + locationYear + 
             comparison:exposition +
             (1|plot), 
           REML = T,
           data = tbi)
simulateResiduals(m8, plot = T)

### b comparison -----------------------------------------------------------------------------------------

AICcmodavg::aictab(cand.set = list("m1a" = m1a, "m1b" = m1b, "m1c" = m1c))
AICcmodavg::aictab(cand.set = list("m1" = m1, "m2" = m2, "m3" = m3, "m4" = m4, "m5" = m5, "m6" = m6, "m7" = m7, "m8" = m8))
car::Anova(m8, type = 3)
sjPlot::plot_model(m8, type = "emm", terms = c("exposition", "comparison"), show.data = F, jitter = .7)
ggsave(here("outputs/figures/figure_tbi_bc_presence_exposition_comparison_(800dpi_9x10cm).tiff"), dpi = 800, width = 9, height = 10, units = "cm")
sjPlot::plot_model(m8)
ggsave(here("outputs/figures/figure_tbi_bc_presence_(800dpi_16x10cm).tiff"), dpi = 800, width = 16, height = 10, units = "cm")
dotwhisker::dwplot(list(m8, m4, m7), show_intercept = T)
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
car::vif(m2) # all < 3 (Zuur et al. 2010 Methods Ecol Evol)


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
MuMIn::r.squaredGLMM(m2) #R2m = 0.373, R2c = 0.373
VarCorr(m2)
sjPlot::plot_model(m2, type = "re", show.values = T)
car::Anova(m2, type = 2)
anova(m2)
### * Effect sizes ####
(emm <- emmeans(m8, revpairwise ~ comparison | exposition, type = "response"))
plot(emm, comparison = T)
(emm <- emmeans(m4, revpairwise ~ side, type = "response"))

### * Save ####
table <- broom::tidy(car::Anova(m2, type = 3))
write.csv(table, here("outputs/statistics/table_anova_tbi_bc_presence.csv"))

