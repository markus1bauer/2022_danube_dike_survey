# Model for PERMDISP ####
# Markus Bauer



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Preparation ################################################################################################################
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
sites <- read_csv2("data_processed_sites.csv", col_names = T, na = "na", col_types = 
                     cols(
                       .default = col_guess(),
                       id = col_factor(),
                       location = col_factor(),
                       block = col_factor(),
                       plot = col_factor(),
                       exposition = col_factor(),
                       ffh = col_factor(),
                       changeType = col_factor(c("better", "change", "worse", "any-FFH", "FFH6510", "non-FFH"))
                     )) %>%
  select(permdisp, surveyYear, constructionYear, id, block, location, plot, exposition, plotAge, PC1, PC2, PC3, changeType, ffh, graminoidCov) %>%
  mutate(n = permdisp) %>%
  mutate(plotAge = scale(plotAge, scale = T, center = T)) %>%
  mutate(graminoidCov = scale(graminoidCov, scale = T, center = T)) %>%
  mutate(surveyYear = scale(surveyYear, scale = T, center = T)) %>%
  mutate(surveyYearF = as_factor(surveyYear)) %>%
  mutate(constructionYearF = as_factor(constructionYear)) %>%
  filter(ffh != "6210", changeType != "any-FFH", !is.na(changeType))



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Statistics ################################################################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


### 1 Data exploration #####################################################################################

#### a Graphs ---------------------------------------------------------------------------------------------
#simple effects:
ggplot(sites, aes(y = n, x = location)) + geom_boxplot() 
ggplot(sites, aes(y = n, x = surveyYearF)) + geom_boxplot() + geom_quasirandom()
ggplot(sites, aes(y = n, x = constructionYearF)) + geom_boxplot() + geom_quasirandom()
ggplot(sites, aes(y = n, x = constructionYear)) + geom_point() + geom_smooth(method = "loess") 
ggplot(sites, aes(y = n, x = plotAge)) + geom_point() + geom_smooth(method = "lm") 
ggplot(sites, aes(y = n, x = exposition)) + geom_boxplot() + geom_quasirandom()
ggplot(sites, aes(y = n, x = PC1)) + geom_point() + geom_smooth(method = "lm")
ggplot(sites, aes(y = n, x = PC2)) + geom_point() + geom_smooth(method = "lm")
ggplot(sites, aes(y = n, x = PC3)) + geom_point() + geom_smooth(method = "lm")
ggplot(sites, aes(y = n, x = graminoidCov)) + geom_point() + geom_smooth(method = "lm")
ggplot(sites, aes(y = n, x = changeType)) + geom_boxplot() + geom_quasirandom()
ggplot(sites, aes(y = n, x = ffh)) + geom_boxplot() + geom_quasirandom()
#2way
ggplot(sites, aes(x = PC1, y = n, color = changeType)) + 
  geom_point() + geom_smooth(method = "lm") + 
  facet_wrap(~changeType)
ggplot(sites, aes(x = PC2, y = n, color = changeType)) +
  geom_point() + 
  geom_smooth(method = "lm") + 
  facet_wrap(~changeType)
ggplot(sites, aes(x = PC3, y = n, color = changeType)) + 
  geom_smooth(aes(group = plot), colour = "grey40", method = "lm", se = F, show.legend = F) + 
  geom_point() + 
  geom_smooth(method = "lm") + 
  facet_wrap(~changeType)
ggplot(sites, aes(x = plotAge, y = n, color = changeType)) + 
  geom_smooth(aes(group = plot), colour = "grey40", method = "lm", se = F, show.legend = F) + 
  geom_point() + 
  geom_smooth(method = "lm") + 
  facet_wrap(~changeType)
ggplot(sites, aes(x = exposition, y = n, color = surveyYearF)) + 
  geom_boxplot()
ggplot(sites, aes(x = plotAge, y = n, color = surveyYearF)) + 
  geom_point() + 
  geom_smooth(method = "lm") + 
  facet_wrap(~surveyYearF)
ggplot(sites, aes(x = PC1, y = n, color = surveyYearF)) + 
  geom_point() + 
  geom_smooth(method = "lm") + 
  facet_wrap(~surveyYearF)
ggplot(sites, aes(x = PC2, y = n, color = surveyYearF)) + 
  geom_point() + 
  geom_smooth(method = "lm") + 
  facet_wrap(~surveyYearF)
ggplot(sites, aes(x = PC3, y = n, color = surveyYearF)) + 
  geom_point() + geom_smooth(method = "lm") + 
  facet_wrap(~surveyYearF)


##### b Outliers, zero-inflation, transformations? -----------------------------------------------------
par(mfrow = c(2,2))
dotchart((sites$n), groups = factor(sites$plotAge), main = "Cleveland dotplot")
par(mfrow=c(1,1));
boxplot(sites$n);#identify(rep(1, length(edata$rgr13)), edata$rgr13, labels = c(edata$n))
plot(table((sites$n)), type = "h", xlab = "Observed values", ylab = "Frequency")
ggplot(sites, aes(n)) + geom_density()
ggplot(sites, aes(sqrt(n))) + geom_density()


## 2 Model building ################################################################################

#### a models ----------------------------------------------------------------------------------------
#random structure
m1a <- lmer(n ~ 1 + (surveyYear|location/block/plot), sites, REML = F)
m1b <- lmer(n ~ 1 + (surveyYear|block/plot), sites, REML = F)
m1c <- lmer(n ~ 1 + (surveyYear|location/plot), sites, REML = F)
m1d <- lmer(n ~ 1 + (surveyYear|plot), sites, REML = F)
m1e <- lmer(n ~ 1 + (1|location/plot), sites, REML = F)
VarCorr(m1d)
VarCorr(m1e)
anova(m1a, m1b, m1c, m1d, m1e)
#2w-model
m2 <- lmer((n) ~ (location + surveyYearF + exposition + PC1 + PC2 + PC3 + changeType) + 
             (surveyYear|plot), sites, REML = F)
simulateResiduals(m2, plot = T)
isSingular(m2)
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
plotResiduals(main = "surveyYearF", simulationOutput$scaledResiduals, sites$surveyYearF)
plotResiduals(main = "location", simulationOutput$scaledResiduals, sites$location)
plotResiduals(main = "plot", simulationOutput$scaledResiduals, sites$plot)
plotResiduals(main = "plotAge", simulationOutput$scaledResiduals, sites$plotAge)
plotResiduals(main = "exposition", simulationOutput$scaledResiduals, sites$exposition)
plotResiduals(main = "PC1", simulationOutput$scaledResiduals, sites$PC1)
plotResiduals(main = "PC2", simulationOutput$scaledResiduals, sites$PC2)
plotResiduals(main = "PC3", simulationOutput$scaledResiduals, sites$PC3)
plotResiduals(main = "changeType", simulationOutput$scaledResiduals, sites$changeType)


## 3 Chosen model output ################################################################################

### Model output ---------------------------------------------------------------------------------------------
#m2
MuMIn::r.squaredGLMM(m2) #R2m = 0.445, R2c = 0.527 2021-03-30
VarCorr(m2)
sjPlot::plot_model(m2, type = "re", show.values = T)
car::Anova(m2, type = 3)

### Effect sizes -----------------------------------------------------------------------------------------
(emm <- emmeans(m2, revpairwise ~ changeType, type = "response"))
(emm <- emmeans(m2, revpairwise ~ exposition, type = "response"))
plot(emm, comparison = T)
contrast(emmeans(m2, ~ changeType, type = "response"), "trt.vs.ctrl", ref = 4)
contrast(emmeans(m2, ~ exposition, type = "response"), "trt.vs.ctrl", ref = 1)

### Save ###
table <- tidy(car::Anova(m2, type = 3))
setwd(here("data/tables"))
write.csv2(table, "table_anova_permdisp.csv")
