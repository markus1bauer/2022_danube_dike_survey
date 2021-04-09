# Model for species richness ####
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
                       changeType = col_factor(c("FFH6510", "better", "change", "worse", "any-FFH", "non-FFH"))
                     )) %>%
  select(fdisAbuSla, fdisAbuSeedmass, fdisAbuHeight, surveyYear, constructionYear, id, block, location, plot, exposition, plotAge, PC1, PC2, PC3, changeType, ffh) %>%
  pivot_longer(c("fdisAbuSla", "fdisAbuSeedmass", "fdisAbuHeight"), names_to = "type", values_to = "n") %>%
  mutate(plotAge = scale(plotAge, scale = T, center = T)) %>%
  mutate(surveyYear = scale(surveyYear, scale = T, center = T)) %>%
  mutate(surveyYearF = as_factor(surveyYear)) %>%
  mutate(constructionYearF = as_factor(constructionYear)) %>%
  filter(ffh != "6210", 
         changeType != "any-FFH", 
         !is.na(changeType),
         exposition != "east",
         exposition != "west") %>%
  mutate(changeType = fct_collapse(changeType, "change/worse" = c("change", "worse")))



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Statistics ################################################################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


### 1 Data exploration #####################################################################################

#### a Graphs ---------------------------------------------------------------------------------------------
#simple effects:
ggplot(sites, aes(y = n, x = type)) + geom_boxplot() 
#2way
ggplot(sites, aes(x = constructionYear, y = n, color = type)) + 
  geom_point() + 
  geom_smooth(method = "lm") + 
  facet_wrap(~type)
ggplot(sites, aes(x = type, y = n, color = exposition)) + 
  geom_boxplot() + 
  geom_quasirandom(dodge.width = .7, groupOnX = T)
ggplot(sites, aes(x = PC1, y = n, color = type)) + 
  geom_point() + geom_smooth(method = "lm") + 
  facet_wrap(~type)
ggplot(sites, aes(x = PC2, y = n, color = type)) +
  geom_point() + 
  geom_smooth(method = "lm") + 
  facet_wrap(~type)
ggplot(sites, aes(x = PC3, y = n, color = type)) + 
  geom_smooth(aes(group = plot), colour = "grey40", method = "lm", se = F, show.legend = F) + 
  geom_point() + 
  geom_smooth(method = "lm") + 
  facet_wrap(~type)
ggplot(sites, aes(x = plotAge, y = n, color = type)) + 
  geom_smooth(aes(group = plot), colour = "grey40", method = "lm", se = F, show.legend = F) + 
  geom_point() + 
  geom_smooth(method = "lm") + 
  facet_wrap(~type)
ggplot(sites, aes(x = type, y = n, color = changeType)) + 
  geom_boxplot() + 
  geom_quasirandom(dodge.width = .7, groupOnX = T)
ggplot(sites, aes(x = type, y = n, color = ffh)) + 
  geom_boxplot() + 
  geom_quasirandom(dodge.width = .7, groupOnX = T)



##### b Outliers, zero-inflation, transformations? -----------------------------------------------------
par(mfrow = c(2,2))
dotchart((sites$n), groups = factor(sites$changeType), main = "Cleveland dotplot")
par(mfrow=c(1,1));
boxplot(sites$n);#identify(rep(1, length(edata$rgr13)), edata$rgr13, labels = c(edata$n))
plot(table((sites$n)), type = "h", xlab = "Observed values", ylab = "Frequency")
ggplot(sites, aes(n)) + geom_density()
ggplot(sites, aes(sqrt(n))) + geom_density()


## 2 Model building ################################################################################

#### a models ----------------------------------------------------------------------------------------
#random structure
m1a <- lmer(n ~ 1 + (surveyYear|location/plot), sites, REML = F)
isSingular(m1a)
m1b <- lmer(n ~ 1 + (surveyYear|plot), sites, REML = F)
isSingular(m1b)
m1c <- lmer(n ~ 1 + (1|location/plot), sites, REML = F)
isSingular(m1c)
m1d <- lmer(n ~ 1 + (1|plot), sites, REML = F)
isSingular(m1d)
VarCorr(m1a)
VarCorr(m1b)
VarCorr(m1c)
VarCorr(m1d)
#1w-model
m2 <- lmer(n ~ type * (exposition + PC1 + PC2 + PC3 + ffh + changeType) + 
             (surveyYear|location/plot), sites, REML = F)
simulateResiduals(m2, plot = T)
isSingular(m2)
#2w-model
m3 <- lmer(n ~ (exposition + PC1 + PC2 + PC3 + ffh + changeType)^2 +
             (surveyYear|location/plot), sites, REML = F)
isSingular(m3)
simulateResiduals(m3, plot = T)

#### b comparison -----------------------------------------------------------------------------------------
anova(m2,m3)
rm(m3)

#### c model check -----------------------------------------------------------------------------------------
simulationOutput <- simulateResiduals(m2, plot = T)
par(mfrow=c(2,2));
plotResiduals(main = "year", simulationOutput$scaledResiduals, sites$year)
plotResiduals(main = "type", simulationOutput$scaledResiduals, sites$type)
plotResiduals(main = "block", simulationOutput$scaledResiduals, sites$block)


## 3 Chosen model output ################################################################################

### Model output ---------------------------------------------------------------------------------------------
MuMIn::r.squaredGLMM(m2)
0.060 / 0.248
VarCorr(m2)
sjPlot::plot_model(m2, type = "re", show.values = T)
car::Anova(m2, type = 2)

### Effect sizes -----------------------------------------------------------------------------------------
(emm <- emmeans(m2, revpairwise ~ changeType, type = "response"))
plot(emm, comparison = T)
contrast(emmeans(m2, ~ changeType, type = "response"), "trt.vs.ctrl", ref = 4)
(emm <- emmeans(m2, revpairwise ~ exposition, type = "response"))

### Save ###
table <- tidy(car::Anova(m2, type = 2))
setwd(here("data/tables"))
write.csv2(table, "table_anova_fdisAbuSingle.csv")
