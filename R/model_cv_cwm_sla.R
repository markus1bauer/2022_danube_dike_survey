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
sites <- read_csv2("data_processed_sites.csv", col_names = T, na = c("na", "NA"), col_types = 
                     cols(
                       .default = "?",
                       id = "f",
                       locationAbb = "f",
                       block = "f",
                       plot = "f",
                       exposition = "f",
                       side = "f",
                       ffh = "f",
                       changeType = col_factor(c("FFH6510", "any-FFH", "better", "change", "worse", "non-FFH"))
                     )) %>%
  select(cwmAbuSla, id, surveyYear, constructionYear, plotAge, locationAbb, block, plot, side, exposition, PC1, PC2, PC3, changeType, ffh) %>%
  mutate(n = cwmAbuSla) %>%
  mutate(locationAbb = factor(locationAbb, levels = unique(locationAbb[order(constructionYear)]))) %>%
  mutate(plotAge = scale(plotAge, scale = T, center = F)) %>%
  mutate(surveyYear = scale(surveyYear, scale = T, center = F)) %>%
  mutate(surveyYearF = as_factor(surveyYear)) %>%
  mutate(constructionYear = scale(constructionYear, scale = T, center = F)) %>%
  mutate(constructionYearF = as_factor(constructionYear)) %>%
  filter(exposition != "east",
         exposition != "west") %>%
  group_by(plot) %>%
  mutate(count = n()) %>%
  ungroup() %>%
  slice_max(count, n = 1) %>%
  select(-count) %>%
  mutate(exposition = factor(exposition)) %>%
  mutate(locationAbb = factor(locationAbb)) %>%
  mutate(locationAbb = factor(locationAbb, levels = unique(locationAbb[order(constructionYear)]))) %>%
  mutate(plot = factor(plot)) %>%
  mutate(block = factor(block)) %>%
  group_by(plot, exposition, side, changeType, PC1, PC2, PC3, locationAbb, constructionYearF, constructionYear) %>%
  summarise(mean = mean(n), sd = sd(n)) %>%
  mutate(n = sd / mean) #Calculate Coefficient of variance (CV)



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Statistics ################################################################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


### 1 Data exploration #####################################################################################

#### a Graphs ---------------------------------------------------------------------------------------------
#simple effects:
ggplot(sites, aes(y = n, x = locationAbb)) + geom_boxplot() + geom_quasirandom()
ggplot(sites, aes(y = n, x = constructionYear)) + geom_point() + geom_smooth(method = "lm") 
ggplot(sites, aes(y = n, x = constructionYearF)) + geom_boxplot() + geom_quasirandom()
ggplot(sites, aes(y = n, x = exposition)) + geom_boxplot() + geom_quasirandom()
ggplot(sites, aes(y = n, x = side)) + geom_boxplot() + geom_quasirandom()
ggplot(sites, aes(y = n, x = PC1)) + geom_point() + geom_smooth(method = "lm")
ggplot(sites, aes(y = n, x = PC2)) + geom_point() + geom_smooth(method = "lm")
ggplot(sites, aes(y = n, x = PC3)) + geom_point() + geom_smooth(method = "lm")
ggplot(sites, aes(y = n, x = changeType)) + geom_boxplot() + geom_quasirandom()
sites %>%
  group_by(changeType) %>%
  summarise(mean = mean(n), sd = sd(n))
#2way
ggplot(sites, aes(x = constructionYearF, y = n, color = exposition)) + 
  geom_boxplot() +
  geom_quasirandom()

##### b Outliers, zero-inflation, transformations? -----------------------------------------------------
dotchart((sites$n), groups = factor(sites$exposition), main = "Cleveland dotplot")
boxplot(sites$n);#identify(rep(1, length(edata$rgr13)), edata$rgr13, labels = c(edata$n))
plot(table((sites$n)), type = "h", xlab = "Observed values", ylab = "Frequency")
ggplot(sites, aes(n)) + geom_density()
ggplot(sites, aes(sqrt(n))) + geom_density()


## 2 Model building ################################################################################

#### a models ----------------------------------------------------------------------------------------
#1w-model
m1 <- lm(sqrt(n) ~ (exposition + side + PC1 + PC2 + PC3), sites)
simulateResiduals(m1, plot = T)
summary(m1)

#### b comparison -----------------------------------------------------------------------------------------
anova(m2,m3)
rm(m3)

#### c model check -----------------------------------------------------------------------------------------
simulationOutput <- simulateResiduals(m1, plot = T)
par(mfrow=c(2,2));
plotResiduals(simulationOutput$scaledResiduals, sites$locationAbb)
plotResiduals(simulationOutput$scaledResiduals, sites$side)
plotResiduals(simulationOutput$scaledResiduals, sites$exposition)
plotResiduals(simulationOutput$scaledResiduals, sites$PC1)
plotResiduals(simulationOutput$scaledResiduals, sites$PC2)
plotResiduals(simulationOutput$scaledResiduals, sites$PC3)


## 3 Chosen model output ################################################################################

### Model output ---------------------------------------------------------------------------------------------
MuMIn::r.squaredGLMM(m1)
0.095 / 0.332
VarCorr(m2)
sjPlot::plot_model(m2, type = "re", show.values = T)
car::Anova(m1, type = 2)

### Effect sizes -----------------------------------------------------------------------------------------
(emm <- emmeans(m2, revpairwise ~ changeType, type = "response"))
plot(emm, comparison = T)
contrast(emmeans(m2, ~ changeType, type = "response"), "trt.vs.ctrl", ref = 1)
(emm <- emmeans(m2, revpairwise ~ ffh, type = "response"))

### Save ###
table <- tidy(car::Anova(m2, type = 2))
setwd(here("data/tables"))
write.csv2(table, "table_anova_cwmAbuSla.csv")
