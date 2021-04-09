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
                       .default = col_guess(),
                       id = "f",
                       location = "f",
                       block = "f",
                       plot = "f",
                       exposition = "f",
                       side = "f",
                       ffh = "f",
                       changeType = col_factor(c("FFH6510", "any-FFH", "better", "change", "worse", "non-FFH"))
                     )) %>%
  select(targetRichness, id, surveyYear, constructionYear, plotAge, location, block, plot, side, exposition, PC1, PC2, PC3, changeType, ffh) %>%
  mutate(n = targetRichness) %>%
  mutate(location = factor(location, levels = unique(location[order(constructionYear)]))) %>%
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
  mutate(exposition = factor(exposition)) %>%
  mutate(location = factor(location)) %>%
  mutate(plot = factor(plot)) %>%
  mutate(block = factor(block)) %>%
  select(-count)



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Statistics ################################################################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


### 1 Data exploration #####################################################################################

#### a Graphs ---------------------------------------------------------------------------------------------
#simple effects:
ggplot(sites, aes(y = n, x = location)) + geom_boxplot() + geom_quasirandom()
ggplot(sites, aes(y = n, x = constructionYear)) + geom_point() + geom_smooth(method = "lm") 
ggplot(sites, aes(y = n, x = constructionYearF)) + geom_boxplot() + geom_quasirandom()
ggplot(sites, aes(y = n, x = plotAge)) + geom_point() + geom_smooth(method = "lm") 
ggplot(sites, aes(y = n, x = surveyYearF)) + geom_boxplot() + geom_quasirandom()
ggplot(sites, aes(y = n, x = exposition)) + geom_boxplot() + geom_quasirandom()
ggplot(sites, aes(y = n, x = PC1)) + geom_point() + geom_smooth(method = "lm")
ggplot(sites, aes(y = n, x = PC2)) + geom_point() + geom_smooth(method = "lm")
ggplot(sites, aes(y = n, x = PC3)) + geom_point() + geom_smooth(method = "lm")
ggplot(sites, aes(y = n, x = changeType)) + geom_boxplot() + geom_quasirandom()
ggplot(sites, aes(y = n, x = ffh)) + geom_boxplot() + geom_quasirandom()
sites %>%
  group_by(changeType) %>%
  summarise(mean = mean(n), sd = sd(n))
#2way
ggplot(sites, aes(x = constructionYear, y = n, color = changeType)) + 
  geom_point() + 
  geom_smooth(method = "lm") + 
  facet_wrap(~changeType)
ggplot(sites, aes(x = plotAge, y = n, color = changeType)) + 
  geom_smooth(aes(group = plot), colour = "grey40", method = "lm", se = F, show.legend = F) + 
  geom_point() + 
  geom_smooth(method = "lm") + 
  facet_wrap(~changeType)


##### b Outliers, zero-inflation, transformations? -----------------------------------------------------
dotchart((sites$n), groups = factor(sites$exposition), main = "Cleveland dotplot")
boxplot(sites$n);#identify(rep(1, length(edata$rgr13)), edata$rgr13, labels = c(edata$n))
plot(table((sites$n)), type = "h", xlab = "Observed values", ylab = "Frequency")
ggplot(sites, aes(n)) + geom_density()
ggplot(sites, aes(sqrt(n))) + geom_density()


## 2 Model building ################################################################################

#### a models ----------------------------------------------------------------------------------------
#random structure
m1a <- lmer(n ~ 1 + (surveyYear|location/plot), sites, REML = F)
VarCorr(m1a) # convergence problems
m1b <- lmer(n ~ 1 + (surveyYear|plot), sites, REML = F)
VarCorr(m1b) # convergence problems
m1c <- lmer(n ~ 1 + (1|location/plot), sites, REML = F)
VarCorr(m1c)
m1d <- lmer(n ~ 1 + (1|plot), sites, REML = F)
VarCorr(m1d)
#fixed effects
m2 <- lmer((n) ~ (exposition + side + PC1 + PC2 + PC3) + 
             (1|plot), sites, REML = F)
simulateResiduals(m2, plot = T)
#fixed and site and year effects
m3 <- lmer((n) ~ (exposition + side + PC1 + PC2 + PC3 + surveyYearF + location) +
             (1|plot), sites, REML = F);
simulateResiduals(m3, plot = T)
isSingular(m3)
#plotAge instead of location
m4 <- lmer((n) ~ (exposition + side + PC1 + PC2 + PC3 + surveyYearF + plotAge) + 
             (1|plot), sites, REML = F);
simulateResiduals(m4, plot = T)
isSingular(m4)

#### b comparison -----------------------------------------------------------------------------------------
anova(m2, m3, m4)
rm(m1a, m1b, m1c, m1d, m4)

#### c model check -----------------------------------------------------------------------------------------
simulationOutput <- simulateResiduals(m3, plot = F)
plotResiduals(simulationOutput$scaledResiduals, sites$surveyYearF)
plotResiduals(simulationOutput$scaledResiduals, sites$location)
plotResiduals(simulationOutput$scaledResiduals, sites$block)
plotResiduals(simulationOutput$scaledResiduals, sites$plot)
plotResiduals(simulationOutput$scaledResiduals, sites$side)
plotResiduals(simulationOutput$scaledResiduals, sites$exposition)
plotResiduals(simulationOutput$scaledResiduals, sites$PC1)
plotResiduals(simulationOutput$scaledResiduals, sites$PC2)
plotResiduals(simulationOutput$scaledResiduals, sites$PC3)


## 3 Chosen model output ################################################################################

### Model output ---------------------------------------------------------------------------------------------
MuMIn::r.squaredGLMM(m3)
0.544 /  0.712
MuMIn::r.squaredGLMM(m2) # to see which R2m is due to 'location' and 'surveyYearF'
0.265 / 0.661
rm(m2)
VarCorr(m3)
sjPlot::plot_model(m3, type = "re", show.values = T)
car::Anova(m3, type = 2)

### Effect sizes -----------------------------------------------------------------------------------------
(emm <- emmeans(m3, revpairwise ~ exposition, type = "response"))
plot(emm, comparison = T)
(emm <- emmeans(m3, revpairwise ~ side, type = "response"))

### Save ###
table <- tidy(car::Anova(m3, type = 2))
setwd(here("data/tables"))
write.csv2(table, "table_anova_targetRichness.csv")
