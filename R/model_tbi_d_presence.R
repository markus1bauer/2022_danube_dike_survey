# Beta diversity on dike grasslands
# Total temporal beta diversity (TBI) (presence-absence data) ####
# Markus Bauer
# 2022-01-11



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Preparation #########################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


### Packages ###
library(here)
library(tidyverse)
library(ggbeeswarm)
library(blme)
library(DHARMa)
library(emmeans)

### Start ###
rm(list = ls())
setwd(here("data", "processed"))

### Load data ###
sites <- read_csv("data_processed_sites_temporal.csv",
  col_names = TRUE,
  na = c("", "na", "NA"), col_types =
    cols(
      .default = "?",
      plot = "f",
      block = "f",
      comparison = "f",
      exposition = "f",
      side = "f",
      locationYear = "f"
    )
) %>%
  filter(comparison == "1718" | comparison == "1819" | comparison == "1921") %>%
  mutate(
    y = D_presence,
    comparison = factor(comparison)
  )

data_collinearity <- sites %>%
  select(where(is.numeric))

sites <- sites %>%
  mutate(across(c("longitude", "latitude", "riverkm", "distanceRiver"), scale))



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Statistics ##########################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


## 1 Data exploration ###################################################

#### * Graphs #####
mean(sites$y)
median(sites$y)
sd(sites$y)
Rmisc::CI(sites$y, ci = .95)
quantile(sites$y, probs = c(0.05, 0.95), na.rm = TRUE)
# main
ggplot(sites, aes(x = D_abundance, y = y)) +
  geom_point() +
  geom_smooth(method = "loess")
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
ggplot(sites, aes(x = (distanceRiver), y = log(y))) +
  geom_point() +
  geom_smooth(method = "lm")
ggplot(sites, aes(x = constructionYear, y = y)) +
  geom_point() +
  geom_smooth(method = "loess")
ggplot(sites, aes(x = PC1soil, y = (y))) +
  geom_point() +
  geom_smooth(method = "lm")
ggplot(sites, aes(x = exp(PC2soil), y = y)) +
  geom_point() +
  geom_smooth(method = "lm")
# 2way
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
dotchart(log(sites$y), groups = factor(sites$exposition),
         main = "Cleveland dotplot")
sites %>% count(locationYear)
sites %>%
  count(plot) %>%
  count(n)
boxplot(sites$y)
# identify(rep(1, length(data$rgr13)), data$rgr13, labels = c(data$n))
plot(table((sites$y)), type = "h", xlab = "Observed values", ylab = "Frequency")
ggplot(sites, aes(y)) +
  geom_density()
ggplot(sites, aes(log(y))) +
  geom_density()

#### * check collinearity ####
GGally::ggpairs(data_collinearity, lower = list(continuous = "smooth_loess"))
#--> riverkm ~ longitude/latitude has r > 0.7 (Dormann et al. 2013 Ecography)
rm(data_collinearity)


## 2 Model building #####################################################

### a models ------------------------------------------------------------

### * random structure ####
m1a <- blmer(log(y) ~ 1 + (1 | locationYear), data = sites, REML = TRUE)
m1b <- blmer(log(y) ~ 1 + (1 | locationYear / plot), data = sites, REML = TRUE)
m1c <- blmer(log(y) ~ 1 + (1 | plot), data = sites, REML = TRUE)
MuMIn::AICc(m1a, m1b, m1c) # m1c most parsimonous

#### * fixed effects ####
m1 <- blmer(log(y) ~
              (comparison + exposition + PC1soil)^2 + PC2soil + PC3soil +
              side + distanceRiver + locationYear + D_abundance +
              (1 | plot),
REML = FALSE,
control = lmerControl(optimizer = "Nelder_Mead"),
cov.prior = wishart,
data = sites
)
simulateResiduals(m1, plot = TRUE)
m2 <- blmer(log(y) ~ comparison + exposition * PC1soil + PC2soil + PC3soil +
              side + distanceRiver + locationYear + D_abundance +
              (1 | plot),
            REML = FALSE,
            control = lmerControl(optimizer = "Nelder_Mead"),
            cov.prior = wishart,
            data = sites)
simulateResiduals(m2, plot = TRUE)
m3 <- blmer(log(y) ~ comparison * exposition + PC1soil + PC2soil + PC3soil + side +
  distanceRiver + locationYear + D_abundance +
  (1 | plot),
  REML = FALSE,
  control = lmerControl(optimizer = "Nelder_Mead"),
  cov.prior = wishart,
  data = sites
)
simulateResiduals(m3, plot = TRUE)
m4 <- blmer(log(y) ~ comparison * PC1soil + exposition + PC2soil + PC3soil + side +
  distanceRiver + locationYear + D_abundance +
  (1 | plot),
  REML = FALSE,
  control = lmerControl(optimizer = "Nelder_Mead"),
  cov.prior = wishart,
  data = sites
)
simulateResiduals(m4, plot = TRUE)
m5 <- blmer(log(y) ~ comparison + exposition + PC1soil + PC2soil + PC3soil + side +
  distanceRiver + locationYear + D_abundance +
  (1 | plot),
  REML = FALSE,
  control = lmerControl(optimizer = "Nelder_Mead"),
  cov.prior = wishart,
  data = sites
)
simulateResiduals(m5, plot = TRUE)

### b comparison --------------------------------------------------------

MuMIn::AICc(m1, m2, m3, m4, m5) # m2 most parsimonious; Use AICc and not AIC since ratio n/K < 40 (Burnahm & Anderson 2002 p. 66)
dotwhisker::dwplot(list(m2, m5),
  show_intercept = FALSE,
  vline = geom_vline(
    xintercept = 0,
    colour = "grey60",
    linetype = 2
  )
) +
  theme_classic()
m <- update(m2, REML = TRUE)
rm(list = setdiff(ls(), c("sites", "m")))

### c model check -------------------------------------------------------

simulationOutput <- simulateResiduals(m, plot = TRUE)
plotResiduals(simulationOutput$scaledResiduals, sites$locationYear)
plotResiduals(simulationOutput$scaledResiduals, sites$plot)
plotResiduals(simulationOutput$scaledResiduals, sites$comparison)
plotResiduals(simulationOutput$scaledResiduals, sites$exposition)
plotResiduals(simulationOutput$scaledResiduals, sites$side)
plotResiduals(simulationOutput$scaledResiduals, sites$PC1soil)
plotResiduals(simulationOutput$scaledResiduals, sites$PC2soil)
plotResiduals(simulationOutput$scaledResiduals, sites$PC3soil)
plotResiduals(simulationOutput$scaledResiduals, sites$distanceRiver)
plotResiduals(simulationOutput$scaledResiduals, sites$D_abundance)
plotResiduals(simulationOutput$scaledResiduals, sites$riverkm)
car::vif(m) # remove riverkm since > 3 oder 10 (Zuur et al. 2010 Methods Ecol Evol)


## 3 Chosen model output ################################################

### * Model output ####
MuMIn::r.squaredGLMM(m) # R2m = 0.375, R2c = 0.534
VarCorr(m)
sjPlot::plot_model(m, type = "re", show.values = TRUE)
dotwhisker::dwplot(m,
  show_intercept = FALSE,
  vline = geom_vline(
    xintercept = 0,
    colour = "grey60",
    linetype = 2
  )
) +
  theme_classic()

### * Effect sizes ####
(emm <- emmeans(m, revpairwise ~ side, type = "response"))
plot(emm, comparison = TRUE)
(emm <- emmeans(m, revpairwise ~ comparison, type = "response"))
sjPlot::plot_model(m, type = "emm", terms = c("PC1soil", "exposition"),
                   show.data = TRUE)
