# Beta diversity on dike grasslands
# Ratio of gains and losses of TBI - All species ####
# Markus Bauer
# 2023-01-16



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Preparation ###############################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



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
                  col_names = TRUE, na = c("", "na", "NA"),
                  col_types =
                    cols(
                      .default = "?",
                      plot = "f",
                      block = "f",
                      comparison = "f",
                      location = "f",
                      location_construction_year = "f",
                      exposition = col_factor(levels = c("south", "north")),
                      orientation = col_factor(levels = c("land", "water"))
                    )) %>%
  filter(
    (comparison == "1718" | comparison == "1819" | comparison == "1921") &
      pool == "all" & presabu == "presence"
    ) %>%
  mutate(
    y = c - b,
    comparison = factor(comparison)
    ) %>%
  mutate(across(
    c("river_km", "river_distance", "biotope_distance", "biotope_area"), scale
    ))



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Statistics ################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



## 1 Data exploration #########################################################


### a Graphs ------------------------------------------------------------------

mean(sites$y)
median(sites$y)
sd(sites$y)
Rmisc::CI(sites$y, ci = .95)
quantile(sites$y, probs = c(0.05, 0.95), na.rm = TRUE)
ggplot(sites, aes(x = comparison, y = y)) +
  geom_quasirandom(color = "grey") + geom_boxplot(fill = "transparent") +
  labs(title = "Comparison of consecutive surveys")
ggplot(sites, aes(x = exposition, y = y)) +
  geom_quasirandom(color = "grey") + geom_boxplot(fill = "transparent") +
  labs(title = "Exposition of dike slopes")
ggplot(sites, aes(x = orientation, y = y)) +
  geom_quasirandom(color = "grey") + geom_boxplot(fill = "transparent") +
  labs(title = "Orientation of dike slopes")
ggplot(sites, aes(x = river_km, y = (y))) +
  geom_point() +  geom_smooth(method = "lm") +
  labs(title = "Position along the river")
ggplot(sites, aes(x = location_construction_year, y = y)) +
  geom_quasirandom(color = "grey") + geom_boxplot(fill = "transparent") +
  labs(title = "Location and construction year of the dike")
ggplot(sites, aes(x = (river_distance), y = (y))) +
  geom_point() + geom_smooth(method = "lm") +
  labs(title = "Distance to river course")
ggplot(sites, aes(x = (biotope_distance), y = (y))) +
  geom_point() + geom_smooth(method = "lm") +
  labs(title = "Distance to closest grassland biotope")
ggplot(sites, aes(x = (biotope_area), y = (y))) +
  geom_point() + geom_smooth(method = "lm") +
  labs(title = "Amount of grassland biotopes with 500 m radius")
ggplot(sites, aes(x = pc1_soil, y = (y))) +
  geom_point() + geom_smooth(method = "lm") +
  labs(title = "PC1 (soil)")
ggplot(sites, aes(x = (pc2_soil), y = y)) +
  geom_point() + geom_smooth(method = "lm") +
  labs(title = "PC2 (soil)")
ggplot(sites, aes(x = comparison, y = y)) +
  geom_quasirandom(color = "grey") + geom_boxplot(fill = "transparent") +
  facet_grid(~exposition) +
  labs(title = "Exposion x Comparison of consecutive surveys")
ggplot(sites, aes(x = pc1_soil, y = y, color = comparison)) +
  geom_point() + geom_smooth(method = "lm") +
  labs(title = "PC1 x Comparison of consecutive surveys")
ggplot(sites, aes(x = pc2_soil, y = y, color = comparison)) +
  geom_point() + geom_smooth(method = "lm") +
  labs(title = "PC2 x Comparison of consecutive surveys")
ggplot(sites, aes(x = pc1_soil, y = y, color = exposition)) +
  geom_point() + geom_smooth(method = "lm") +
  labs(title = "PC1 x Exposition")
ggplot(sites, aes(x = (pc2_soil), y = y, color = exposition)) +
  geom_point() + geom_smooth(method = "lm") +
  labs(title = "PC2 x Exposition")


### b Outliers, zero-inflation, transformations? ------------------------------

sites %>%
  count(location_construction_year)
boxplot(sites$n)
ggplot(sites, aes(x = exposition, y = y)) +
  geom_quasirandom()
ggplot(sites, aes(x = y)) +
  geom_histogram(binwidth = 0.03)
ggplot(sites, aes(x = y)) +
  geom_density()
ggplot(sites, aes(x = log(y))) +
  geom_density()


### c Check collinearity ------------------------------------------------------

data_collinearity <- sites %>%
  select(where(is.numeric), -b, -c, -d, -y) %>%
  GGally::ggpairs(
    data_collinearity,
    lower = list(continuous = "smooth_loess"),
    axisLabels = "internal"
    )
#--> exclude r > 0.7
# Dormann et al. 2013 Ecography
# https://doi.org/10.1111/j.1600-0587.2012.07348.x



## 2 Model building ###########################################################


### a models ------------------------------------------------------------------

### * Random structure ####
m1a <- blmer(y ~ 1 + (1 | location_construction_year),
             data = sites, REML = TRUE)
m1b <- blmer(y ~ 1 + (1 | location_construction_year / plot),
             data = sites, REML = TRUE)
m1c <- blmer(y ~ 1 + (1 | plot), data = sites, REML = TRUE)
MuMIn::AICc(m1a, m1b, m1c) %>%
  arrange(AICc)

### * Fixed effects ####
m1 <- blmer(
  y ~ (comparison + exposition + pc1_soil)^2 + pc2_soil + pc3_soil +
    orientation + river_distance + location_construction_year +
    (1 | plot),
  REML = FALSE,
  control = lmerControl(optimizer = "Nelder_Mead"),
  cov.prior = wishart,
  data = sites
  )
simulateResiduals(m1, plot = TRUE)
m2 <- blmer(
  y ~ comparison + exposition * pc1_soil + pc2_soil + pc3_soil +
    orientation + river_distance + location_construction_year +
    (1 | plot),
  REML = FALSE,
  control = lmerControl(optimizer = "Nelder_Mead"),
  cov.prior = wishart,
  data = sites
)
simulateResiduals(m2, plot = TRUE)
m3 <- blmer(
  y ~ comparison * exposition + pc1_soil + pc2_soil + pc3_soil +
    orientation + river_distance + location_construction_year +
    (1 | plot),
  REML = FALSE,
  control = lmerControl(optimizer = "Nelder_Mead"),
  cov.prior = wishart,
  data = sites
)
simulateResiduals(m3, plot = TRUE)
m4 <- blmer(
  y ~ comparison * pc1_soil + exposition + pc2_soil + pc3_soil +
    orientation + river_distance + location_construction_year +
    (1 | plot),
  REML = FALSE,
  control = lmerControl(optimizer = "Nelder_Mead"),
  cov.prior = wishart,
  data = sites
)
simulateResiduals(m4, plot = TRUE)
m5 <- blmer(
  y ~ comparison + exposition + pc1_soil + pc2_soil + pc3_soil +
    orientation + river_distance + location_construction_year +
    (1 | plot),
  REML = FALSE,
  control = lmerControl(optimizer = "Nelder_Mead"),
  cov.prior = wishart,
  data = sites
)
simulateResiduals(m5, plot = TRUE)
