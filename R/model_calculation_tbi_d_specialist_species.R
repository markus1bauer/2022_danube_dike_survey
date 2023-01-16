# Beta diversity on dike grasslands
# Temporal beta-diversity index (TBI) - Specialist species ####
# Markus Bauer
# 2023-01-16



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Preparation ###############################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



### Packages ###
library(here)
library(tidyverse)
library(ggbeeswarm)
library(patchwork)
library(blme)
library(DHARMa)
library(emmeans)

### Start ###
rm(list = ls())
setwd(here("data", "processed"))

### Load data ###
sites <- read_csv(here("data", "processed", "data_processed_sites_temporal.csv"),
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
      pool == "target" & presabu == "presence"
  ) %>%
  mutate(
    y = d,
    comparison = factor(comparison)
  ) %>%
  mutate(
    river_km_scaled = scale(river_km),
    river_distance_scaled = scale(river_distance),
    biotope_distance_scaled = scale(biotope_distance),
    biotope_area_scaled = scale(biotope_area)
  )



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Statistics ################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



## 1 Data exploration #########################################################


### a Means and deviation -----------------------------------------------------

mean(sites$y)
median(sites$y)
sd(sites$y)
Rmisc::CI(sites$y, ci = .95)
quantile(sites$y, probs = c(0.05, 0.95), na.rm = TRUE)


### b Graphs ------------------------------------------------------------------

plot1 <- ggplot(sites, aes(x = comparison, y = y)) +
  geom_quasirandom(color = "grey") + geom_boxplot(fill = "transparent") +
  labs(title = "Comparison of consecutive surveys")
plot2 <- ggplot(sites, aes(x = exposition, y = y)) +
  geom_quasirandom(color = "grey") + geom_boxplot(fill = "transparent") +
  labs(title = "Exposition of dike slopes")
plot3 <- ggplot(sites, aes(x = orientation, y = y)) +
  geom_quasirandom(color = "grey") + geom_boxplot(fill = "transparent") +
  labs(title = "Orientation of dike slopes")
plot4 <- ggplot(sites, aes(x = river_km, y = (y))) +
  geom_point() +  geom_smooth(method = "lm") +
  labs(title = "Position along the river")
(plot1 + plot2) / (plot3 + plot4)
plot1 <- ggplot(sites, aes(x = location_construction_year, y = y)) +
  geom_quasirandom(color = "grey") + geom_boxplot(fill = "transparent") +
  labs(title = "Location and construction year of the dike")
plot2 <- ggplot(sites, aes(x = (river_distance), y = (y))) +
  geom_point() + geom_smooth(method = "lm") +
  labs(title = "Distance to river course")
plot3 <- ggplot(sites, aes(x = (biotope_distance), y = (y))) +
  geom_point() + geom_smooth(method = "lm") +
  labs(title = "Distance to closest grassland biotope")
plot4 <- ggplot(sites, aes(x = (biotope_area), y = (y))) +
  geom_point() + geom_smooth(method = "lm") +
  labs(title = "Amount of grassland biotopes with 500 m radius")
(plot1 + plot2) / (plot3 + plot4)
plot1 <- ggplot(sites, aes(x = pc1_soil, y = (y))) +
  geom_point() + geom_smooth(method = "lm") +
  labs(title = "PC1 (soil)")
plot2 <- ggplot(sites, aes(x = (pc2_soil), y = y)) +
  geom_point() + geom_smooth(method = "lm") +
  labs(title = "PC2 (soil)")
plot3 <- ggplot(sites, aes(x = (pc3_soil), y = y)) +
  geom_point() + geom_smooth(method = "lm") +
  labs(title = "PC3 (soil)")
plot4 <- ggplot(sites, aes(x = comparison, y = y)) +
  geom_quasirandom(color = "grey") + geom_boxplot(fill = "transparent") +
  facet_grid(~exposition) +
  labs(title = "Exposion x Comparison of consecutive surveys")
(plot1 + plot2) / (plot3 + plot4)
plot1 <- ggplot(sites, aes(x = pc1_soil, y = y, color = comparison)) +
  geom_point() + geom_smooth(method = "lm") +
  labs(title = "PC1 x Comparison of consecutive surveys")
plot2 <- ggplot(sites, aes(x = pc2_soil, y = y, color = comparison)) +
  geom_point() + geom_smooth(method = "lm") +
  labs(title = "PC2 x Comparison of consecutive surveys")
plot3 <- ggplot(sites, aes(x = pc1_soil, y = y, color = exposition)) +
  geom_point() + geom_smooth(method = "lm") +
  labs(title = "PC1 x Exposition")
plot4 <- ggplot(sites, aes(x = (pc2_soil), y = y, color = exposition)) +
  geom_point() + geom_smooth(method = "lm") +
  labs(title = "PC2 x Exposition")
(plot1 + plot2) / (plot3 + plot4)


### c Outliers, zero-inflation, transformations? ------------------------------

sites %>%
  count(location_construction_year)
plot1 <- ggplot(sites, aes(x = exposition, y = y)) +
  geom_quasirandom()
plot2 <- ggplot(sites, aes(x = y)) +
  geom_histogram(binwidth = 0.03)
plot3 <- ggplot(sites, aes(x = y)) +
  geom_density()
plot4 <- ggplot(sites, aes(x = log(y))) +
  geom_density()
(plot1 + plot2) / (plot3 + plot4)


### d Check collinearity ------------------------------------------------------

sites %>%
  select(where(is.numeric), -b, -c, -d, -y, -ends_with("scaled")) %>%
  GGally::ggpairs(
    lower = list(continuous = "smooth_loess")
  ) +
  theme(strip.text = element_text(size = 7))
sites <- sites %>%
  select(-biotope_area)
#--> exclude r > 0.7
# Dormann et al. 2013 Ecography
# https://doi.org/10.1111/j.1600-0587.2012.07348.x



## 2 Model building ###########################################################


### a Random structure ---------------------------------------------------------

m1a <- blmer(y ~ 1 + (1 | location_construction_year),
             data = sites, REML = TRUE)
m1b <- blmer(y ~ 1 + (1 | location_construction_year / plot),
             data = sites, REML = TRUE)
m1c <- blmer(y ~ 1 + (1 | plot), data = sites, REML = TRUE)
MuMIn::AICc(m1a, m1b, m1c) %>%
  arrange(AICc)


### b Fixed effects ------------------------------------------------------------

m1 <- blmer(
  log(y) ~ (comparison + exposition + pc1_soil)^2 + pc2_soil + pc3_soil +
    orientation + river_distance_scaled + river_km_scaled +
    biotope_distance_scaled +
    (1 | plot),
  REML = FALSE,
  control = lmerControl(optimizer = "Nelder_Mead"),
  cov.prior = wishart,
  data = sites
)
simulateResiduals(m1, plot = TRUE)
m2 <- blmer(
  log(y) ~ comparison * exposition + pc1_soil + pc2_soil + pc3_soil +
    orientation + river_distance_scaled + river_km_scaled +
    biotope_distance_scaled +
    (1 | plot),
  REML = FALSE,
  control = lmerControl(optimizer = "Nelder_Mead"),
  cov.prior = wishart,
  data = sites
)
simulateResiduals(m2, plot = TRUE)
m3 <- blmer(
  log(y) ~ comparison + exposition * pc1_soil + pc2_soil + pc3_soil +
    orientation + river_distance_scaled + river_km_scaled +
    biotope_distance_scaled +
    (1 | plot),
  REML = FALSE,
  control = lmerControl(optimizer = "Nelder_Mead"),
  cov.prior = wishart,
  data = sites
)
simulateResiduals(m3, plot = TRUE)
m4 <- blmer(
  log(y) ~ comparison * pc1_soil + exposition + pc2_soil + pc3_soil +
    orientation + river_distance_scaled + river_km_scaled +
    biotope_distance_scaled +
    (1 | plot),
  REML = FALSE,
  control = lmerControl(optimizer = "Nelder_Mead"),
  cov.prior = wishart,
  data = sites
)
simulateResiduals(m4, plot = TRUE)
m5 <- blmer(
  log(y) ~ comparison + exposition + pc1_soil + pc2_soil + pc3_soil +
    orientation + river_distance_scaled + river_km_scaled +
    biotope_distance_scaled +
    (1 | plot),
  REML = FALSE,
  control = lmerControl(optimizer = "Nelder_Mead"),
  cov.prior = wishart,
  data = sites
)
simulateResiduals(m5, plot = TRUE)


### c Save ---------------------------------------------------------------------

save(m1, file = here("outputs", "models", "model_tbi_d_specialist_1.Rdata"))
save(m3, file = here("outputs", "models", "model_tbi_d_specialist_3.Rdata"))
save(m5, file = here("outputs", "models", "model_tbi_d_specialist_5.Rdata"))
