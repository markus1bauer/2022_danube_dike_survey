# Calculate TBI ####
# Markus Bauer
# Citation: Bauer M, Huber J, Kollmann J (202x)



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Preparation ################################################################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


### Packages ###
library(here)
library(tidyverse)
library(lme4)  
library(DHARMa)

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
  filter(comparison %in% c("1718", "1819", "1921")) %>%
  mutate(comparison = factor(comparison)) %>%
  pivot_wider(names_from = "presabu", values_from = "D") %>%
  group_by(plot, comparison, exposition, side, locationYear) %>%
  summarise(across(c(presence, abundance, PC1soil, PC2soil, PC3soil, distanceRiver), ~ max(.x, na.rm = T))) %>%
  ungroup()



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Statistics ################################################################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


## 1 Data exploration #####################################################################################

### a Graphs ---------------------------------------------------------------------------------------------
ggplot(tbi, aes(x = presence, y = abundance)) + 
  geom_point() +
  geom_smooth(method = "lm")
ggplot(data = tbi, aes(x = presence, y = abundance)) +
  geom_point() +
  geom_density_2d_filled(alpha = .5,
                         contour_var = "ndensity",
                         show.legend = F) +
  coord_fixed(ratio = 1, xlim = c(0, 1), ylim = c(0, 1))
  
m <- lmer((abundance) ~ presence + (PC2soil) + comparison + exposition + PC1soil + PC3soil + side + log(distanceRiver) + locationYear +
            presence:PC2soil + presence:PC1soil + 
            (1|plot), 
          REML = T,
          data = tbi);simulateResiduals(m, plot = T)
car::Anova(m, type = 2)
visreg::visreg2d(m, "presence", "PC2soil", trans = log)

