# Model for ordination ####
# Markus Bauer



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Preparation ################################################################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


### Packages ###
library(here)
library(tidyverse)
library(vegan)

### Start ###
rm(list = ls())
setwd(here("data/processed"))

### Load data ###
sites <- read_csv2("data_processed_sites.csv", col_names = T, na = "na", col_types = 
                     cols(
                       .default = col_guess(),
                       id = col_factor(),
                       block = col_factor(),
                       plot = col_factor(),
                       exposition = col_factor(),
                       ffh = col_factor(),
                       changeType = col_factor(c("better", "change", "worse", "any-FFH", "FFH6510", "non-FFH"))
                     )) %>%
  select(permdisp, id, block, plot, exposition, plotAge, phosphorous, sandPerc, topsoilDepth, changeType, ffh) %>%
  filter(ffh != "6210") %>%
  filter(changeType != "any-FFH") %>%
  filter(!is.na(changeType))


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Statistics ################################################################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#### a PERMDISP ----------------------------------------------------------------------------------------
ggplot(aes(x = ffh, y = permdisp), data = sites) +
  geom_boxplot()
ggplot(aes(x = changeType, y = permdisp), data = sites) +
  geom_boxplot()
ggplot(aes(x = exposition, y = permdisp), data = sites) +
  geom_boxplot()
ggplot(aes(x = sandPerc, y = permdisp), data = sites) +
  geom_point() +
  geom_smooth(method = "lm")
ggplot(aes(x = phosphorous, y = permdisp), data = sites) +
  geom_point() +
  geom_smooth(method = "lm")
ggplot(aes(x = topsoilDepth, y = permdisp), data = sites) +
  geom_point() +
  geom_smooth(method = "lm")
ggplot(aes(x = plotAge, y = permdisp), data = sites) +
  geom_smooth(aes(group = plot), colour = "grey40", method = "lm", se = F, show.legend = F) +
  geom_point() +
  geom_smooth(method = "lm")
ggplot(aes(x = plotAge, y = permdisp, color = ffh), data = sites) +
  geom_smooth(aes(group = plot), colour = "grey40", method = "lm", se = F, show.legend = F) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_grid(~ffh)
ggplot(aes(x = plotAge, y = permdisp, color = changeType), data = subset(sites, !(changeType %in% c("NA", "any-FFH")))) +
  geom_smooth(aes(group = plot), colour = "grey40", method = "lm", se = F, show.legend = F) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_grid(~changeType)
ggplot(aes(x = plotAge, y = permdisp, color = exposition), data = sites) +
  geom_smooth(aes(group = plot), colour = "grey40", method = "lm", se = F, show.legend = F) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_grid(~exposition)


