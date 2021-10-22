# Show Figure targetRichness ####
# Markus Bauer



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Preparation ################################################################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


### Packages ###
library(here)
library(tidyverse)

### Start ###
rm(list = ls())

### Load data ###
x = seq(1, 1000, by = 1)

### Functions ###
themeMB <- function(){
  theme(
    panel.background = element_rect(fill = "white"),
    text  = element_text(size = 8, color = "black"),
    strip.text = element_text(size = 10),
    axis.text = element_blank(),
    axis.line = element_line(),
    axis.ticks = element_blank(),
    legend.key = element_rect(fill = "white"),
    legend.position = "right",
    legend.margin = margin(0, 0, 0, 0, "cm"),
    plot.margin = margin(0, 0, 0, 0, "cm")
  )
}



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Plot ################################################################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


# 1 Classic concept ################################################################################################################

set.seed(5)
a1 <- .85
b1 <- .18
y1 = a1 + 
  b1 * log(x) * (1 + rnorm(length(x), sd = 0.001)) + 
  rnorm(length(x), sd = 1.5)
a2 <- 1.2
b2 <- .11
y2 = a2 + 
  b2 * log(x) * (1 + rnorm(length(x), sd = 0.001)) + 
  rnorm(length(x), sd = 1.5)
a3 <- 0.5
b3 <- .08
y3 = a3 + 
  b3 * log(x) * (1 + rnorm(length(x), sd = 0.001)) + 
  rnorm(length(x), sd = 2)
data <- data_frame(x, y1, y2, y3) %>%
  pivot_longer(-x, names_to = "treatment", values_to = "y")
data <- data %>%
  mutate(treatment = fct_recode(treatment, "Treatment 1" = "y1", "Treatment 2" = "y2", "Treatment 3" = "y3"))
ggplot(data, aes(x = x, y = y)) +
  annotate("rect", xmin = 0, xmax = 1000, ymin = 1.85, ymax = 2.4, 
           fill = "white") +
  geom_smooth(aes(color = treatment), method = "lm", se = T, 
              formula = y ~ log(x)) +
  annotate("segment", x = 0, xend = 1000, y = 1.85, yend = 1.85, 
           alpha = 1, linetype = 2) +
  annotate("segment", x = 0, xend = 1000, y = 2.4, yend = 2.4, 
           alpha = 1, linetype = 1) +
  annotate("text", x = 180, y = 2.3, size = 2.5, 
           label = "Defined area of reference") +
  coord_cartesian(ylim = c(0, 2.409), expand = F) +
  scale_colour_manual(values = c("forestgreen", "gold3", "firebrick")) +
  labs(x = "Time since restoration", 
       y = expression(paste(Delta, " Species composition")),
       colour = "") +
  themeMB()

### Save ###
setwd(here("outputs/figures"))
ggsave("figure_concept_classic_(300dpi_12x8cm).tiff",
       dpi = 300, width = 12, height = 8, units = "cm")


# 2 Variation concept ################################################################################################################

set.seed(5)
a1 <- .8
b1 <- .18
y1 = a1 + 
  b1 * log(x) * (1 + rnorm(length(x), sd = 0.001)) + 
  rnorm(length(x), sd = 2)
a2 <- 1.2
b2 <- .11
y2 = a2 + 
  b2 * log(x) * (1 + rnorm(length(x), sd = 0.001)) + 
  rnorm(length(x), sd = 7)
a3 <- 0.5
b3 <- .08
y3 = a3 + 
  b3 * log(x) * (1 + rnorm(length(x), sd = 0.001)) + 
  rnorm(length(x), sd = 8)
data <- data_frame(x, y1, y2, y3) %>%
  pivot_longer(-x, names_to = "treatment", values_to = "y")
data <- data %>%
  mutate(treatment = fct_recode(treatment, "Treatment 1" = "y1", "Treatment 2" = "y2", "Treatment 3" = "y3"))
ggplot(data, aes(x = x, y = y)) +
    annotate("rect", xmin = 0, xmax = 1000, ymin = 1.85, ymax = 2.4, 
             fill = "white") +
    geom_smooth(aes(color = treatment), method = "lm", se = T, 
              formula = y ~ log(x)) +
    annotate("segment", x = 0, xend = 1000, y = 1.85, yend = 1.85, 
             alpha = 1, linetype = 2) +
    annotate("segment", x = 0, xend = 1000, y = 2.4, yend = 2.4, 
            alpha = 1, linetype = 1) +
    annotate("text", x = 180, y = 2.3, size = 2.5, 
            label = "Defined area of reference") +
    coord_cartesian(ylim = c(0, 2.409), expand = F) +
    scale_colour_manual(values = c("forestgreen", "gold3", "firebrick")) +
    labs(x = "Time since restoration", 
         y = expression(paste(Delta, " Species composition")),
         colour = "") +
    themeMB()

### Save ###
setwd(here("outputs/figures"))
ggsave("figure_concept_variation_(300dpi_12x8cm).tiff",
       dpi = 300, width = 12, height = 8, units = "cm")
