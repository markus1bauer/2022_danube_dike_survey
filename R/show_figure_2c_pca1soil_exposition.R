# Beta diversity on dike grasslands
# Plot Fig 2C ####
# Markus Bauer
# 2022-01-11



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Preparation #########################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


### Packages ###
library(here)
library(tidyverse)
library(blme)
library(ggeffects)

### Start ###
rm(list = setdiff(ls(), c("graph_a", "graph_b", "graph_c", "graph_d")))
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
      exposition = col_factor(levels = c("north", "south")),
      side = "f",
      locationYear = "f"
    )
) %>%
  filter(comparison == "1718" | comparison == "1819" | comparison == "1921") %>%
  mutate(
    y = D_presence,
    comparison = factor(comparison)
  ) %>%
  mutate(across(where(is.numeric) & !y, scale))

### * Model ####
m2 <- blmer(log(y) ~ comparison + exposition * PC1soil + PC2soil + PC3soil +
  side + distanceRiver + locationYear + D_abundance +
  (1 | plot),
REML = T,
control = lmerControl(optimizer = "Nelder_Mead"),
cov.prior = wishart,
data = sites
)

### * Functions ####
theme_mb <- function() {
  theme(
    panel.background = element_rect(fill = "white"),
    text = element_text(size = 9, color = "black"),
    strip.text = element_text(size = 10),
    axis.text.y = element_text(angle = 0, hjust = 0.5, size = 9, color = "black"),
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 9, color = "black"),
    axis.title = element_text(angle = 0, hjust = 0.5, size = 9, color = "black"),
    axis.line = element_line(),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    legend.key = element_rect(fill = "white"),
    legend.margin = margin(0, 0, 0, 0, "cm"),
    plot.margin = margin(0, 0, 0, 0, "cm")
  )
}



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Plot ################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

data_model <- ggeffect(m2, type = "emm", c("PC1soil", "exposition"), back.transform = TRUE) %>%
  mutate(
    predicted = exp(predicted),
    conf.low = exp(conf.low),
    conf.high = exp(conf.high),
    group = fct_recode(group, "North" = "north", "South" = "south")
  )

data <- sites %>%
  rename(predicted = y, x = PC1soil, group = exposition) %>%
  mutate(group = fct_recode(group, "North" = "north", "South" = "south"))

(graph_c <- ggplot() +
  geom_ribbon(
    data = data_model,
    aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, fill = group),
    alpha = .3
  ) +
  geom_point(
    data = data,
    aes(x = x, y = predicted, shape = group),
    size = 1, color = "grey70", fill = "grey70"
  ) +
  geom_line(
    data = data_model,
    aes(x = x, y = predicted, group = group)
  ) +
  scale_y_continuous(limits = c(0, .92), breaks = seq(-100, 400, .1)) +
  scale_x_continuous(breaks = seq(-100, 400, 1)) +
  scale_fill_manual(values = c("grey40", "grey70")) +
  scale_shape_manual(values = c("circle filled", "circle open")) +
  annotate("text",
    label = c("Nitrogen", "Sand"),
    x = c(-1.8, 2.8),
    y = c(0, 0),
    size = 2.5
  ) +
  labs(x = expression(PC1[soil]), y = expression(Temporal ~"beta"~ diversity ~ "[" * italic("D")[sor] * "]"), color = "Exposition", fill = "Exposition", shape = "Exposition") +
  theme_mb() +
  theme(legend.position = c(.8, .9)))

### Save ###
ggsave(here("outputs", "figures", "figure_2c_800dpi_8x8cm.tiff"),
  dpi = 800, width = 8, height = 8, units = "cm"
)
