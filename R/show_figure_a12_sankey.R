# Beta diversity on dike grasslands
# Figure S12 ####
# Sankey Diagram of habitat types

# Markus Bauer
# 2024-02-28



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Preparation ################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



### Packages ###
library(here)
library(tidyverse)
library(remotes)
#remotes::install_github("davidsjoberg/ggsankey")
library(ggsankey)

### Start ###
rm(list = setdiff(ls(), c("graph_a", "graph_b", "graph_c", "graph_d")))

### Functions ###
theme_mb <- function() {
  theme(
    panel.background = element_rect(fill = "white"),
    text = element_text(size = 9, color = "black"),
    strip.text = element_text(size = 10),
    axis.text.x = element_text(
      angle = 0, hjust = 0.5, size = 9, color = "black"
      ),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.line = element_blank(),
    legend.position = "none",
    plot.margin = margin(0, 0, 0, 0, "cm")
  )
}

### Load data sites ###
sites <- read_csv(
  here("data", "processed", "data_processed_sites_spatial.csv"),
  col_names = TRUE, na = c("na", "NA", ""), col_types = cols(
    .default = "?",
    id = "f",
    survey_year = "f",
    esy = "f"
  )
) %>%
  select(plot, survey_year,esy) %>%
  pivot_wider(id_cols = plot, names_from = "survey_year", values_from = "esy") %>%
  make_long("2017", "2018", "2019", "2021")



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Plot ######################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



(graph_b <- ggplot(
  aes(
    x = x, next_x = next_x, node = node, next_node = next_node, label = node,
    fill = factor(node)
  ),
  data = sites
) +
  geom_sankey(flow.alpha = 0.5, node.color = 1) +
  geom_sankey_label(size = 3.5, color = 1, fill = "white") +
  scale_fill_manual(
    values = c(
      "R22" = "#00BFC4",
      "R21" = "#00BFC4",
      "R1A" = "#F8766D",
      "R" = "grey30",
      "V38" = "#C77CFF",
      "?" = "grey90"
    )
  ) +
  theme_mb())



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# C Save ######################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



ggsave(
  here("outputs", "figures", "figure_3_800dpi_16.5x10cm.tiff"),
  dpi = 800, width = 16.5, height = 10, units = "cm"
)
