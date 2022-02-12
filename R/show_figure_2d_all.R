# Beta diversity on dike grasslands
# Plot Fig 2D ####
# Markus Bauer
# 2022-01-11



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Preparation ##########################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


### Packages ###
library(here)
library(tidyverse)
library(blme)

### Start ###
rm(list = setdiff(ls(), c("graph_a", "graph_b", "graph_c", "graph_d")))
setwd(here("data", "processed"))


### Load data ###
sites <- read_csv("data_processed_sites_temporal.csv", col_names = TRUE,
                  na = c("", "na", "NA"), col_types =
                  cols(
                    .default = "?",
                    side = col_factor(levels = c("land", "water")),
                    exposition = col_factor(levels = c("south", "north")),
                    plot = "f",
                    block = "f",
                    comparison = "f",
                    locationYear = "f"
                  )) %>%
  rename(y = D_presence) %>%
  mutate(across(where(is.numeric) & !y, scale))

### * Model ####
m2 <- blmer(log(y) ~ comparison + exposition * PC1soil + PC2soil + PC3soil + 
              side + distanceRiver + locationYear + abundance +
              (1|plot), 
            REML = T,
            control = lmerControl(optimizer = "Nelder_Mead"),
            cov.prior = wishart,
            data = sites)

### * Functions ####
themeMB <- function(){
  theme(
    panel.background = element_rect(fill = "white"),
    text  = element_text(size = 9, color = "black"),
    strip.text = element_text(size = 10),
    axis.text.y = element_text(angle = 0, hjust = 1, size = 9, color = "black"),
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 9, color = "black"),
    axis.title.x = element_text(angle = 0, hjust = 0.5, size = 9, color = "black"),
    axis.title.y = element_blank(),
    axis.line = element_line(),
    legend.key = element_rect(fill = "white"),
    legend.position = "none",
    legend.margin = margin(0, 0, 0, 0, "cm"),
    plot.margin = margin(0, 0, 0, 0, "cm")
  )
}



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Plot #################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


(graph_d <- m2 %>% 
   broom.mixed::tidy(conf.int = T, conf.level = .95) %>%
   filter(
     !str_detect(term, "location*") & 
      !str_detect(term, "comparison*") & 
      !str_detect(term, "sd_*") &
      !str_detect(term, "(Intercept)")) %>%
   mutate(cross = if_else(term %in% c("sidewater", "abundance"), "filled", "open"),
          term = fct_relevel(term, c("expositionnorth:PC1soil", "PC3soil", "PC2soil", "PC1soil", 
                                     "distanceRiver", "sidewater", "expositionnorth", "abundance")),
          term = fct_recode(term,
                            "South | North exposition" = "expositionnorth",
                            "Land | Water side" = "sidewater",
                            "Distance to river" = "distanceRiver",
                            "D [bc]" = "abundance",
                            "Exposition:PC1soil" = "expositionnorth:PC1soil")) %>%
   ggplot(aes(x = estimate, y = term, xmin = conf.low, xmax = conf.high)) +
   geom_vline(xintercept = 0, linetype = 2, color = "black") +
   geom_point(aes(shape = cross), size = 2) +
   geom_linerange() +
   scale_shape_manual(values = c("circle", "circle open")) +
   labs(x = expression("Estimate [log("*italic('D')[sor]*")]")) +
   themeMB())

### Save ###
ggsave(here("outputs", "figures", "figure_2d_800dpi_8x8cm.tiff"),
       dpi = 800, width = 8, height = 8, units = "cm")
