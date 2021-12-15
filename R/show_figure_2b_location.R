# Show figure 2 ####
# Markus Bauer
# Citation: Markus Bauer 



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Preparation ################################################################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


### Packages ###
library(here)
library(tidyverse)
library(blme)
library(ggeffects)
library(ggbeeswarm)

### Start ###
rm(list = setdiff(ls(), c("graph_a", "graph_b", "graph_c", "graph_d")))
setwd(here("data/processed"))


### Load data ###
tbi <- read_csv("data_processed_tbi.csv", col_names = T, na = c("", "na", "NA"), col_types = 
                  cols(
                    .default = "?"
                  )) %>%
  filter(comparison %in% c("1718", "1819", "1921")) %>%
  pivot_wider(id_cols = c(plot, comparison, exposition, side, block, 
                          location, locationYear, constructionYear, 
                          longitude, latitude, riverkm, distanceRiver, 
                          PC1soil, PC2soil, PC3soil, conf.low, conf.high), 
              names_from = "presabu", 
              values_from = "D") %>%
  mutate(plot = factor(plot),
         block = factor(block),
         comparison = factor(comparison),
         exposition = factor(exposition),
         side = factor(side),
         constructionYear = factor(constructionYear),
         locationYear = factor(locationYear)) %>%
  rename(y = presence) %>%
  mutate(across(where(is.numeric) & !y, scale))

### * Model ####
m2 <- blmer(log(y) ~ comparison + exposition * PC1soil + PC2soil + PC3soil + 
              side + distanceRiver + locationYear + abundance +
            (1|plot), 
            REML = T,
            control = lmerControl(optimizer = "Nelder_Mead"),
            cov.prior = wishart,
            data = tbi)

### * Functions ####
themeMB <- function(){
  theme(
    panel.background = element_rect(fill = "white"),
    text  = element_text(size = 9, color = "black"),
    strip.text = element_text(size = 10),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0, size = 9, color = "black"),
    axis.line.x = element_line(),
    axis.line.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.key = element_rect(fill = "white"),
    legend.position = "none",
    legend.margin = margin(0, 0, 0, 0, "cm"),
    plot.margin = margin(0, 0, 0, 0, "cm")
  )
}



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Plot ##############################################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

data_model <- ggeffect(m2, type = "emm", c("locationYear"), back.transform = T) %>%
  mutate(predicted = exp(predicted),
         conf.low = exp(conf.low),
         conf.high = exp(conf.high),
         cross = if_else(x %in% c("PFE-2008", "IRL-2003"), "filled", "open"))

data <- tbi %>%
  rename(predicted = y, x = locationYear) 

(graph_b <- ggplot() +
    geom_quasirandom(data = data, 
                     aes(x = x, predicted),
                     dodge.width = .6, size = 1, shape = 16, color = "grey70") + 
    geom_hline(yintercept = c(mean(tbi$y), mean(tbi$y) + 0.5 * sd(tbi$y), mean(tbi$y) - 0.5 * sd(tbi$y)), 
               linetype = c(1, 2, 2),
               color = "grey70") +
    geom_errorbar(data = data_model, 
                  aes(x, predicted, ymin = conf.low, ymax = conf.high), 
                  width = 0.0, size = 0.4) +
    geom_point(data = data_model,
               aes(x, predicted, shape = cross),
               size = 2) +
    scale_y_continuous(limits = c(0, .92), breaks = seq(0, 400, .1)) +
    scale_shape_manual(values = c("circle", "circle open")) +
    labs(x = "", y = expression(Temporal~"beta"~diversity~"["*italic('D')[sor]*"]")) +
    themeMB())

### Save ###
ggsave(here("outputs/figures/figure_2b_(800dpi_8x8cm).tiff"),
       dpi = 800, width = 8, height = 8, units = "cm")
