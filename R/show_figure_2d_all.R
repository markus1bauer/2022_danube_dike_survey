# Show figure 2 ####
# Markus Bauer
# Citation: Markus Bauer 



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Preparation ################################################################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


### Packages ###
library(here)
library(tidyverse)
library(lme4)
library(dotwhisker)

### Start ###
rm(list = setdiff(ls(), c("graph_a", "graph_b", "graph_c", "graph_d")))
setwd(here("data/processed"))


### Load data ###
tbi <- read_csv("data_processed_tbi.csv", col_names = T, na = c("", "na", "NA"), col_types = 
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
  filter(comparison %in% c("1718", "1819", "1921") & presabu == "presence") %>%
  mutate(comparison = factor(comparison),
         locationYear = factor(locationYear),
         block = factor(block),
         plot = factor(plot),
         exposition = factor(exposition)
  ) %>%
  rename(y = D) %>%
  mutate(across(where(is.numeric) & !y, scale))

### * Model ####
m5 <- lmer(log(y) ~ comparison + exposition + side + PC1soil + PC2soil + PC3soil + distanceRiver + locationYear + 
             exposition:PC2soil +
             (1|plot), 
           REML = T,
           data = tbi)

### * Functions ####
themeMB <- function(){
  theme(
    panel.background = element_rect(fill = "white"),
    text  = element_text(size = 9, color = "black"),
    strip.text = element_text(size = 10),
    axis.text.y = element_text(angle = 0, hjust = 0.5, size = 9, color = "black"),
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



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Plot ##############################################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


(graph_d <- m5 %>% 
   broom.mixed::tidy(conf.int = T, conf.level = .95) %>%
   filter(
     !str_detect(term, "location*") & 
      !str_detect(term, "comparison*") & 
      !str_detect(term, "sd_*") &
      !str_detect(term, "expositionnorth:PC2soil") &
      !str_detect(term, "(Intercept)")) %>%
   mutate(cross = if_else(term %in% c("sideland"), "filled", "open"),
          term = fct_relevel(term, c("PC3soil", "PC2soil", "PC1soil", "distanceRiver", "sideland", "expositionnorth")),
          term = fct_recode(term,
                            "South | North exposition" = "expositionnorth",
                            "Water | Land side" = "sideland",
                            "Distance to river" = "distanceRiver")) %>%
   ggplot(aes(x = estimate, y = term, xmin = conf.low, xmax = conf.high)) +
   geom_vline(xintercept = 0, linetype = 2, color = "black") +
   geom_point(aes(shape = cross), size = 2) +
   geom_linerange() +
   scale_shape_manual(values = c("circle", "circle open")) +
   labs(x = expression("log("*italic('D')[sor]*")")) +
   themeMB())

### Save ###
ggsave(here("outputs/figures/figure_2d_(800dpi_8x8cm).tiff"),
       dpi = 800, width = 8, height = 8, units = "cm")
