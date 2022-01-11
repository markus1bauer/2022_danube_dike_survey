# Beta diversity on dike grasslands
# Plot Fig 3C ####
# Markus Bauer
# 2022-01-11
# Citation: 
## Bauer M, Huber J, Kollmann J (submitted) 
## Balanced turnover is a main aspect of biodiversity on restored dike grasslands: not only deterministic environmental effects, but also non-directional year and site effects drive spatial and temporal beta diversity.
## Unpublished data.



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Preparation ################################################################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


### Packages ###
library(here)
library(tidyverse)
library(blme)
library(dotwhisker)

### Start ###
rm(list = setdiff(ls(), c("graph_a", "graph_b", "graph_c", "graph_d")))
setwd(here("data/processed"))


### Load data ###
sites <- read_csv("data_processed_sites_temporal.csv", col_names = T, na = c("", "na", "NA"), col_types = 
                  cols(
                    .default = "?",
                    plot = "f",
                    block = "f",
                    comparison = "f",
                    locationYear = "f",
                    exposition = col_factor(levels = c("south", "north")),
                    side = col_factor(levels = c("land", "water"))
                  )) %>%
  mutate(across(c("longitude", "latitude", "riverkm", "distanceRiver"), scale)) %>%
  mutate(y = C_presence - B_presence)

### * Model ####
m3 <- blmer(y ~ comparison * exposition + PC1soil + PC2soil + PC3soil + 
              side + distanceRiver + locationYear + 
              (1|plot), 
            REML = F,
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



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Plot ##############################################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


(graph_c <- m3 %>% 
   broom.mixed::tidy(conf.int = T, conf.level = .95) %>%
   filter(
     !str_detect(term, "location*") & 
       !str_detect(term, "comparison") & 
       !str_detect(term, "sd_*") &
       !str_detect(term, "(Intercept)")) %>%
   mutate(term = fct_relevel(term, c("PC3soil", "PC2soil", "PC1soil", 
                                     "distanceRiver", "sidewater", "expositionnorth")),
          term = fct_recode(term,
                            "South | North exposition" = "expositionnorth",
                            "Water | Land side" = "sidewater",
                            "Distance to river" = "distanceRiver")) %>%
   ggplot(aes(x = estimate, y = term, xmin = conf.low, xmax = conf.high)) +
   geom_vline(xintercept = 0, linetype = 2, color = "black") +
   geom_point(size = 2, shape = "circle open") +
   geom_linerange() +
   labs(x = expression(Estimate~"["*italic('C')[sor]-italic('B')[sor]*"]")) +
   themeMB())

### Save ###
ggsave(here("outputs/figures/figure_3c_(800dpi_8x8cm).tiff"),
       dpi = 800, width = 8, height = 8, units = "cm")
