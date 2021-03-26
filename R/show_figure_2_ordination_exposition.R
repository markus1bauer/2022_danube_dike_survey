# Show Figure ordination 1984-1993-2018 ####
# Markus Bauer



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Preparation ################################################################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


### Packages ###
library(here)
library(tidyverse)
library(vegan)
library(ggrepel)

### Start ###
rm(list = ls())
setwd(here("data/processed"))

### Load data ###
sites <- read_csv2("data_processed_sites.csv", col_names = T, na = "na", col_types = 
                     cols(
                       .default = col_double(),
                       id = col_factor(),
                       location = col_factor(),
                       block = col_factor(),
                       plot = col_factor(),
                       position = col_skip(),
                       side = col_factor(),
                       exposition = col_factor(c("north", "east", "south", "west")),
                       phosphorousClass = col_factor(),
                       biotopeType = col_factor(),
                       baykompv = col_factor(),
                       ffh = col_factor(levels = c("6210", "6510", "non-FFH")),
                       min8 = col_skip(),
                       min9 = col_skip()
                     )        
)

species <- read_csv2("data_processed_species.csv", col_names = T, na = "na", col_types = 
                       cols(
                         .default = col_double(),
                         name = col_factor()
                       )        
)

traits <- read_csv2("data_processed_traits.csv", col_names = T, na = "na", col_types = 
                      cols(
                        .default = col_double(),
                        name = col_factor(),
                        abb = col_factor(),
                        family = col_factor(),
                        rlg = col_factor(),
                        rlb = col_factor(),
                        targetHerb = col_factor(),
                        targetGrass = col_factor(),
                        ffh6510 = col_factor(),
                        ffh6210 = col_factor(),
                        ruderalIndicator = col_factor(),
                        leanIndicator = col_factor(),
                        target = col_factor()
                      )) %>% 
  right_join(species, by = "name") %>%
  select(name, family, abb) %>%
  mutate(family = if_else(family == "Poaceae" | family == "Cyperaceae" | family == "Juncaeae", 
                          "Graminoids", if_else(family == "Fabaceae", 
                                                "Legumes", "Forbs")))

species <- species %>%  
  pivot_longer(-name, "site", "value") %>%
  pivot_wider(site, name) %>%
  column_to_rownames("site")

#### a Choosen model ----------------------------------------------------------------------------------------
(ordi <- metaMDS(species, try = 99, previous.best = T, na.rm = T))
(ef <- envfit(ordi ~  vegetationCov + plotAge + topsoilDepth + phosphorous + pH + NtotalPerc + sandPerc + siltPerc, 
              data = sites, permu = 999, na.rm = T)) #model: ef_vector2



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Plot ################################################################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


themeMB <- function(){
  theme(
    panel.background = element_rect(fill = "white"),
    text  = element_text(size = 10, color = "black"),
    axis.line.y = element_line(),
    axis.line.x = element_line(),
    axis.ticks.x = element_blank(),
    legend.key = element_rect(fill = "white"),
    legend.position = "bottom",
    legend.margin = margin(-.4, 0, 0, 0, "cm"),
    plot.margin = margin(.1, .15, 0, 0, "cm")
  )
}

data.scores <- as.data.frame(scores(ordi)) #input of model
data.scores$site <- rownames(data.scores)
data.scores$variable <- sites$exposition #Write data and 1. variable
### Create Ellipses ###
veganCovEllipse <- function(cov, center = c(0,0), scale = 1, npoints = 100)
{
  theta <- (0:npoints) * 2 * pi / npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

df_ell <- data.frame()
for(g in levels(data.scores$variable)){
  df_ell <- rbind(df_ell, cbind(as.data.frame(with(data.scores[data.scores$variable == g,],
                                                   veganCovEllipse(cov.wt(cbind(NMDS1,NMDS2), wt = rep(1 / length(NMDS1), length(NMDS1)))$cov, center = c(mean(NMDS1), mean(NMDS2)))))
                                ,variable = g))
}

data.ef <- as.data.frame(ef$vectors$arrows * ((sqrt(ef$vectors$r)))) * 2
data.ef$variables <- rownames(data.ef)
set.seed(4393)
### Plot ###
ggplot(data = data.scores) +
  geom_path(aes(x = NMDS1, y = NMDS2, color = variable), data = df_ell,
            size = 1, linetype = 1) +
  geom_point(aes(x = NMDS1, y = NMDS2, shape = variable, color = variable), data = data.scores, 
             size = 3) +
  geom_segment(aes(x = 0, xend = (NMDS1), y = 0, yend = (NMDS2)), data = data.ef, 
               size = 1, colour = "black",
               arrow = arrow(length = unit(0.2, "cm"))) +
  geom_text_repel(aes(x = NMDS1, y = NMDS2), data = data.ef, 
                  label = row.names(data.ef), 
                  colour = "black",
                  nudge_x = .6,
                  min.segment.length = .1,
                  max.overlaps = Inf,
                  segment.linetype = 4,
                  segment.color = "grey80",
                  direction = "y") +
  scale_color_manual(values = c("grey20","grey40","grey60","grey80")) +
  scale_shape_manual(values = c("square", "circle", "diamond", "triangle"))+
  annotate("text", x = -.35, y = .65, label = "2D stress = 0.29") +
  coord_fixed(ratio = 1, xlim = c(-.6, .8), ylim = c(-.9, .7)) +
  guides(color = guide_legend(reverse = F, title = ""), 
         shape = guide_legend(reverse = F, title = "")) +
  themeMB()

### Save ###
setwd(here("outputs/figures"))
ggsave("figure_2_ordination_sites_exposition1_(800dpi_10x10cm).tiff",
       dpi = 800, width = 10, height = 10, units = "cm")

### Ordination 2 ####
ggplot(data = data.scores) +
  geom_density_2d_filled(aes(x = NMDS1, y = NMDS2), data = data.scores,
                         alpha = .5,
                         contour_var = "ndensity",
                         show.legend = F) +
  geom_path(aes(x = NMDS1, y = NMDS2), data = df_ell,
            size = .5, linetype = 1) +
  geom_point(aes(x = NMDS1, y = NMDS2), data = data.scores, 
             size = 1) +
  facet_wrap(~ variable) +
  annotate("text", x = -.3, y = .6, size = 2, label = "2D stress = 0.29") +
  coord_fixed(ratio = 1, xlim = c(-.6, .8), ylim = c(-.9, .7)) +
  themeMB()

### Save ###
setwd(here("outputs/figures"))
ggsave("figure_2_ordination_sites_exposition2_(800dpi_16x8cm).tiff",
       dpi = 800, width = 16, height = 8, units = "cm")
