# Show Figure 2 NMDS ####
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
species <- read_csv("data_processed_species.csv", col_names = T, na = "na", col_types = 
                      cols(
                        .default = "d",
                        name = "f"
                      )) %>%  
  pivot_longer(-name, "id", "value") %>%
  pivot_wider(id, name) %>%
  arrange(id) 

sites <- read_csv("data_processed_sites.csv", col_names = T, na = "na", col_types = 
                    cols(
                      .default = "d",
                      id = "f",
                      ffh = "f"
                    )) %>%
  select(id, NMDS1, NMDS2, ffh, vegetationCov, targetCov, graminoidCovratio, speciesRichness, targetRichness, shannon, eveness) %>%
  semi_join(species, by = "id") %>%
  rename(group = ffh)
data <- sites

species <- species %>%
  column_to_rownames("id")

#### a Choosen model ----------------------------------------------------------------------------------------
set.seed(1)
(ordi <- metaMDS(species, dist = "bray", binary = F,
                 try = 99, previous.best = T, na.rm = T))
### ef_vector2 ###
(ef <- envfit(ordi ~  targetRichness + eveness + graminoidCovratio, 
                      data = sites, permu = 999, na.rm = T))
rownames(ef$vectors$arrows) <- c("target species richness", "eveness", "graminoid coverage")



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Plot ################################################################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


themeMB <- function(){
  theme(
    panel.background = element_rect(fill = "white"),
    text  = element_text(size = 11, color = "black"),
    axis.line = element_line(),
    legend.title = element_text(size = 11,  face = "bold"),
    legend.key = element_rect(fill = NA),
    legend.background = element_rect(fill = NA),
    legend.position = c(0.05, 0.05),
    legend.justification = c("left", "bottom"),
    legend.direction = "vertical",
    legend.margin = margin(-.4, 0, 0, 0, "cm"),
    legend.box = "horizontal",
    plot.margin = margin(0, 0, 0, 0, "cm")
  )
}
veganCovEllipse <- function(cov, center = c(0,0), scale = 1, npoints = 100)
{
  theta <- (0:npoints) * 2 * pi / npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}


### 1 Normal ordination #####################################################################################

### Create Ellipses ###
data.ell <- data.frame()
for(g in levels(data$group)){
  data.ell <- rbind(data.ell, cbind(as.data.frame(with(data[data$group == g,], veganCovEllipse(cov.wt(cbind(NMDS1, NMDS2), wt = rep(1 / length(NMDS1), length(NMDS1)))$cov, center = c(mean(NMDS1), mean(NMDS2))))), group = g))
}
### Environmental factors ###
data.ef <- as.data.frame(ef$vectors$arrows * ((sqrt(ef$vectors$r))) * 1) %>%
  rownames_to_column(var = "envFactors") %>%
  as_tibble()
### Plot ###
set.seed(1)
ggplot(data = sites) +
  geom_point(aes(x = NMDS1, y = NMDS2, color = group, shape = group, size = speciesRichness), data = sites, 
             stroke = .7) +
  #geom_path(aes(x = NMDS1, y = NMDS2, colour = group), data = data.ell,
  #          size = 1, linetype = 1) +
  geom_segment(aes(x = 0, xend = (NMDS1), y = 0, yend = (NMDS2)), data = data.ef, 
               size = .5, colour = "black",
               arrow = arrow(length = unit(0.2, "cm"))) +
  geom_text_repel(aes(x = NMDS1, y = NMDS2, label = envFactors), data = data.ef, 
                  colour = "black",
                  min.segment.length = .1,
                  max.overlaps = Inf,
                  segment.linetype = 4,
                  segment.color = "grey95") +
  scale_size(breaks = c(10, 20, 30, 40)) +
  scale_color_manual(values = c("black","grey40","grey95")) +
  scale_shape_manual(values = c("circle open", "circle", "circle")) +
  annotate("text", x = -.65, y = .65, label = "2D stress = 0.29") +
  coord_fixed(ratio = 1, xlim = c(-.8, .8), ylim = c(-.9, .7)) +
  guides(size = guide_legend(reverse = F, title = "Species richness"),
         colour = guide_legend(reverse = F, title="FFH type"), 
         shape = guide_legend(reverse = F, title = "FFH type")) +
  themeMB()

### Save ###
setwd(here("outputs/figures"))
ggsave("figure_4_betadiversity_nmds_(800dpi_16x10cm).tiff",
       dpi = 800, width = 16, height = 10, units = "cm")

