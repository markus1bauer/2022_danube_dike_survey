# Show Figure ordination ####
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
                       position = col_factor(),
                       dataset = col_factor(),
                       side = col_factor(),
                       exposition = col_factor(),
                       phosphorousClass = col_factor(c("A", "B", "C", "D", "E")),
                       biotopeType = col_factor(),
                       baykompv = col_factor(),
                       ffh = col_factor(c("6210", "6510", "non-FFH")),
                       min8 = col_factor(),
                       min9 = col_factor()
                     ))
data <- sites %>%
  select(id, plot, ffh, NMDS1, NMDS2) %>%
  mutate(group = ffh)

species <- read_csv2("data_processed_species.csv", col_names = T, na = "na", col_types = 
                       cols(
                         .default = col_double(),
                         name = col_factor()
                       )) %>%  
  pivot_longer(-name, "id", "value") %>%
  pivot_wider(id, name) %>%
  arrange(id) %>%
  column_to_rownames("id")

#### a Choosen model ----------------------------------------------------------------------------------------
set.seed(1)
ordi <- metaMDS(species, try = 99, previous.best = T, na.rm = T)
### ef_vector2 ###
(ef <- envfit(ordi ~  vegetationCov + plotAge + topsoilDepth + phosphorous + NtotalConc + NtotalPerc + cnRatio + sandPerc + ffh6510Richness + ffh6210Richness, 
                      data = sites, permu = 999, na.rm = T))
rm(species)



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
  data.ell <- rbind(data.ell,cbind(as.data.frame(with(data[data$group == g,], veganCovEllipse(cov.wt(cbind(NMDS1, NMDS2), wt = rep(1 / length(NMDS1), length(NMDS1)))$cov, center = c(mean(NMDS1), mean(NMDS2))))), group = g))
}
### Environmental factors ###
data.ef <- as.data.frame(ef$vectors$arrows * ((sqrt(ef$vectors$r))) * 1) %>%
  rownames_to_column(var = "envFactors") %>%
  as_tibble()
### Plot ###
set.seed(1)
ggplot(data = data) +
  geom_path(aes(x = NMDS1, y = NMDS2, colour = group), data = data.ell,
            size = 1, linetype = 1) +
  geom_point(aes(x = NMDS1, y = NMDS2, shape = group, colour = group), data = data, 
             size = 3) +
  geom_segment(aes(x = 0, xend = (NMDS1), y = 0, yend = (NMDS2)), data = data.ef, 
               size = .5, colour = "black",
               arrow = arrow(length = unit(0.2, "cm"))) +
  geom_text_repel(aes(x = NMDS1, y = NMDS2, label = envFactors), data = data.ef, 
                   colour = "black",
                   min.segment.length = .1,
                   max.overlaps = Inf,
                   segment.linetype = 4,
                   segment.color = "grey80") +
  scale_color_manual(values = c("grey20","grey40","grey60")) +
  scale_shape_manual(values = c("square", "circle", "circle open"))+
  annotate("text", x = -.35, y = .65, label = "2D stress = 0.29") +
  coord_fixed(ratio = 1, xlim = c(-.8, .8), ylim = c(-.9, .7)) +
  guides(colour = guide_legend(reverse = F, title="FFH type:"), 
         shape = guide_legend(reverse = F, title = "FFH type:")) +
  themeMB()

### Save ###
setwd(here("outputs/figures"))
ggsave("figure_2_ordination_sites_ffh_normal_(800dpi_16x16cm).tiff",
       dpi = 800, width = 16, height = 16, units = "cm")
