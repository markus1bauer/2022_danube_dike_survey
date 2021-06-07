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
                       .default = col_guess(),
                       id = col_factor(),
                       block = col_factor(),
                       plot = col_factor(),
                       ffh = col_factor(c("6210", "6510", "non-FFH")),
                       changeType = col_factor(c("better", "change", "worse", "any-FFH", "FFH6510", "non-FFH"))
                     )) %>%
  select(id, block, plot, ffh, changeType, NMDS1, NMDS2) 
data <- sites %>%
  select(id, plot, ffh, NMDS1, NMDS2) %>%
  mutate(group = ffh)
sites <- sites %>%
  filter(!is.na(changeType))


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

### Create Ellipses ###
data.ell <- data.frame()
for(g in levels(data$group)){
  data.ell <- rbind(data.ell,cbind(as.data.frame(with(data[data$group == g,], veganCovEllipse(cov.wt(cbind(NMDS1, NMDS2), wt = rep(1 / length(NMDS1), length(NMDS1)))$cov, center = c(mean(NMDS1), mean(NMDS2))))), group = g))
}


### 2 Change ordination #####################################################################################

ggplot() +
  geom_polygon(aes(x = NMDS1, y = NMDS2, fill = changeType, group = plot), data = sites,
               size = .5, linetype = 1, alpha = .5, show.legend = F) +
  geom_path(aes(x = NMDS1, y = NMDS2, colour = group), data = data.ell,
            size = 1, linetype = 1) +
  geom_point(aes(x = NMDS1, y = NMDS2, shape = ffh, colour = ffh), data = sites, 
             size = 1) +
  facet_wrap(~changeType) +
  scale_color_manual(values = c("grey20","grey35","grey50")) +
  scale_shape_manual(values = c("square", "circle", "circle open"))+
  annotate("text", x = -.3, y = .65, size = 2, label = "2D stress = 0.29") +
  coord_fixed(ratio = 1, xlim = c(-.6, .8), ylim = c(-.9, .7)) +
  labs(shape = "FFH type", colour = "FFH type", fill = "") +
  themeMB()

### Save ###
setwd(here("outputs/figures"))
ggsave("figure_2_ordination_ffh_plots_(800dpi_14x12cm).tiff",
       dpi = 800, width = 14, height = 12, units = "cm")
