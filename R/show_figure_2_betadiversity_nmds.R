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
sites <- read_csv("data_processed_sites.csv", col_names = T, na = c("na", "NA", ""), col_types = 
                    cols(
                      .default = "d",
                      id = "f",
                      locationAbb = "f",
                      block = "f",
                      plot = "f",
                      exposition = "f"
                    )) %>%
  select(id, plot, block, locationAbb, surveyYear, exposition, vegetationCov, targetCov, graminoidCovratio, speciesRichness, targetRichness, shannon, eveness, accumulatedCov) %>%
  mutate(surveyYearF = as_factor(surveyYear)) %>%
  filter(accumulatedCov > 0)

species <- read_csv("data_processed_species.csv", col_names = T, na = c("na", "NA", ""), col_types = 
                      cols(
                        .default = "d",
                        name = "f"
                      )) %>%  
  mutate(across(where(is.numeric), ~replace(., is.na(.), 0))) %>%
  pivot_longer(-name, "id", "value") %>%
  pivot_wider(id, name) %>%
  arrange(id) %>%
  semi_join(sites, by = "id") %>%
  column_to_rownames("id")

#### a Choosen model ----------------------------------------------------------------------------------------
set.seed(1)
(ordi <- metaMDS(species, dist = "bray", binary = F,
                 try = 99, previous.best = T, na.rm = T))
### ef_vector ###
(ef <- envfit(ordi ~  targetRichness + eveness + graminoidCovratio, 
                      data = sites, permu = 999, na.rm = T))
rownames(ef$vectors$arrows) <- c("target species\nrichness", "eveness", "graminoid coverage")



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
    legend.background = element_rect(fill = "white"),
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

data.ef <- as.data.frame(ef$vectors$arrows * ((sqrt(ef$vectors$r)))) * 1.5
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
  geom_text_repel(data = subset(data.ef, variables == "eveness"),
                  aes(x = NMDS1, y = NMDS2, label = variables),
                  nudge_y = .3,
                  #nudge_x = .1,
                  segment.color = "grey60") +
  geom_text_repel(data = subset(data.ef, variables == "target species\nrichness"),
                  aes(x = NMDS1, y = NMDS2, label = variables),
                  nudge_y = -.1,
                  segment.color = "grey60") +
  geom_text_repel(data = subset(data.ef, variables == "graminoid coverage"),
                  aes(x = NMDS1, y = NMDS2, label = variables),
                  nudge_y = -.3,
                  segment.color = "grey60") +
  coord_fixed(ratio = 1, xlim = c(-.8, .9), ylim = c(-.7, .6)) +
  scale_x_continuous(breaks = seq(-100, 100, 0.5)) +
  scale_y_continuous(breaks = seq(-100, 100, 0.5)) +
  scale_color_manual(breaks = c("north", "south", "west", "east"), values = c("grey10","grey40","grey70","grey70")) +
  scale_shape_manual(breaks = c("north", "south", "west", "east"), values = c("square", "circle", "plus", "cross"))+
  annotate("text", x = .65, y = .6, label = "2D stress = .28") +
  guides(color = guide_legend(reverse = F, title = ""), 
         shape = guide_legend(reverse = F, title = "")) +
  themeMB()

### Save ###
ggsave(here("outputs/figures/figure_2_betadiversity_nmds_(800dpi_16x10cm).tiff"),
       dpi = 800, width = 16, height = 10, units = "cm")
