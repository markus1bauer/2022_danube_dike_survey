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

#### a Load data ----------------------------------------------------------------------------------------
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
                     )        
)

species <- read_csv2("data_processed_species.csv", col_names = T, na = "na", col_types = 
                       cols(
                         .default = col_double(),
                         name = col_factor()
                       )) %>%  
  pivot_longer(-name, "id", "value") %>%
  pivot_wider(id, name) %>%
  arrange(id) %>%
  column_to_rownames("id")

#### b Load functions ----------------------------------------------------------------------------------------
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

#### c Choosen model for NMDS ----------------------------------------------------------------------------------------
ordi <- metaMDS(species, try = 99, previous.best = T, na.rm = T)
data <- ordi %>%
  scores() %>%
  as.data.frame() %>%
  rownames_to_column(var = "id") %>%
  as_tibble() %>%
  left_join(sites, by = "id") %>%
  select(id, plot, NMDS1, NMDS2, ffh, surveyYear) %>%
  mutate(group = ffh, .keep = "unused")
### Create Ellipses ###
df_ell <- data.frame()
for(g in levels(data$group)){
  df_ell <- rbind(df_ell,cbind(as.data.frame(with(data[data$group == g,], veganCovEllipse(cov.wt(cbind(NMDS1, NMDS2), wt = rep(1 / length(NMDS1), length(NMDS1)))$cov, center = c(mean(NMDS1), mean(NMDS2))))), group = g))
}
data.ell <- df_ell

#### d Choosen model for environmental factors ----------------------------------------------------------------------------------------
(ef <- envfit(ordi ~  vegetationCov + plotAge + topsoilDepth + phosphorous + NtotalConc + NtotalPerc + cnRatio + sandPerc + ffh6510Richness + ffh6210Richness, 
              data = sites, permu = 999, na.rm = T))
data.ef <- as.data.frame(ef$vectors$arrows * ((sqrt(ef$vectors$r))) * 1) %>%
  rownames_to_column(var = "envFactors") %>%
  as_tibble()


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Plot ################################################################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

### 1 Density ordination #####################################################################################

ggplot(data = data) +
  geom_density_2d_filled(aes(x = NMDS1, y = NMDS2), data = data,
                         alpha = .5,
                         contour_var = "ndensity",
                         show.legend = F) +
  geom_path(aes(x = NMDS1, y = NMDS2), data = data.ell,
            size = .3, linetype = 1) +
  geom_point(aes(x = NMDS1, y = NMDS2), data = data, 
             size = .3) +
  facet_wrap(~ group, ncol = 2) +
  annotate("text", x = -.2, y = -.75, size = 2, label = "2D stress = 0.29") +
  coord_fixed(ratio = 1, xlim = c(-.6, .8), ylim = c(-.9, .7)) +
  themeMB()

### Save ###
setwd(here("outputs/figures"))
ggsave("figure_2_ordination_sites_ffh_density_(800dpi_16x8cm).tiff",
       dpi = 800, width = 8, height = 8, units = "cm")
