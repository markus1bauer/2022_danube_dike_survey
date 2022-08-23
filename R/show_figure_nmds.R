# Show Figure 2 NMDS ####
# Markus Bauer



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Preparation ###############################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


### Packages ###
library(here)
library(tidyverse)
library(vegan)
library(ggrepel)

### Start ###
rm(list = ls())
setwd(here("data", "processed"))

### Load data ###
sites <- read_csv("data_processed_sites.csv", col_names = TRUE,
                  na = c("na", "NA", ""), col_types =
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

#### a Choosen model ----------------------------------------------------------
set.seed(1)
(ordi <- metaMDS(species, dist = "bray", binary = FALSE,
                 try = 99, previous.best = TRUE, na.rm = TRUE))
### ef_vector ###
(ef <- envfit(ordi ~  targetRichness + eveness + graminoidCovratio, 
                      data = sites, permu = 999, na.rm = T))
rownames(ef$vectors$arrows) <- c("target species\nrichness",
                                 "eveness",
                                 "graminoid coverage")



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Plot ######################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


### * Functions ####
theme_mb <- function() {
  theme(
    panel.background = element_rect(fill = "white"),
    text = element_text(size = 9, color = "black"),
    strip.text = element_text(size = 10),
    axis.text = element_text(angle = 0, hjust = 0.5, size = 9,
                             color = "black"),
    axis.title = element_text(angle = 0, hjust = 0.5, size = 9,
                              color = "black"),
    axis.line = element_line(),
    legend.key = element_rect(fill = "white"),
    legend.position = "bottom",
    legend.margin = margin(0, 0, 0, 0, "cm"),
    plot.margin = margin(0, 0, 0, 0, "cm")
  )
}

veganCovEllipse <- function(cov, center = c(0, 0), scale = 1, npoints = 100) {
  theta <- (0:npoints) * 2 * pi / npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Plot ######################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


ellipses <- tibble()

data1 <-  sites %>%
  mutate(group_type = str_c(surveyYear_fac, exposition, targetType,
                            sep = "."),
         group_type = factor(group_type))

for(group in levels(data1$group_type)) {
  
  data2 <- data1 %>%
    filter(group_type == group) %>%
    with(
      cov.wt(
        cbind(NMDS1, NMDS2),
        wt = rep(1 / length(NMDS1), length(NMDS1))
      )
    )
  
  ellipses <- 
    veganCovEllipse(cov = data2$cov, center = data2$center) %>%
    as_tibble() %>%
    bind_cols(group = group) %>%
    bind_rows(ellipses)
  
  data <- ellipses %>%
    separate(
      group,
      sep = "\\.",
      c("surveyYear_fac", "exposition", "targetType")
    )
  
}


(graph_a <- ggplot() +
    geom_point(
      aes(y = NMDS2, x = NMDS1, color = surveyYear_fac,
          shape = surveyYear_fac, alpha = surveyYear_fac),
      data = sites,
      cex = 2
    ) +
    geom_path(
      aes(x = NMDS1, y = NMDS2, color = surveyYear_fac),
      data = data,
      size = 1
    ) +
    facet_grid(
      exposition ~ targetType,
      labeller = as_labeller(
        c(south = "South", north = "North",
          "dry_grassland" = "Dry grassland", "hay_meadow" = "Hay meadow")
      )
    ) +
    coord_fixed() +
    scale_color_manual(
      values = c(
        "orange1", "firebrick2", "deeppink3", "mediumpurple4",
        "royalblue", "black"
      )
    ) +
    scale_shape_manual(
      values = c(
        "circle open", "circle open", "circle open", "circle open",
        "square", "square open"
      )
    ) +
    scale_alpha_manual(values = c(.3, .3, .3, .3, .7, .6)) +
    labs(
      x = "NMDS1", y = "NMDS2", fill = "", color = "", shape = "", alpha = ""
    ) +
    theme_mb())


### Save ###
ggsave(here("outputs", "figures", "figure_nmds_800dpi_16.5x16cm.tiff"),
       dpi = 800, width = 16.5, height = 16, units = "cm")
