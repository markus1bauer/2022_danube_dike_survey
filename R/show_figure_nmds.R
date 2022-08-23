# Show Figure NMDS ####
# Markus Bauer



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Preparation ###############################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


### Packages ###
library(here)
library(tidyverse)
library(vegan)

### Start ###
rm(list = ls())
setwd(here("data", "processed"))

### * Load data sites ####
sites_dikes <- read_csv("data_processed_sites_spatial_nmds.csv",
                        col_names = TRUE, na = c("na", "NA", ""), col_types =
                          cols(
                            .default = "?",
                            id = "f"
                          )) %>%
  select(id, survey_year, target_richness) %>%
  mutate(survey_year_factor = as_factor(survey_year),
         target_richness_group = if_else(
           target_richness < 10, "<10", if_else(
             target_richness >= 10 & target_richness < 20, "10-19", if_else(
               target_richness >= 20 & target_richness < 30, "20-29", if_else(
                 target_richness >= 30, ">30", "warning"
                 )
               )
             )
           ),
         target_richness_group = fct_relevel(target_richness_group,
                                             "<10", "10-19", "20-29", ">30"))

sites_splot <- read_csv("data_processed_sites_splot.csv", col_names = TRUE,
                        na = c("na", "NA", ""), col_types =
                          cols(
                            .default = "?"
                          ))

sites <- sites_dikes %>%
  bind_rows(sites_splot) %>%
  mutate(esy = if_else(
    is.na(esy), "Dike plots", if_else(
      esy == "E12a", "Dry grassland", if_else(
        esy == "E22", "Hay meadow", "warning"
      )
    
  )))

### * Load data species ####
species_dikes <- read_csv("data_processed_species.csv", col_names = TRUE,
                          na = c("na", "NA", ""), col_types =
                            cols(
                              .default = "d",
                              name = "f"
                            ))

species_splot <- read_csv("data_processed_species_splot.csv", col_names = TRUE,
                          na = c("na", "NA", ""), col_types =
                            cols(
                              .default = "?"
                            ))

species <- species_dikes %>%
  full_join(species_splot, by = "name") %>%
  mutate(across(where(is.numeric), ~replace(., is.na(.), 0))) %>%
  pivot_longer(cols = -name, names_to = "id", values_to = "value") %>%
  pivot_wider(names_from = "name", values_from = "value") %>%
  arrange(id) %>%
  semi_join(sites, by = "id") %>%
  column_to_rownames("id")

rm(list = setdiff(ls(), c("sites", "species")))

#### * Choosen model ####
set.seed(1)
(ordi <- metaMDS(species, dist = "bray", binary = FALSE,
                 try = 99, previous.best = TRUE, na.rm = TRUE))



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
    legend.position = "right",
    legend.margin = margin(0, 0, 0, 0, "cm"),
    plot.margin = margin(0, 0, 0, 0, "cm")
  )
}

veganCovEllipse <- function(cov, center = c(0, 0), scale = 1, npoints = 100) {
  theta <- (0:npoints) * 2 * pi / npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}


### * Preparation ####
ellipses <- tibble()

data_nmds <-  sites %>%
  select(id, esy, target_richness_group) %>% # modify group
  mutate(group_type = as_factor(esy), # modify group
         NMDS1 = ordi$points[,1],
         NMDS2 = ordi$points[,2]) %>%
  group_by(group_type) %>%
  mutate(mean1 = mean(NMDS1),
         mean2 = mean(NMDS2))

for(group in levels(data_nmds$group_type)) {
  
  ellipses_calc <- data_nmds %>%
    filter(group_type == group) %>%
    with(
      cov.wt(
        cbind(NMDS1, NMDS2),
        wt = rep(1 / length(NMDS1), length(NMDS1))
      )
    )
  
  ellipses <- 
    veganCovEllipse(cov = ellipses_calc$cov, center = ellipses_calc$center) %>%
    as_tibble() %>%
    bind_cols(group_type = group) %>%
    bind_rows(ellipses)
  
  data <- ellipses
  
}


(graph_a <- ggplot() +
    geom_point(
      aes(y = NMDS2, x = NMDS1, color = target_richness_group,
          shape = group_type),
      data = data_nmds,
      cex = 2
    ) +
    geom_path(
      aes(x = NMDS1, y = NMDS2, linetype = group_type),
      data = data,
      size = 1,
      show.legend = FALSE
    ) +
    geom_text(
      aes(x = mean1, y = mean2, label = group_type),
      data = data_nmds
    ) +
    coord_fixed() +
    scale_color_brewer(na.value = "black") +
    scale_linetype_manual(values = c(1, 1, 1)) +
    labs(
      x = "NMDS1", y = "NMDS2", color = "Target species\nrichness", shape = ""
    ) +
    theme_mb())


### Save ###
ggsave(here("outputs", "figures", "figure_nmds_800dpi_16.5x16cm.tiff"),
       dpi = 800, width = 16.5, height = 16, units = "cm")
