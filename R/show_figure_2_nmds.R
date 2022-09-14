b# Beta diversity on dike grasslands
# Plot Fig. 2 ####
# Markus Bauer
# 2022-09-12



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

### Functions ###
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
    legend.text = element_text(size = 10),
    plot.margin = margin(0, 0, 0, 0, "cm")
  )
}

veganCovEllipse <- function(cov, center = c(0, 0), scale = 1, npoints = 100) {
  theta <- (0:npoints) * 2 * pi / npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

#### * Load data sites ####

sites_dikes <- read_csv("data_processed_sites_spatial.csv",
                        col_names = TRUE, na = c("na", "NA", ""), col_types =
                          cols(
                            .default = "?",
                            id = "f"
                          )) %>%
  select(id, survey_year, orientation, exposition, esy,
         species_richness, eveness, shannon,
         ellenberg_richness, ellenberg_cover_ratio,
         accumulated_cover, graminoid_cover_ratio, ruderal_cover) %>%
  mutate(survey_year_factor = as_factor(survey_year))

sites_splot <- read_csv("data_processed_sites_splot.csv", col_names = TRUE,
                        na = c("na", "NA", ""), col_types =
                          cols(
                            .default = "?"
                          ))

sites <- sites_dikes %>%
  bind_rows(sites_splot) %>%
  mutate(
    reference = as.character(survey_year),
    reference = if_else(
      esy == "E12a", "Dry grassland", if_else(
        esy == "E22", "Hay meadow", reference
      )
    ),
    esy = if_else(
      esy == "E22", "R22-ref", if_else(
        esy == "E12a", "R1A-ref", if_else(
          esy == "?", "no", if_else(
            esy == "R21", "R", esy
          )
        )
      )
    )
  ) %>%
  select(-givd_id, -longitude, -latitude)

#### * Load data species ####

species_dikes <- read_csv("data_processed_species.csv", col_names = TRUE,
                          na = c("na", "NA", ""), col_types =
                            cols(
                              .default = "d",
                              name = "f"
                            )) %>%
  select(name, all_of(sites_dikes$id))

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

rm(list = setdiff(ls(), c("sites", "species", "theme_mb", "veganCovEllipse")))

#### * Choosen model ####

set.seed(10)
(ordi <- metaMDS(species, binary = TRUE,
                 try = 99, previous.best = TRUE, na.rm = TRUE))

(data_envfit <- envfit(ordi ~ graminoid_cover_ratio + ruderal_cover +
                         ellenberg_richness,
                      data = sites,
                      perm = 999,
                      na.rm = TRUE))



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Plot ######################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



### * Preparation ####

data_envfit <- data_envfit %>%
  scores(display = "vectors") %>%
  as_tibble(rownames = NA) %>%
  rownames_to_column(var = "variable") %>%
  mutate(
    variable = as_factor(variable),
    variable = fct_recode(
      variable,
      "Graminoid cover" = "graminoid_cover_ratio",
      "Ruderal cover" = "ruderal_cover",
      "Specialist richness" = "ellenberg_richness"
      )
    )

ellipses <- tibble()

data_nmds <-  sites %>%
  select(id, esy) %>% # modify group
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
  
  data_ellipses <- ellipses
  
}

#### * Plot ####

(graph_a <- ggplot() +
   geom_point(
     aes(y = NMDS2, x = NMDS1, color = group_type, shape = group_type),
     data = data_nmds,
     cex = 2
   ) +
   geom_path(
     aes(x = NMDS1, y = NMDS2, linetype = group_type, color = group_type),
     data = data_ellipses %>% filter(group_type != "no"),
     size = 1,
     show.legend = FALSE
   ) +
   ggrepel::geom_label_repel(
     aes(x = NMDS1, y = NMDS2, label = variable),
     data = data_envfit %>% filter(NMDS2 < 0),
     fill = alpha("white", .7),
     size = 3,
     nudge_y = -.1,
     min.segment.length = Inf
   ) +
   ggrepel::geom_label_repel(
     aes(x = NMDS1, y = NMDS2, label = variable),
     data = data_envfit %>% filter(NMDS2 > 0),
     fill = alpha("white", .7),
     size = 3,
     nudge_y = .1,
     min.segment.length = Inf
   ) +
   geom_segment(
     data = data_envfit,
     aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
     arrow = arrow(length = unit(0.25, "cm")),
     color = "black",
     size = 1
   ) +
   geom_label(
     aes(x = mean1, y = mean2, label = group_type, fill = group_type),
     data = data_nmds %>%
       filter(group_type != "no") %>%
       group_by(group_type) %>%
       slice(1),
     color = "white",
     size = 3,
     show.legend = FALSE
   ) +
   coord_fixed() +
   scale_color_manual(
     labels = c(
       "R22-ref" = "R22-ref:\nHay meadow reference",
       "R1A-ref" = "R1A-ref:\nDry grassland reference",
       "R22" = "R22: Hay meadow",
       "R1A" = "R1A: Dry grassland",
       "R" = "R: General grassland",
       "V38" = "V38:\nDry anthropogenic\nvegetation",
       "no" = "no classification"
     ),
     values = c(
       "R22-ref" = "#00BFC4",
       "R1A-ref" = "#F8766D",
       "R22" = "#00BFC4",
       "R1A" = "#F8766D",
       "R" = "grey30",
       "V38" = "#C77CFF",
       "no" = "grey90"
       )
     ) +
   scale_fill_manual(
     values = alpha(c(
       "R22-ref" = "#00BFC4",
       "R1A-ref" = "#F8766D",
       "R22" = "#00BFC4",
       "R1A" = "#F8766D",
       "R" = "grey30",
       "V38" = "#C77CFF",
       "no" = "grey30"
     ), alpha = 0.6)
   ) +
   scale_shape_manual(
     labels = c(
       "R22-ref" = "R22-ref:\nHay meadow reference",
       "R1A-ref" = "R1A-ref:\nDry grassland reference",
       "R22" = "R22: Hay meadow",
       "R1A" = "R1A: Dry grassland",
       "R" = "R: General grassland",
       "V38" = "V38:\nDry anthropogenic\nvegetation",
       "no" = "no classification"
     ),
     values = c(
       "R22-ref" = 1,
       "R1A-ref" = 1,
       "R22" = 16,
       "R1A" = 16,
       "R" = 3,
       "V38" = 15,
       "no" = 4
     )
   ) +
   scale_linetype_manual(values = c(1, 1, 1, 1, 1, 1, 1)) +
   labs(x = "NMDS1", y = "NMDS2", shape = "", color = "") +
   theme_mb() +
   theme(legend.key.height = unit(.9, "cm")))


### Save ###
ggsave(here("outputs", "figures", "figure_2_800dpi_16.5x11cm.tiff"),
       dpi = 800, width = 16.5, height = 11, units = "cm")
