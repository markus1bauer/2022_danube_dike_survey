# Beta diversity on dike grasslands
# Table 3 ####
# Markus Bauer
# 2022-01-11
# Citation: 
## Bauer M, Huber J, Kollmann J (submitted) 
## Balanced turnover is a main aspect of biodiversity on restored dike grasslands: not only deterministic environmental effects, but also non-directional year and site effects drive spatial and temporal beta diversity.
## Unpublished data.


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Preparation ################################################################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

### Packages ###
library(here)
library(tidyverse)
library(gt)

### Start ###
rm(list = ls())
setwd(here("outputs/statistics"))

### Load data ###
data <- read_csv("pca_soil.csv", col_names = T, na = c("na", "NA", " "), col_types = 
                     cols(
                       .default = "?")) %>%
  relocate(variables) %>%
  select(-PC4) %>%
  filter(variables != "Eigenvalues" & variables != "Proportion Explained" & variables != "Cumulative Proportion") %>%
  mutate(across(where(is.numeric), ~round(., digits = 2))) %>%
  mutate(variables = fct_relevel(variables, c("topsoilDepth",
                                              "sandPerc",
                                              "siltPerc",
                                              "clayPerc",
                                              "pH",
                                              "calciumcarbonatPerc",
                                              "humusPerc",
                                              "cnRatio",
                                              "NtotalPerc",
                                              "NtotalConc",
                                              "phosphorous",
                                              "potassium",
                                              "magnesium")),
         variables = fct_recode(variables, 
                                "Topsoil depth" = "topsoilDepth",
                                "Sand" = "sandPerc",
                                "Silt" = "siltPerc",
                                "Clay" = "clayPerc",
                                "Humus" = "humusPerc",
                                "CaCO3" = "calciumcarbonatPerc",
                                "C:N ratio" = "cnRatio",
                                "N" = "NtotalPerc",
                                "N concentration" = "NtotalConc",
                                "P" = "phosphorous",
                                "K" = "potassium",
                                "Mg2+" = "magnesium")) %>%
  arrange(variables)



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Plot with gt ################################################################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


(table <- data %>%
    gt() %>%

### General settings ####
    opt_table_lines(
       extent = "none"
    ) %>%
    tab_options( 
       table.font.style = "Arial",
       table.font.size = px(12),
       table.font.color = "black",
       data_row.padding = px(4),
       table.align = "left",
       column_labels.border.top.style = "solid",
       table_body.border.bottom.style = "solid",
       table_body.border.top.style = "solid",
       column_labels.border.top.color = "black",
       table_body.border.bottom.color = "black",
       table_body.border.top.color = "black",
       column_labels.border.top.width = px(2),
       table_body.border.bottom.width = px(2),
       table_body.border.top.width = px(1)
    ) %>%

### Alignment of first column ####
    tab_style(
       style = cell_text(
          align = "left"
       ),
       locations = cells_column_labels(
          columns = "variables"
       )
    ) %>%
    tab_style(
       style = cell_text(
          align = "left"
       ),
       locations = cells_body(
          columns = "variables"
       )
    ) %>%
   
### Heading ####
    tab_header( 
       title = ""
    ) %>%
    tab_style(
       style = cell_text(
          align = "left"
       ),
       locations = cells_title()
    ) %>%
   
### Spanners ####
    #tab_spanner( 
    # label = "Treatments",
    #columns = vars(
    # brickRatio, acid, soilFertility
    #)) %>%
    #styles of spanners
    #tab_style( 
    # style = cell_text(
    #  weight = "bold"
    #),
 #locations = cells_column_spanners(
 # spanners = T
 #  )) %>%
 #tab_style(
 # style = cell_borders(
 #  sides = c("bottom"),
 # color = "black",
 # style = "solid",
 # weight = px(1)
 #),
 #locations = cells_column_spanners(
 #    spanners = T
 #)) %>%

### Column labels ####
 cols_label( 
    variables = md(""),
    PC1 = md("**PC1 (43%)**"),
    PC2 = md("**PC2 (21%)**"),
    PC3 = md("**PC3 (10%)**"),
 ) %>%
  
### Footnote ####
    #tab_source_note( 
    #   source_note = "PCA for old dikes"
    #) %>%

### Highlight certain cells ####
    tab_style(
       style = list(
          cell_fill(color = "grey"),
          cell_text(weight = "bold")
       ),
       locations = cells_body(
          columns = PC1,
          rows = PC1 >= 1.2 & PC1 < 2 | PC1 <= -1.16
       )) %>%
    tab_style(
       style = list(
          cell_fill(color = "grey"),
          cell_text(weight = "bold")
       ),
       locations = cells_body(
          columns = PC2,
          rows = PC2 >= .9 & PC2 < 2 | PC2 <= -1
       )) %>%
    tab_style(
       style = list(
          cell_fill(color = "grey"),
          cell_text(weight = "bold")
       ),
       locations = cells_body(
          columns = PC3,
          rows = PC3 >= .9 & PC3 < 1.15 | PC3 <= -1
       )) %>%
  
  ### Make subscripts ####
   text_transform(
    locations = cells_body(
      columns = c(variables),
      rows = c(6, 9, 10)
    ),
    fn = function(x){
      x <- c("3", "total", "concentration")
      text <- c("CaCO", "N", "N")
      glue::glue("{text}<sub>{x}</sub>")
    }
  ) %>%
  text_transform(
    locations = cells_body(
      columns = c(variables),
      rows = c(13)
    ),
    fn = function(x){
      x <- c("2+")
      text <- c("Mg")
      glue::glue("{text}<sup>{x}</sup>")
    }
  )
   
)



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# C Save ################################################################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


gtsave(table, here("outputs/tables/table_3_pca_soil.png"))
