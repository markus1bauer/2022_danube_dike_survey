# Beta diversity on dike grasslands
# Table A5 ####
# Markus Bauer
# 2022-09-08



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Preparation ###############################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



### Packages ###
library(here)
library(tidyverse)
library(gt)

### Start ###
rm(list = ls())
setwd(here("outputs", "statistics"))

### Load data ###
data <- read_csv("pca_soil.csv",
                 col_names = TRUE,
                 na = c("na", "NA", " "), col_types =
                   cols(
                     .default = "?"
                   )) %>%
  select(-PC4) %>%
  mutate(across(where(is.numeric), ~ round(., digits = 2)),
         variables = fct_relevel(
           variables, c(
           "ph",
           "topsoil_depth",
           "clay",
           "silt",
           "sand",
           "calciumcarbonat",
           "humus",
           "cn_ratio",
           "n_total_ratio",
           "n_total_concentration",
           "phosphorus",
           "potassium",
           "magnesium"
         )),
         variables = fct_recode(
           variables,
           "pH" = "ph",
           "Topsoil depth" = "topsoil_depth",
           "Clay" = "clay",
           "Silt" = "silt",
           "Sand" = "sand",
           "Humus" = "humus",
           "CaCO3" = "calciumcarbonat",
           "C:N ratio" = "cn_ratio",
           "N" = "n_total_ratio",
           "N concentration" = "n_total_concentration",
           "P" = "phosphorus",
           "K" = "potassium",
           "Mg2+" = "magnesium"
         )) %>%
  arrange(variables) %>%
  mutate(unit = c(
    "",
    "cm",
    "wt%",
    "wt%",
    "wt%",
    "wt%",
    "wt%",
    "",
    "wt%",
    "kg/m²",
    "mg/100g",
    "mg/100g",
    "mg/100g",
    "", "", ""
  )) %>%
  relocate(variables, unit)



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Plot with gt ##############################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



(table <- data %>%
   gt() %>%

### 1 General settings ####
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

### 2 Alignment of specific columns ####
 tab_style(
   locations = cells_body(columns = "variables"),
   style = cell_text(align = "left")
 ) %>%
   tab_style(
     locations = cells_column_labels(columns = "unit"),
     style = cell_text(align = "center")
   ) %>%
   tab_style(
     locations = cells_body(columns = "unit"),
     style = cell_text(align = "center")
   ) %>%

### 3 Column labels ####
 cols_label(
   variables = md(""),
   unit = md("**Unit**"),
   PC1 = md("**PC1**"),
   PC2 = md("**PC2**"),
   PC3 = md("**PC3**")
 ) %>%

### 4 Highlight certain cells ####
 tab_style(
   style = list(
     cell_fill(color = "grey"),
     cell_text(weight = "bold")
     ),
   locations = cells_body(
     columns = PC1,
     rows = PC1 >= 1.2 & PC1 < 2 | PC1 <= -1.16
     )
   ) %>%
   tab_style(
     style = list(
       cell_fill(color = "grey"),
       cell_text(weight = "bold")
     ),
     locations = cells_body(
       columns = PC2,
       rows = PC2 >= .91 & PC2 < 2 | PC2 <= -1
     )
   ) %>%
   tab_style(
     style = list(
       cell_fill(color = "grey"),
       cell_text(weight = "bold")
     ),
     locations = cells_body(
       columns = PC3,
       rows = PC3 >= .9 & PC3 < 1.15 | PC3 <= -.88
     )
   ) %>%

### 5 Make subscripts ####
 text_transform(
   locations = cells_body(
     columns = c(variables),
     rows = c(6, 9, 10)
   ),
   fn = function(x) {
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
     fn = function(x) {
       x <- c("2+")
       text <- c("Mg")
       glue::glue("{text}<sup>{x}</sup>")
     }
   ))


### Save ###
gtsave(table, here("outputs", "tables", "table_a5_pca_soil.png"))
