# Beta diversity on dike grasslands
# Table A7 ####

# Markus Bauer
# 2022-09-13



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Preparation ###############################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



### Packages ###
library(here)
library(tidyverse)
library(gt)

### Start ###
rm(list = ls())

### Load data ###
data <- read_csv(
  here("outputs", "statistics", "vegetation_classification.csv"),
  col_names = TRUE, na = c("na", "NA", " "), col_types = cols(.default = "?")
) %>%
  mutate(
    esy = as_factor(esy),
    esy = fct_relevel(esy, "no class", after = 4),
    esy = fct_relevel(esy, "R1A", after = 2),
    across(where(is.numeric), ~ replace(., is.na(.), 0))
  ) %>%
  arrange(esy)



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Plot ######################################################################
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
     column_labels.font.weight = "bold",
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
     locations = cells_body(columns = "esy"),
     style = cell_text(align = "left")
     ) %>%
   tab_style(
     locations = cells_column_labels(columns = "esy"),
     style = cell_text(align = "left")
     ) %>%
### 3 Column labels ####
   cols_label(
     "esy" = "Class"
     ))

### Save ###
# gtsave(
#   table, here("outputs", "tables", "table_a7_vegetation_classification.png")
#   )
