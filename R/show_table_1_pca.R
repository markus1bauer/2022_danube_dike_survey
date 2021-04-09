### Show table PCA ###
### Markus Bauer ###


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Preparation ################################################################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

### Packages ###
library(here)
library(tidyverse)
library(gt)

### Start ###
rm(list = ls())
setwd(here("data/processed"))

### Load data ###
data <- read_csv2("data_processed_pca.csv", col_names = T, na = "na", col_types = 
                     cols(
                       .default = col_guess())) %>%
   relocate(variables)



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Plot with gt ################################################################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


(table <- data %>%
    gt() %>%
    opt_table_lines(
       extent = "none"
    ) %>%
    tab_options( #general settings
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
    #general alignments
    tab_style( 
       style = cell_text(
          v_align = "top",
          align = "center"
       ),
       locations = cells_column_labels(
          columns = T
       )
    ) %>%
    tab_style(
       style = cell_text(
          align = "center",
          v_align = "middle"
       ),
       locations = cells_body(
          columns = T
       )
    ) %>%
    #alignment of first column
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
    #heading
    tab_header( 
       title = "Table 1"
    ) %>%
    tab_style(
       style = cell_text(
          align = "left"
       ),
       locations = cells_title()
    ) %>%
    #spanners
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
 #column labels
 cols_label( 
    variables = md("**Variable**")
 ) %>%
    #footnote
    tab_source_note( 
       source_note = "PCA for old dikes"
    ) %>%
    #Highlight certain cells   
    tab_style(
       style = list(
          cell_fill(color = "grey"),
          cell_text(weight = "bold")
       ),
       locations = cells_body(
          columns = "PC1",
          rows = PC1 >= 1.2 & PC1 < 2 | PC1 <= -1.2
       )) %>%
    tab_style(
       style = list(
          cell_fill(color = "grey"),
          cell_text(weight = "bold")
       ),
       locations = cells_body(
          columns = "PC2",
          rows = PC2 >= .9 & PC2 < 2 | PC2 <= -1
       )) %>%
    tab_style(
       style = list(
          cell_fill(color = "grey"),
          cell_text(weight = "bold")
       ),
       locations = cells_body(
          columns = "PC3",
          rows = PC3 >= 1 & PC3 < 1.15 | PC3 <= -1
       ))
)

### Save ###
setwd(here("outputs/tables"))
gtsave(table, "table_1_pca.png")
