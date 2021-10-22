### Show table PERMDISP ###
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
setwd(here("data/tables"))

### Load data ###
data <- read_csv2("table_anova_permdisp.csv", col_names = T, na = "na", col_types = 
                     cols(
                       .default = col_guess()
                       )) %>%
   mutate(across(where(is.numeric), round, 3)) %>%
   filter(term != "(Intercept)")



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
       columns = "term"
     )
   ) %>%
   tab_style(
     style = cell_text(
       align = "left"
     ),
     locations = cells_body(
       columns = "term"
     )
   ) %>%
#heading
   tab_header( 
     title = "Table 3"
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
     df = md("**df**")
   ) %>%
#footnote
   tab_source_note( 
     source_note = "ANOVA of permdisp for old dikes"
   ) %>%
#Highlight certain cells   
   tab_style(
      style = list(
         cell_fill(color = "grey"),
         cell_text(weight = "bold")
      ),
      locations = cells_body(
         columns = "p.value",
         rows = p.value < .05
      ))
)


### Save ###
setwd(here("outputs/tables"))
gtsave(table, "table_anova_permdisp.png")
