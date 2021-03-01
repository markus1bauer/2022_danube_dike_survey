# Show figure 2 ####
# Markus Bauer
# Citation: Markus Bauer, Jonathan Kiefer & Harald Albrecht (2020) Berichte der Bayerischen Botanischen Gesellschaft 90, 43-66



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Preparation ################################################################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Explanation: #
# edataT = 1984-vs.-2018 dataset,  edataB = 2002-vs.-2018 dataset


### Packages ###
library(tidyverse)
library(adespatial)

### Load data ###
# CONTINUE directly after the script prepare_data_TBI.R #

### Models ###
resTpa <- TBI(vdataTpa[1:40,], vdataTpa[41:80,], method = "sorensen", 
              nperm = 9999, test.t.perm = T, clock = T)
resBpa <- TBI(vdataBpa[1:42,], vdataBpa[43:84,], method = "sorensen",
              nperm = 9999, test.t.perm = T, clock = T)



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# C Plots ################################################################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


par(mfrow = c(1,2))
plot(resTpa, type = "BC", xlim = c(0,0.6), ylim = c(0,0.5), 
     col.bg = "grey50", main = "Artenzahl 1984-2018")
plot(resBpa, type = "BC", xlim = c(0,0.6), ylim = c(0,0.5), 
     col.bg = "grey70", main = "Artenzahl 2003-2018")

