# Calculate TBI ####
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



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Statistics ################################################################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


## 1 Dataset 1984-2018 #####################################################################################

### 1.1 Presence-absence data 1984-2018 --------------------------------------------------------------------------------------------
####Total
resTpa <- TBI(vdataTpa[1:40,], vdataTpa[41:80,], method = "sorensen", 
              nperm = 9999, test.t.perm = T, clock = T)
####Grasses-Herbs-Legumes
resTpaG <- TBI(vdataTpaG[1:40,], vdataTpaG[41:80,], method = "sorensen", 
               nperm = 9999, test.t.perm = T, clock = T)
resTpaH <- TBI(vdataTpaH[1:40,], vdataTpaH[41:80,], method = "sorensen",
               nperm = 9999, test.t.perm = T, clock = T)
resTpaL <- TBI(vdataTpaL[1:40,], vdataTpaL[41:80,], method = "sorensen",
               nperm = 9999, test.t.perm = T, clock = T)
####Ellenberg N=2 vs. N=3
resTpa2 <- TBI(vdataTpa2[1:40,], vdataTpa2[41:80,], method = "sorensen",
               nperm = 9999, test.t.perm = T, clock = T)
resTpa3 <- TBI(vdataTpa3[1:40,], vdataTpa3[41:80,], method = "sorensen", 
               nperm = 9999, test.t.perm = T, clock = T)
####Ellenberg T=5 vs. T=6
resTpa5 <- TBI(vdataTpa5[1:40,], vdataTpa5[41:80,], method = "sorensen", 
               nperm = 9999, test.t.perm = T, clock = T)
resTpa6 <- TBI(vdataTpa6[1:40,], vdataTpa6[41:80,], method = "sorensen", 
               nperm = 9999, test.t.perm = T, clock = T)
####Probeplots
par(mfrow = c(1,1))
plot(resTpa, type = "BC")
resTpa#sig.-
par(mfrow = c(1,3))
plot(resTpaG, type = "BC", xlim = c(0,0.8), ylim = c(0,0.8))
plot(resTpaH, type = "BC", xlim = c(0,0.8), ylim = c(0,0.8))
plot(resTpaL, type = "BC", xlim = c(0,0.8), ylim = c(0,0.8))
resTpaG#n.s.
resTpaH#sig.-
resTpaL#n.s.
par(mfrow = c(1,2))
plot(resTpa2, type = "BC", xlim = c(0,0.8), ylim = c(0,0.8))
plot(resTpa3, type = "BC", xlim = c(0,0.8), ylim = c(0,0.8))
resTpa2
resTpa3
par(mfrow = c(1,2))
plot(resTpa5, type = "BC", xlim = c(0,0.8), ylim = c(0,0.8))
plot(resTpa6, type = "BC", xlim = c(0,0.8), ylim = c(0,0.8))
resTpa5
resTpa6
### 1.2 Abundance data 1984-2018--------------------------------------------------------------------------------------------
####Total
resTabu <- TBI(vdataTabu[1:40,], vdataTabu[41:80,], method = "%diff", 
               nperm = 9999, test.t.perm = T, clock = T)
####Grasses-Herbs-Legumes
resTabuG <- TBI(vdataTabuG[1:40,], vdataTabuG[41:80,], method = "%diff",
                nperm = 9999, test.t.perm = T, clock = T)
resTabuH <- TBI(vdataTabuH[1:40,], vdataTabuH[41:80,], method = "%diff",
                nperm = 9999, test.t.perm = T, clock = T)
resTabuL <- TBI(vdataTabuL[1:40,], vdataTabuL[41:80,], method = "%diff",
                nperm = 9999, test.t.perm = T, clock = T)
####Ellenberg N=2 vs. N=3
resTabu2 <- TBI(vdataTabu2[1:40,], vdataTabu2[41:80,], method = "%diff",
                nperm = 9999, test.t.perm = T, clock = T)
resTabu3 <- TBI(vdataTabu3[1:40,], vdataTabu3[41:80,], method = "%diff",
                nperm = 9999, test.t.perm = T, clock = T)
####Ellenberg T=5 vs. T=6
resTabu5 <- TBI(vdataTabu5[1:40,], vdataTabu5[41:80,], method = "%diff",
                nperm = 9999, test.t.perm = T, clock = T)
resTabu6 <- TBI(vdataTabu6[1:40,], vdataTabu6[41:80,], method = "%diff", 
                nperm = 9999, test.t.perm = T, clock = T)
#### Plotten
par(mfrow = c(1,1))
plot(resTabu, type = "BC")
resTabu#sig.-
par(mfrow = c(1,3))
plot(resTabuG, type = "BC", xlim = c(0,0.8), ylim = c(0,0.8))
plot(resTabuH, type = "BC", xlim = c(0,0.8), ylim = c(0,0.8))
plot(resTabuL, type = "BC", xlim = c(0,0.8), ylim = c(0,0.8))
resTabuG#n.s.
resTabuH#sig.-
resTabuL#sig.-
par(mfrow = c(1,2))
plot(resTabu2, type = "BC", xlim = c(0,0.8), ylim = c(0,0.8))
plot(resTabu3, type = "BC", xlim = c(0,0.8), ylim = c(0,0.8))
resTabu2
resTabu3
par(mfrow = c(1,2))
plot(resTabu5, type = "BC", xlim = c(0,0.8), ylim = c(0,0.8))
plot(resTabu6, type = "BC", xlim = c(0,0.8), ylim = c(0,0.8))
resTabu5
resTabu6

## 2 Dataset 2003-2018 #########################################################################################################

### 2.1 Presence-absence data 2003-2018--------------------------------------------------------------------------------------------
####Total
resBpa <- TBI(vdataBpa[1:42,], vdataBpa[43:84,], method = "sorensen",
              nperm = 9999, test.t.perm = T, clock = T)
####Grasses-Herbs-Legumes
resBpaG <- TBI(vdataBpaG[1:42,], vdataBpaG[43:84,], method = "sorensen", 
               nperm = 9999, test.t.perm = T, clock = T)
resBpaH <- TBI(vdataBpaH[1:42,], vdataBpaH[43:84,], method = "sorensen", 
               nperm = 9999, test.t.perm = T, clock = T)
resBpaL <- TBI(vdataBpaL[1:42,], vdataBpaL[43:84,], method = "sorensen", 
               nperm = 9999, test.t.perm = T, clock = T)
####Ellenberg N=2 vs. N=3
resBpa2 <- TBI(vdataBpa2[1:42,], vdataBpa2[43:84,], method="sorensen",
               nperm = 9999, test.t.perm = T, clock = T)
resBpa3 <- TBI(vdataBpa3[1:42,], vdataBpa3[43:84,], method = "sorensen",
               nperm = 9999, test.t.perm = T, clock = T)
####Ellenberg T=5 vs. T=6
resBpa5 <- TBI(vdataBpa5[1:42,], vdataBpa5[43:84,], method = "sorensen",
               nperm = 9999, test.t.perm = T, clock = T)
resBpa6 <- TBI(vdataBpa6[1:42,], vdataBpa6[43:84,], method = "sorensen",
               nperm = 9999, test.t.perm = T, clock = T)
####North vs. South
resBpaN <- TBI(vdataBpa[1:24,], vdataBpa[43:66,], method = "sorensen", 
               nperm = 9999, test.t.perm = T, clock = T)
resBpaS <- TBI(vdataBpa[25:42,], vdataBpa[67:84,], method = "sorensen",
               nperm = 9999, test.t.perm = T, clock = T)
####Plots
par(mfrow = c(1,1))
plot(resBpa, type = "BC", xlim = c(0,0.7), ylim = c(0,0.3))
resBpa #n.s.
par(mfrow = c(1,3))
plot(resBpaG, type = "BC", xlim = c(0,0.7), ylim = c(0,0.3))
plot(resBpaH, type = "BC", xlim = c(0,0.7), ylim = c(0,0.3))
plot(resBpaL, type = "BC", xlim = c(0,0.7), ylim = c(0,0.3))
resBpaG#sig. +
resBpaH#n.s.
resBpaL#sig. -
par(mfrow = c(1,2))
plot(resBpa2, type = "BC", xlim = c(0,0.7), ylim = c(0,0.3))
plot(resBpa3, type = "BC", xlim = c(0,0.7), ylim = c(0,0.3))
resBpa2
resBpa3
par(mfrow = c(1,2))
plot(resBpa5, type = "BC", xlim = c(0,0.7), ylim = c(0,0.3))
plot(resBpa6, type = "BC", xlim = c(0,0.7), ylim = c(0,0.3))
resBpa5
resBpa6
par(mfrow = c(1,2))
plot(resBpaN, type = "BC", xlim = c(0,0.7), ylim = c(0,0.3))
plot(resBpaS, type = "BC", xlim = c(0,0.7), ylim = c(0,0.3))
resBpaN #n.s.
resBpaS #n.s.
### 2.2 Abundance data 2003-2018--------------------------------------------------------------------------------------------
####Total
resBabu <- TBI(vdataBabu[1:42,], vdataBabu[43:84,], method = "%diff",
               nperm = 9999, test.t.perm = T, clock = T)
####Grasses-Herbs-Legumes
resBabuG <- TBI(vdataBabuG[1:42,], vdataBabuG[43:84,], method = "%diff",
                nperm = 9999, test.t.perm = T, clock = T)
resBabuH <- TBI(vdataBabuH[1:42,], vdataBabuH[43:84,], method = "%diff",
                nperm = 9999, test.t.perm = T, clock = T)
resBabuL <- TBI(vdataBabuL[1:42,], vdataBabuL[43:84,], method = "%diff",
                nperm = 9999, test.t.perm = T, clock = T)
####Ellenberg N=2 vs. N=3
resBabu2 <- TBI(vdataBabu2[1:42,],vdataBabu2[43:84,], method = "%diff",
                nperm = 9999, test.t.perm = T, clock = T)
resBabu3 <- TBI(vdataBabu3[1:42,], vdataBabu3[43:84,], method = "%diff",
                nperm = 9999, test.t.perm = T, clock = T)
####Ellenberg T=5 vs. T=6
resBabu5 <- TBI(vdataBabu5[1:42,], vdataBabu5[43:84,], method = "%diff",
                nperm = 9999, test.t.perm = T, clock = T)
resBabu6 <- TBI(vdataBabu6[1:42,], vdataBabu6[43:84,], method = "%diff",
                nperm = 9999, test.t.perm = T, clock = T)
####North vs. South
resBabuN <- TBI(vdataBabu[1:24,], vdataBabu[43:66,], method = "%diff",
                nperm = 9999, test.t.perm = T, clock = T)
resBabuS <- TBI(vdataBabu[25:42,], vdataBabu[67:84,], method = "%diff",
                nperm = 9999, test.t.perm = T, clock = T)
####Grasses of north vs. grasses of south
resBabuNg <- TBI(vdataBabu[c(1:24),c(2,11:13,18:21,26,45,46)],
                 vdataBabu[c(43:66),c(2,11:13,18:21,26,45,46)], 
                 method = "%diff",
                 nperm = 9999, test.t.perm = T, clock = T)
resBabuSg <- TBI(vdataBabu[c(25:42),c(2,11:13,18:21,26,45,46)],
                 vdataBabu[c(67:84),c(2,11:13,18:21,26,45,46)], 
                 method = "%diff",
                 nperm = 9999, test.t.perm = T, clock = T)
####Plots
par(mfrow = c(1,1))
plot(resBabu, type = "BC", xlim=c(0,0.7), ylim=c(0,0.3))
resBabu#sig.-
par(mfrow = c(1,3))
plot(resBabuG, type = "BC", xlim = c(0,0.7), ylim = c(0,0.3))
plot(resBabuH, type = "BC", xlim = c(0,0.7), ylim = c(0,0.3))
plot(resBabuL, type = "BC", xlim = c(0,0.7), ylim = c(0,0.3))
resBabuG#sig.+
resBabuH#sig.-
resBabuL#n.s.
par(mfrow = c(1,2))
plot(resBabu2, type = "BC", xlim = c(0,0.7), ylim = c(0,0.3))
plot(resBabu3, type = "BC", xlim = c(0,0.7), ylim = c(0,0.3))
resBabu2
resBabu3
par(mfrow = c(1,2))
plot(resBabu5, type = "BC", xlim = c(0,0.7), ylim = c(0,0.3))
plot(resBabu6, type = "BC", xlim = c(0,0.7), ylim = c(0,0.3))
resBabu5
resBabu6
par(mfrow = c(1,2))
plot(resBabuN, type = "BC",xlim = c(0,0.7), ylim = c(0,0.3))
plot(resBabuS, type = "BC", xlim = c(0,0.7), ylim = c(0,0.3))
resBabuN#sig.-
resBabuS#sig.-
par(mfrow = c(1,2))
plot(resBabuNg ,type = "BC", xlim = c(0,0.7), ylim = c(0,0.3))
plot(resBabuSg, type = "BC", xlim = c(0,0.7), ylim = c(0,0.3))
resBabuNg
resBabuSg