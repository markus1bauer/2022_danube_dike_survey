# Calculate TBI ####
# Markus Bauer
# Citation: Bauer M, Huber J, Kollmann J (202x)



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Preparation ################################################################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


### Packages ###
library(here)
library(tidyverse)
library(adespatial)
library(lme4)  
library(DHARMa)
library(emmeans)

### Start ###
rm(list = ls())
setwd(here("data/processed"))

### Load data ###
sites <- read_csv("data_processed_sites.csv", col_names = T, na = c("na", "NA"), col_types = 
                    cols(
                      .default = "?",
                      id = "f",
                      locationAbb = "f",
                      block = "f",
                      plot = "c",
                      locationYear = "f",
                      exposition = "f",
                      side = "f"
                    )) %>%
  filter(vegetationCov > 0) %>%
  select(id, plot, block, locationAbb, surveyYear, constructionYear, locationYear, exposition, side, PC1soil, PC2soil, PC3soil) %>%
  mutate(surveyYearF = as_factor(surveyYear)) %>%
  mutate(constructionYearF = as_factor(constructionYear)) %>%
  add_count(plot) %>%
  filter(n == max(n)) %>%
  select(-n) 

species17 <- read_csv("data_processed_species17.csv", col_names = T, na = "na", col_types = 
                        cols(
                          .default = "d",
                          plot = "f"
                        )) %>%
  column_to_rownames(var = "plot")
species18 <- read_csv("data_processed_species18.csv", col_names = T, na = "na", col_types = 
                        cols(
                          .default = "d",
                          plot = "f"
                        )) %>%
  column_to_rownames(var = "plot")
species19 <- read_csv("data_processed_species19.csv", col_names = T, na = "na", col_types = 
                        cols(
                          .default = "d",
                          plot = "f"
                        )) %>%
  column_to_rownames(var = "plot")
species21 <- read_csv("data_processed_species21.csv", col_names = T, na = "na", col_types = 
                        cols(
                          .default = "d",
                          plot = "f"
                        )) %>%
  column_to_rownames(var = "plot")



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Statistics ################################################################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


## 1 Calculate TBI #####################################################################################

### a 2017 vs. 2018 ---------------------------------------------------------------------------------------------
res1718 <- TBI(species17, species18, method = "%diff", 
               nperm = 9999, test.t.perm = T, clock = T)
res1718$BCD.summary #B = 0.213, C = 0.260, D = 0.473 (45.0% vs. 54.9%)
res1718$t.test_B.C # p.perm = 0.1756
tbi1718 <- as_tibble(res1718$BCD.mat) 
tbi1718 <- sites %>%
  filter(surveyYearF == "2017") %>%
  select(plot, block, surveyYearF, locationYear, exposition, side, PC1soil, PC2soil, PC3soil) %>%
  add_column(tbi1718) %>%
  mutate(plot = factor(plot)) %>%
  mutate(comparison = factor(str_replace(surveyYearF, "2017", "1718"))) %>%
  as_tibble()
colnames(tbi1718) <- c("plot", "block", "surveyYearF", "locationYear", "exposition", "side", "PC1soil", "PC2soil", "PC3soil", "B", "C", "D", "change", "comparison")
#### Test plot
plot(res1718, type = "BC")

### b 2018 vs. 2019 ---------------------------------------------------------------------------------------------
res1819 <- TBI(species18, species19, method = "%diff", 
               nperm = 9999, test.t.perm = T, clock = T)
res1819$BCD.summary #B = 0.167, C = 0.302, D = 0.470 (35.7% vs. 64.2%)
res1819$t.test_B.C # p.perm = 1e-04
tbi1819 <- as_tibble(res1819$BCD.mat) 
tbi1819 <- sites %>%
  filter(surveyYearF == "2017") %>%
  select(plot, block, surveyYearF, locationYear, exposition, side, PC1soil, PC2soil, PC3soil) %>%
  add_column(tbi1819) %>%
  mutate(plot = factor(plot)) %>%
  mutate(comparison = factor(str_replace(surveyYearF, "2017", "1819"))) %>%
  as_tibble()
colnames(tbi1819) <- c("plot", "block", "surveyYearF", "locationYear", "exposition", "side", "PC1soil", "PC2soil", "PC3soil", "B", "C", "D", "change", "comparison")
#### Test plot
plot(res1819, type = "BC")

### c 2019 vs. 2021 ---------------------------------------------------------------------------------------------
res1921 <- TBI(species19, species21, method = "%diff", 
               nperm = 9999, test.t.perm = T, clock = T)
res1921$BCD.summary #B = 0.331, C = 0.168, D = 0.499 (66.3% vs. 33.6%)
res1921$t.test_B.C # p.perm = 1e-04
tbi1921 <- as_tibble(res1921$BCD.mat) 
tbi1921 <- sites %>%
  filter(surveyYearF == "2017") %>%
  select(plot, block, surveyYearF, locationYear, exposition, side, PC1soil, PC2soil, PC3soil) %>%
  add_column(tbi1921) %>%
  mutate(plot = factor(plot)) %>%
  mutate(comparison = factor(str_replace(surveyYearF, "2017", "1921"))) %>%
  as_tibble()
colnames(tbi1921) <- c("plot", "block", "surveyYearF", "locationYear", "exposition", "side", "PC1soil", "PC2soil", "PC3soil", "B", "C", "D", "change", "comparison")
#### Test plot
plot(res1921, type = "BC")

### d Combine datasets ---------------------------------------------------------------------------------------------
data <- bind_rows(tbi1718, tbi1819, tbi1921) %>%
  select(-change) %>%
  pivot_longer(-c(plot:PC3soil, comparison), names_to = "index", values_to = "tbi") %>%
  mutate(index = factor(index),
         exposition = factor(exposition),
         block = factor(block),
         locationYear = factor(locationYear)) %>%
  filter(exposition == "south" | exposition == "north")


## 2 Modelling #############################################################################################


### 1 Data exploration #####################################################################################

#### a Graphs ---------------------------------------------------------------------------------------------
#2way
ggplot(data, aes(x = comparison, y = tbi)) + 
  geom_boxplot() +
  geom_quasirandom() + 
  facet_wrap(~index)
ggplot(data, aes(x = exposition, y = tbi)) + 
  geom_boxplot() +
  geom_quasirandom() + 
  facet_wrap(~index)
ggplot(data, aes(x = side, y = tbi)) + 
  geom_boxplot() +
  geom_quasirandom() + 
  facet_wrap(~index)
ggplot(data, aes(x = PC1soil, y = tbi)) + 
  geom_smooth(method = "lm") +
  geom_point() + 
  facet_wrap(~index)
ggplot(data, aes(x = exp(PC2soil), y = tbi)) + 
  geom_smooth(method = "lm") +
  geom_point() + 
  facet_wrap(~index)
ggplot(data, aes(x = PC3soil, y = tbi)) + 
  geom_smooth(method = "lm") +
  geom_point() + 
  facet_wrap(~index)
#3way
ggplot(data, aes(x = exposition, y = tbi, color = comparison)) + 
  geom_boxplot() +
  facet_wrap(~index)
ggplot(data, aes(x = comparison, y = tbi, color = exposition)) + 
  geom_boxplot() +
  facet_wrap(~index)
ggplot(data, aes(x = PC2soil, y = tbi, color = comparison)) + 
  geom_smooth(method = "lm") +
  geom_point() +
  facet_wrap(~index)

#### b Outliers, zero-inflation, transformations? -----------------------------------------------------
dotchart((data$tbi), groups = factor(data$exposition), main = "Cleveland dotplot")
data %>% count(locationYear)
boxplot(data$tbi);#identify(rep(1, length(edata$rgr13)), edata$rgr13, labels = c(edata$n))
plot(table((data$tbi)), type = "h", xlab = "Observed values", ylab = "Frequency")
ggplot(data, aes(tbi)) + geom_density()
ggplot(data, aes(sqrt(tbi))) + geom_density()
ggplot(data, aes(PC2soil)) + geom_density()
ggplot(data, aes(sqrt(PC2soil))) + geom_density()


### 2 Model building ################################################################################

#### a models ----------------------------------------------------------------------------------------
#random structure
m1a <- lmer(tbi ~ 1 + (1|plot) + (1|locationYear), data, REML = F)
VarCorr(m1a) 
m1b <- lmer(tbi ~ 1 + (1|plot), data, REML = F)
VarCorr(m1b)
m1c <- lmer(tbi ~ 1 + (1 + exposition|plot), data, REML = F)
VarCorr(m1c);isSingular(m1c)
#fixed effects
m2 <- lmer(sqrt(tbi) ~ index * comparison * (exposition + PC1soil) + PC2soil + side + PC3soil + locationYear +
             (1|plot), 
           data = data)
isSingular(m2);simulateResiduals(m2, plot = T)
m3 <- lmer(sqrt(tbi) ~ index * comparison * (exposition + PC2soil) + PC1soil + side + PC3soil + locationYear +
             (1|plot), 
           data = data)
isSingular(m3);simulateResiduals(m3, plot = T)
m4 <- lmer(sqrt(tbi) ~ index * comparison + exposition + PC2soil + PC1soil + side + PC3soil + locationYear +
             (1|plot), 
           data = data)
isSingular(m4);simulateResiduals(m4, plot = T)
m5 <- lmer(sqrt(tbi) ~ index + comparison + exposition + PC1soil + PC2soil + side + PC3soil + locationYear +
             (1 + exposition|plot), 
           data = data)
isSingular(m5);simulateResiduals(m5, plot = T)
m6 <- lmer(sqrt(tbi) ~ index * comparison * exposition + PC1soil + PC2soil + side + PC3soil + locationYear +
             (1|plot), 
           data = data)
isSingular(m6);simulateResiduals(m6, plot = T)



#### b comparison -----------------------------------------------------------------------------------------
anova(m2, m3, m4, m5, m6)
rm(m1a, m1b, m1c, m2, m3, m4, m5)

#### c model check -----------------------------------------------------------------------------------------
simulationOutput <- simulateResiduals(m6, plot = F)
plotResiduals(simulationOutput$scaledResiduals, data$locationYear)
plotResiduals(simulationOutput$scaledResiduals, data$block)
plotResiduals(simulationOutput$scaledResiduals, data$plot)
plotResiduals(simulationOutput$scaledResiduals, data$index)
plotResiduals(simulationOutput$scaledResiduals, data$comparison)
plotResiduals(simulationOutput$scaledResiduals, data$exposition)
plotResiduals(simulationOutput$scaledResiduals, data$side)
plotResiduals(simulationOutput$scaledResiduals, data$PC1soil)
plotResiduals(simulationOutput$scaledResiduals, data$PC2soil)
plotResiduals(simulationOutput$scaledResiduals, data$PC3soil)


### 3 Chosen model output ################################################################################

### * Model output ####
MuMIn::r.squaredGLMM(m6) #R2m = 0.470, R2c = 0.603
VarCorr(m6)
sjPlot::plot_model(m6, type = "re", show.values = T)
car::Anova(m6, type = 3)

### * Effect sizes ####
(emm <- emmeans(m6, revpairwise ~ comparison | index | exposition, type = "response"))
plot(emm, comparison = T)

### * Save ####
table <- broom::tidy(car::Anova(m6, type = 3))
write.csv(table, here("outputs/statistics/table_anova_tbi_abundance_year.csv"))

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# C Plotten ################################################################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

library(ggbeeswarm)
library(ggeffects)

themeMB <- function(){
  theme(
    panel.background = element_rect(fill = "white"),
    text  = element_text(size = 9, color = "black"),
    axis.text  = element_text(size = 9, color = "black"),
    legend.text  = element_text(size = 8, color = "black"),
    strip.text = element_text(size = 9, color = "black"),
    axis.text.y = element_text(angle = 0, hjust = 0.5),
    axis.line.y = element_line(),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.key = element_rect(fill = "transparent"),
    legend.justification = c("right", "top"),
    legend.direction = "vertical",
    legend.position = c(.95, .17),
    legend.background = element_rect(fill = "transparent"),
    legend.margin = margin(0, 0, 0, 0, "cm"),
    plot.margin = margin(0, 0, 0, 0, "cm")
  )
}

pdata <- ggemmeans(m6, terms = c("comparison", "index","exposition"), type = "fe") %>%
  mutate(group = fct_recode(group, Losses = "B", Gains = "C", Total = "D")) %>%
  mutate(x = fct_recode(x, "'17 vs '18" = "1718", "'18 vs '19" = "1819", "'19 vs '21" = "1921")) %>%
  mutate(facet = fct_recode(facet, South = "south", North = "north"))
pd <- position_dodge(.1)
(graph <- ggplot(pdata, 
                 aes(x, predicted, ymin = conf.low, ymax = conf.high, shape = group)) +
    #geom_quasirandom(data = sites, aes(x, predicted), 
    #                 color = "grey70", dodge.width = .6, size = 0.7) +
    geom_errorbar(position = pd, width = 0.0, size = 0.4) +
    geom_point(position = pd, size = 2.5, fill = "white") +
    facet_grid(~ facet) +
    #annotate("text", label = "n.s.", x = 2.2, y = 0.6) +
    scale_y_continuous(limits = c(0, 0.63), breaks = seq(-100, 100, 0.1)) +
    scale_shape_manual(values = c(25, 24, 16)) +
    labs(x = "", y = "Temporal beta diversity", shape = "") +
    themeMB() +
    theme()
)

### Save ###
ggsave(here("outputs/figures/figure_4_tbi_abundance_year_(800dpi_9x5cm).tiff"),
       dpi = 800, width = 9, height = 5, units = "cm")
