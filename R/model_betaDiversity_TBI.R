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
                      side = "f",
                      ffh = "f",
                      vegetationCov = "d",
                      locationYear = "f"
                    )) %>%
  select(id, plot, block, locationAbb, surveyYear, constructionYear, locationYear, exposition, side, PC1, PC2, PC3) %>%
  mutate(surveyYearF = as_factor(surveyYear)) %>%
  mutate(constructionYearF = as_factor(constructionYear)) %>%
  add_count(plot) %>%
  filter(n == max(n)) %>%
  select(-n) 

species <- read_csv("data_processed_species.csv", col_names = T, na = "na", col_types = 
                      cols(
                        .default = "d",
                        name = "f"
                      )) %>%  
  pivot_longer(-name, "id", "value") %>%
  pivot_wider(id, name) %>%
  arrange(id) %>%
  semi_join(sites, by = "id") %>%
  column_to_rownames("id") %>%
  select(which(!colSums(., na.rm = T) %in% 0)) %>%
  rownames_to_column(var = "id") %>%
  mutate(year = factor(str_match(id, "\\d\\d\\d\\d"))) %>%
  mutate(plot = factor(str_match(id, "\\d\\d"))) %>%
  select(-id)

species17 <- species %>%
  filter(year == 2017) %>%
  column_to_rownames("plot") %>%
  select(-year)
species18 <- species %>%
  filter(year == 2018) %>%
  column_to_rownames("plot") %>%
  select(-year)
species19 <- species %>%
  filter(year == 2019) %>%
  column_to_rownames("plot") %>%
  select(-year)



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Statistics ################################################################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


## 1. 2017 vs. 2019 #####################################################################################

res1719 <- TBI(species17, species19, method = "%diff", 
               nperm = 9999, test.t.perm = T, clock = T)
res1719$BCD.summary #B = 0.216, C = 0.375, D = 0.591 (36% vs. 63%)
res1719$t.test_B.C # p = 1e-04
tbi1719 <- as_tibble(res1719$BCD.mat) 
tbi1719 <- sites %>%
  filter(surveyYearF == "2018") %>%
  select(plot, block, surveyYearF, locationYear, exposition, side, PC1, PC2, PC3) %>%
  add_column(tbi1719) %>%
  mutate(plot = factor(plot)) %>%
  mutate(comparison = factor(str_replace(surveyYearF, "2018", "1719"))) %>%
  as_tibble()
colnames(tbi1719) <- c("plot", "block", "surveyYearF", "locationYear", "exposition", "side", "PC1", "PC2", "PC3", "B", "C", "D", "change", "comparison")
#### Test plot
plot(res1719, type = "BC")


## 2. 2017 vs. 2018 #####################################################################################

res1718 <- TBI(species17, species18, method = "%diff", 
               nperm = 9999, test.t.perm = T, clock = T)
res1718$BCD.summary #B = 0.215, C = 0.255, D = 0.470 (46% vs. 54%)
res1718$t.test_B.C # p = 1.1e-01
tbi1718 <- as_tibble(res1718$BCD.mat) 
tbi1718 <- sites %>%
  filter(surveyYearF == "2017") %>%
  select(plot, block, surveyYearF, locationYear, exposition, side, PC1, PC2, PC3) %>%
  add_column(tbi1718) %>%
  mutate(plot = factor(plot)) %>%
  mutate(comparison = factor(str_replace(surveyYearF, "2017", "1718"))) %>%
  as_tibble()
colnames(tbi1718) <- c("plot", "block", "surveyYearF", "locationYear", "exposition", "side", "PC1", "PC2", "PC3", "B", "C", "D", "change", "comparison")
#### Test plot
plot(res1718, type = "BC")

## 3. 2018 vs. 2019 #####################################################################################

res1819 <- TBI(species18, species19, method = "%diff", 
              nperm = 9999, test.t.perm = T, clock = T)
res1819$BCD.summary #B = 0.175, C = 0.296, D = 0.471 (37% vs. 63%)
res1819$t.test_B.C # p = 1e-04
tbi1819 <- as_tibble(res1819$BCD.mat)
tbi1819 <- sites %>%
  filter(surveyYearF == "2019") %>%
  select(plot, block, surveyYearF, locationYear, exposition, side, PC1, PC2, PC3) %>%
  add_column(tbi1819) %>%
  mutate(plot = factor(plot)) %>%
  mutate(comparison = factor(str_replace(surveyYearF, "2019", "1819"))) %>%
  as_tibble()
colnames(tbi1819) <- c("plot", "block", "surveyYearF", "locationYear", "exposition", "side", "PC1", "PC2", "PC3", "B", "C", "D", "change", "comparison")
#### Test plot
plot(res1819, type = "BC")


## 3. Combined data sets #####################################################################################

tbi <- bind_rows(tbi1719, tbi1718, tbi1819) %>%
  select(-change) %>%
  pivot_longer(-c(plot:PC3, comparison), names_to = "index", values_to = "tbi") %>%
  mutate(index = factor(index)) %>%
  filter(comparison != 1719) %>%
  filter(exposition == "south" | exposition == "north")
  

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# C Plotten ################################################################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

library(lme4)  
library(DHARMa)
library(emmeans)
library(ggbeeswarm)
library(ggeffects)
m1 <- lmer(tbi ~ index * comparison * (exposition + side + PC1 + PC2 + PC3) +
             index:exposition:comparison + index:side:comparison +
            (1 + exposition|locationYear), data = tbi)
simulationOutput <- simulateResiduals(m1, plot = T)
MuMIn::r.squaredGLMM(m1)
car::Anova(m1, type = 3)
#rename B C D to declines increases and TBI
#data2 <- rename(data, predicted = tbi, x = exposition, group = surveyYearF)
ggplot(tbi, aes(y = tbi, x = PC2, colour = comparison)) +
  geom_point() +
  geom_smooth(method = "lm", se = T) +
  facet_grid(~index) +
  themeMB()

setwd(here("outputs/figures"))
ggsave("figure_4_tbi_PC2_(800dpi_16x10cm).tiff",
       dpi = 800, width = 16, height = 10, units = "cm")

ggplot(tbi, aes(y = tbi, x = exposition)) +
  geom_quasirandom() +
  geom_boxplot() +
  facet_grid(~index) +
  themeMB()

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
    legend.position = c(.95, .45),
    legend.background = element_rect(fill = "transparent"),
    legend.margin = margin(0, 0, 0, 0, "cm"),
    plot.margin = margin(0, 0, 0, 0, "cm")
  )
}

pdata <- ggemmeans(m1, terms = c("exposition", "comparison", "index"), type = "fe") %>%
  mutate(facet = fct_recode(facet, Declines = "B", Increases = "C", Total = "D")) %>%
  mutate(group = fct_recode(group, "2017/normal vs. 2018/dry" = "1718", "2018/dry vs. 2019/dry" = "1819")) %>%
  mutate(x = fct_recode(x, South = "south", North = "north"))
pd <- position_dodge(.6)
(graph <- ggplot(pdata, 
                 aes(x, predicted, ymin = conf.low, ymax = conf.high, shape = group)) +
    #geom_quasirandom(data = data2, aes(x, predicted), 
     #                color = "grey70", dodge.width = .6, size = 0.7) +
    geom_errorbar(position = pd, width = 0.0, size = 0.4) +
    geom_point(position = pd, size = 2.5, fill = "white") +
    facet_grid(~ facet) +
    annotate("text", label = "n.s.", x = 2.2, y = 0.6) +
    scale_y_continuous(limits = c(0, 0.6), breaks = seq(-100, 100, 0.1)) +
    scale_shape_manual(values = c(21, 16)) +
    labs(x = "Exposition", y = "Temporal beta-diversity", shape = "") +
    themeMB() +
    theme()
)

### Save ###
setwd(here("outputs/figures"))
ggsave("figure_4_tbi_(800dpi_9x5cm).tiff",
       dpi = 800, width = 9, height = 5, units = "cm")


### 1 Data exploration #####################################################################################

#### a Graphs ---------------------------------------------------------------------------------------------
#simple effects:
ggplot(sites, aes(y = n, x = locationAbb)) + geom_boxplot() + geom_quasirandom()
ggplot(sites, aes(y = n, x = constructionYear)) + geom_point() + geom_smooth(method = "lm") 
ggplot(sites, aes(y = n, x = constructionYearF)) + geom_boxplot() + geom_quasirandom()
ggplot(sites, aes(y = n, x = plotAge)) + geom_point() + geom_smooth(method = "lm") 
ggplot(sites, aes(y = n, x = surveyYearF)) + geom_boxplot() + geom_quasirandom()
ggplot(sites, aes(y = n, x = exposition)) + geom_boxplot() + geom_quasirandom()
ggplot(sites, aes(y = n, x = side)) + geom_boxplot() + geom_quasirandom()
ggplot(sites, aes(y = n, x = PC1)) + geom_point() + geom_smooth(method = "lm")
ggplot(sites, aes(y = n, x = PC2)) + geom_point() + geom_smooth(method = "lm")
ggplot(sites, aes(y = n, x = PC3)) + geom_point() + geom_smooth(method = "lm")
ggplot(sites, aes(y = n, x = changeType)) + geom_boxplot() + geom_quasirandom()
ggplot(sites, aes(y = n, x = ffh)) + geom_boxplot() + geom_quasirandom()
ggplot(sites, aes(y = n, x = baykompv)) + geom_boxplot() + geom_quasirandom()
ggplot(sites, aes(y = n, x = biotopeType)) + geom_boxplot() + geom_quasirandom()
#2way
ggplot(sites, aes(x = locationAbb, y = n, color = exposition)) + 
  geom_boxplot() +
  geom_quasirandom()
ggplot(sites, aes(x = locationAbb, y = n, color = exposition)) + 
  geom_boxplot() +
  facet_wrap(~exposition)
ggplot(sites, aes(x = locationAbb, y = n, color = surveyYearF)) + 
  geom_boxplot()
ggplot(sites, aes(x = PC2, y = n, color = surveyYearF)) + 
  geom_point() + 
  geom_smooth(method = "lm") + 
  facet_wrap(~surveyYearF)
ggplot(sites, aes(x = plotAge, y = n, color = changeType)) + 
  geom_smooth(aes(group = plot), colour = "grey40", method = "lm", se = F, show.legend = F) + 
  geom_point() + 
  geom_smooth(method = "lm") + 
  facet_wrap(~changeType)


##### b Outliers, zero-inflation, transformations? -----------------------------------------------------
dotchart((sites$n), groups = factor(sites$exposition), main = "Cleveland dotplot")
sites %>% count(locationAbb)
boxplot(sites$n);#identify(rep(1, length(edata$rgr13)), edata$rgr13, labels = c(edata$n))
plot(table((sites$n)), type = "h", xlab = "Observed values", ylab = "Frequency")
ggplot(sites, aes(n)) + geom_density()
ggplot(sites, aes(sqrt(n))) + geom_density()


## 2 Model building ################################################################################

#### a models ----------------------------------------------------------------------------------------
#random structure
m1a <- lmer(n ~ 1 + (surveyYear|locationAbb/plot), sites, REML = F)
VarCorr(m1a) # convergence problems
m1b <- lmer(n ~ 1 + (surveyYear|plot), sites, REML = F)
VarCorr(m1b) # convergence problems
m1c <- lmer(n ~ 1 + (1|locationAbb/plot), sites, REML = F)
VarCorr(m1c)
m1d <- lmer(n ~ 1 + (1|plot), sites, REML = F)
VarCorr(m1d)
#fixed effects
m2 <- lmer((n) ~ (exposition + side + PC1 + PC2 + PC3) + 
             (1|plot), sites, REML = F)
simulateResiduals(m2, plot = T)
#fixed and site and year effects
m3 <- lmer((n) ~ (exposition + side + PC1 + PC2 + PC3 + surveyYearF + locationAbb) +
             (1|plot), sites, REML = F);
simulateResiduals(m3, plot = T)
isSingular(m3)
#plotAge instead of location
m4 <- lmer((n) ~ (exposition + side + PC1 + PC2 + PC3 + surveyYearF + plotAge) + 
             (1|plot), sites, REML = F);
simulateResiduals(m4, plot = T)
isSingular(m4)

m5 <- lmer((n) ~ (exposition + side + PC1 + PC2 + PC3 + surveyYearF + locationAbb) +
             locationAbb:exposition + locationAbb:surveyYearF +
             (1|plot), sites, REML = F);
simulateResiduals(m5, plot = T)
isSingular(m5)


#### b comparison -----------------------------------------------------------------------------------------
anova(m2, m3, m4)
rm(m1a, m1b, m1c, m1d, m4)

#### c model check -----------------------------------------------------------------------------------------
simulationOutput <- simulateResiduals(m3, plot = F)
plotResiduals(simulationOutput$scaledResiduals, sites$surveyYearF)
plotResiduals(simulationOutput$scaledResiduals, sites$locationAbb)
plotResiduals(simulationOutput$scaledResiduals, sites$block)
plotResiduals(simulationOutput$scaledResiduals, sites$plot)
plotResiduals(simulationOutput$scaledResiduals, sites$side)
plotResiduals(simulationOutput$scaledResiduals, sites$exposition)
plotResiduals(simulationOutput$scaledResiduals, sites$PC1)
plotResiduals(simulationOutput$scaledResiduals, sites$PC2)
plotResiduals(simulationOutput$scaledResiduals, sites$PC3)


## 3 Chosen model output ################################################################################

### Model output ---------------------------------------------------------------------------------------------
MuMIn::r.squaredGLMM(m1)
0.384 /  0.500
VarCorr(m1)
sjPlot::plot_model(m1, type = "re", show.values = T)
car::Anova(m1, type = 2)

### Effect sizes -----------------------------------------------------------------------------------------
(emm <- emmeans(m1, revpairwise ~ exposition * comparison | index, type = "response"))
plot(emm, comparison = T)
(emm <- emmeans(m1, revpairwise ~ side, type = "response"))

### Save ###
table <- broom::tidy(car::Anova(m3, type = 2))
setwd(here("data/tables"))
write.csv2(table, "table_anova_cwmAbuSla.csv")
