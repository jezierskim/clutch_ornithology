#### DATA ANALYSIS CODE
#### Manuscript 'Birds that nest exclusively on islands have smaller clutches'. By M.T.Jezierski.
#### Versioned 17th December 2023.

#### Required packages.
library(tidyverse)
library(ape)
library(phylolm)
library(nlme)
library(car)
library(MuMIn)
library(caper)
library(ggdist)
library(gghalves)
library(cowplot)

#### SET WORKING DIRECTORY
setwd('path/to/this/directory')

##### LOAD IN THE DATA REQUIRED: ====
## Load clutch dataset; create clutch-midpoint variable; remove outliers for all three values.
clutch <- read.csv('analysis_dataset.csv') # load in the clutch dataset.
clutch <- clutch %>% mutate(migratory = ifelse(is.na(clutch$migratory), 'resident', clutch$migratory))
seabirds <- read_csv('seabirds.csv') %>% unite(binomial, genus:species, sep = ' ') # seabird dataset.
clutch <- clutch %>% mutate(seabird = ifelse(binomial %in% seabirds$binomial, 'seabird', 'landbird')) %>% drop_na() # add it to the analysis dataset.

##### REMOVE OUTLIERS AND PREPARE THE SEA AND COMPLETE DATASETS: ====
Q <- quantile(clutch$clu_geom, probs=c(.25, .75), na.rm = FALSE) # calculate 25% and 75% quantiles.
IQR <- IQR(clutch$clu_geom) # calculate IQR
clutch <- clutch %>% # remove outliers from both ends of the distribution.
  filter(clu_geom < (Q[2]+1.5*IQR)) %>%
  filter(clu_geom > (Q[1] -1.5*IQR)) # END RESULT: dataframe of 4,530 species to be retained in the analysis.
rm(IQR, Q) # not needed anymore.

## Create landbird only, landbird and seabird and seabird only dataset (the latter two used at the end).
clutch_land_and_sea <- clutch # complete dataset.
clutch_sea <- clutch %>% filter(seabird == 'seabird') # seabird only dataset.
clutch_land <- clutch %>% filter(seabird == 'landbird') # landbird only dataset.

##### LOAD PHYLOGENETIC TREES: ====
trees <- read.tree('100_trees_Hackett_genetic.tre') # random sample of 100 trees pruned to the dataset.

##### MODELLING PREPARATION -- LANDBIRD-ONLY DATASET ====
### Firstly, prepare variables: calculate log of area; calculate absolute latitude, check VIF
### run different evolutionary models to see which one fits the data best; run a PGLS model with the best evolutionary model
### for each metric of clutch size to examine model diagnostics.

## Prepare a dataset to include log of area and absolute latitude. 
clutch_land <- clutch_land %>% mutate( # log maximum body mass 
                            area = log10(area), # log area 
                            latitude = abs(latitude)) %>% # force latitude to absolute
  column_to_rownames(var = 'treelabels') %>% drop_na() # move phylogeny names to rownames for downstream use with phylolm.

## Check VIF using basic lm() in R.
vif(lm(clu_geom ~ is*latitude + migratory + mode + mass_geom, data = clutch_land), type = 'predictor')

## Examine the best evolutionary model using different phylolm() run of the model.
evo_fit <- function(x) {model.sel(A <- phylolm(formula = x, data = clutch_land, as.phylo(trees[[1]]), model = 'BM'),
                       B <- phylolm(formula = x, data = clutch_land, phy = as.phylo(trees[[1]]), model = 'OUrandomRoot'),
                       C <- phylolm(formula = x, data = clutch_land, as.phylo(trees[[1]]), model = 'OUfixedRoot'),
                       D <- phylolm(formula = x, data = clutch_land, as.phylo(trees[[1]]), model = 'lambda'),
                       E <- phylolm(formula = x, data = clutch_land, as.phylo(trees[[1]]), model = 'kappa'),
                       G <- phylolm(formula = x, data = clutch_land, as.phylo(trees[[1]]), model = 'delta'))
} # function that test 6 different evolutionary models (see METHODS)
evo_fit(clu_geom ~ is*latitude + migratory + mode + mass_geom) # for clumin lambda is the best; warning of 534 missing taxa is correct; filtered out relative
# to the tree.

##### MODEL EXAMINATION -- LANDBIRD-ONLY DATASET ====
### Check the health of the model and the linear model assumptions. 
### WARNING: the code run on nlme takes significant time (possibly overnigh on a Macbook Air 2020).
cd <- comparative.data(clutch_land %>% rownames_to_column(var = 'treelabels'), phy = trees[[1]], names.col = treelabels) # created phylogeny-ordered data.frame
nlmePagelclutch <- gls(clu_geom ~ is*latitude + migratory + mode + mass_geom, # run gls with Pagel; NB: THIS STEP TAKES A WHILE.
                       correlation = corPagel(1, phy = cd$phy),
                       data = cd$data, method = "ML")
plot(fitted(nlmePagelclutch), resid(nlmePagelclutch, type = 'normalized')) # plots normalised residuals against fitted. 
lines(predict(lm(fitted(nlmePagelclutch)~resid(nlmePagelclutch, type = 'normalized'))),col='green') # prediction line to show the trend line.
abline(h = 0, col = "darkorange", lwd = 2) # 0 line to compare
qqnorm(resid(nlmePagelclutch)) # distribution of quantiles 
qqline(resid(nlmePagelclutch)) # quantile line for the distribution.

## NB: From now on code is much less intensive and all functions can run in < 1.5h.

##### MODELLING -- LANDBIRD-ONLY DATASET -- RUN ACROSS 100 TREES. ======
clutch_land_estimates <- data.frame()
for (i in 1:length(trees)) {
  m <- phylolm(clu_geom ~ is*latitude + migratory + mode + mass_geom, data = clutch_land, phy = trees[[i]], model = 'lambda')
  d <- as.data.frame(summary(m)$coefficients) %>% rownames_to_column(var = 'predictor') %>%
    mutate(lambda = m$optpar) %>% # extract Pagel's lambda for each model
    mutate(r2 = m$adj.r.squared) # extract R^squared for each model
  clutch_land_estimates <- rbind(clutch_land_estimates, d)
}
write_csv(clutch_land_estimates, 'clutch_land_estimates.csv') # save the output of the model for further checks.

##### RESULTS -- GENERAL STATISTICS AND MODEL RESULTS ON LANDBIRD-ONLY DATASET ====
clutch_land %>% group_by(is) %>% summarise(mean = mean(clu_geom)) # mean geometric clutch size between continental/island-breeding endemics. 

## Median of beta coefficients - landbird-only dataset.
round(median(clutch_land_estimates[clutch_land_estimates$predictor == 'isY', 2]),3)
round(median(clutch_land_estimates[clutch_land_estimates$predictor == 'latitude', 2]),3) # NB: result significant for FIGURE 2.
round(median(clutch_land_estimates[clutch_land_estimates$predictor == 'modeprecocial', 2]),3)
round(median(clutch_land_estimates[clutch_land_estimates$predictor == 'migratoryresident', 2]),3)
round(median(clutch_land_estimates[clutch_land_estimates$predictor == 'mass_geom', 2]), 3)
round(median(clutch_land_estimates[clutch_land_estimates$predictor == 'isY:latitude', 2]), 3) # NB: result significant for FIGURE 2.
### Quantiles of beta coefficients - landbird-only dataset.
round(quantile(clutch_land_estimates[clutch_land_estimates$predictor == 'isY', 2], c(0.05, 0.95)),3)
round(quantile(clutch_land_estimates[clutch_land_estimates$predictor == 'latitude', 2], c(0.05, 0.95)),3) # NB: result significant for FIGURE 2.
round(quantile(clutch_land_estimates[clutch_land_estimates$predictor == 'modeprecocial', 2], c(0.05, 0.95)),3)
round(quantile(clutch_land_estimates[clutch_land_estimates$predictor == 'migratoryresident', 2], c(0.05, 0.95)),3)
round(quantile(clutch_land_estimates[clutch_land_estimates$predictor == 'mass_geom', 2], c(0.05, 0.95)),3)
round(quantile(clutch_land_estimates[clutch_land_estimates$predictor == 'isY:latitude', 2], c(0.05, 0.95)) ,3)# NB: result significant for FIGURE 2.
### Medians of standard errors - landbird-only dataset
round(median(clutch_land_estimates[clutch_land_estimates$predictor == 'isY', 3]),3)
round(median(clutch_land_estimates[clutch_land_estimates$predictor == 'latitude', 3]),3) # NB: result significant for FIGURE 2.
round(median(clutch_land_estimates[clutch_land_estimates$predictor == 'modeprecocial', 3]),3)
round(median(clutch_land_estimates[clutch_land_estimates$predictor == 'migratoryresident', 3]),3)
round(median(clutch_land_estimates[clutch_land_estimates$predictor == 'mass_geom', 3]), 3)
round(median(clutch_land_estimates[clutch_land_estimates$predictor == 'isY:latitude', 3]), 3) # NB: result significant for FIGURE 2.
### Quantiles of standard errors - landbird-only dataset.
round(quantile(clutch_land_estimates[clutch_land_estimates$predictor == 'isY', 3], probs = c(0.05, 0.95)),3)
round(quantile(clutch_land_estimates[clutch_land_estimates$predictor == 'latitude', 3], c(0.05, 0.95)),3) # NB: result significant for FIGURE 2.
round(quantile(clutch_land_estimates[clutch_land_estimates$predictor == 'modeprecocial', 3], c(0.05, 0.95)),3)
round(quantile(clutch_land_estimates[clutch_land_estimates$predictor == 'migratoryresident', 3], c(0.05, 0.95)),3)
round(quantile(clutch_land_estimates[clutch_land_estimates$predictor == 'mass_geom', 3], c(0.05, 0.95)),3)
round(quantile(clutch_land_estimates[clutch_land_estimates$predictor == 'isY:latitude', 3], c(0.05, 0.95)),3) # NB: result significant for FIGURE 2.
### Medians of p-values - landbird-only dataset.
round(median(clutch_land_estimates[clutch_land_estimates$predictor == 'isY', 5]),3)
round(median(clutch_land_estimates[clutch_land_estimates$predictor == 'latitude', 5]),3) # NB: result significant for FIGURE 2.
round(median(clutch_land_estimates[clutch_land_estimates$predictor == 'modeprecocial', 5]),3)
round(median(clutch_land_estimates[clutch_land_estimates$predictor == 'migratoryresident', 5]),3)
round(median(clutch_land_estimates[clutch_land_estimates$predictor == 'mass_geom', 5]), 3)
round(median(clutch_land_estimates[clutch_land_estimates$predictor == 'isY:latitude', 5]), 3) # NB: result significant for FIGURE 2.
### Quantiles of p-values - landbird-only dataset.
round(quantile(clutch_land_estimates[clutch_land_estimates$predictor == 'isY', 5], probs = c(0.05, 0.95)),3)
round(quantile(clutch_land_estimates[clutch_land_estimates$predictor == 'latitude', 5], c(0.05, 0.95)),3) # NB: result significant for FIGURE 2.
round(quantile(clutch_land_estimates[clutch_land_estimates$predictor == 'modeprecocial', 5], c(0.05, 0.95)),3)
round(quantile(clutch_land_estimates[clutch_land_estimates$predictor == 'migratoryresident', 5], c(0.05, 0.95)),3)
round(quantile(clutch_land_estimates[clutch_land_estimates$predictor == 'mass_geom', 5], c(0.05, 0.95)),3)
round(quantile(clutch_land_estimates[clutch_land_estimates$predictor == 'isY:latitude', 5], c(0.05, 0.95)),3) # NB: result significant for FIGURE 2.
### R-squared and lambda medians - landbird-only dataset.
round(median(clutch_land_estimates[,6]),3) # lambda
round(median(clutch_land_estimates[,7]),3) # R-squared
### R-squared and lambda quantiles - landbird-only dataset.
round(quantile(clutch_land_estimates[,6], c(0.05, 0.95)),3) # lambda
round(quantile(clutch_land_estimates[,7], c(0.05, 0.95)),3) # R-squared


##### RESULTS -- FIGURE 1. ====
### FIGURE 1A - raincloud plots of geometric clutch size between continental and island-breeding endemics.
(Fig1A <- ggplot(data = clutch_land, aes(x = is, y = clu_geom)) + geom_boxplot(width = 0.25) +
  geom_half_point(side = 'l', range_scale = 0.2, aes(col = is)) + ylab('Clutch size\n') +
  stat_slab(width = 0.3, justification = -0.5, fill = 'black', adjust = 0.5) + ylim(0,8) + theme_minimal() +
  scale_colour_manual(values = c("#E69F00", "#0072B2")) +
  scale_x_discrete(labels = c('N' = 'Continental', 'Y' = 'Breeding island endemic')) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(family = 'serif', size = 16, color = 'black'),
        axis.title.y = element_text(family = 'serif', size = 20, color = 'black'),
        axis.text.y = element_text(family = 'serif', size = 16, color = 'black'),
        legend.position = 'none') + geom_blank())

### Prepare datasets needed to plot distributions of i) beta-coefficients, and ii) p-values across 100 models. 
clutch_land_estimates$predictor <- factor(clutch_land_estimates$predictor, levels=c("mass_geom", "modeprecocial", "migratoryresident", 
                                                                          "isY:latitude", "latitude", "isY"))
## FIGURE 1B - distributions of beta-coefficients.
(Fig1B <- ggplot(data = clutch_land_estimates %>% filter(predictor != '(Intercept)'), # removing intercept for clarity of scale.
       aes(y = predictor, x = Estimate)) + 
  stat_pointinterval() + 
  geom_vline(xintercept = 0, col = 'red') + theme_bw() + xlab('Beta coefficients') +
  scale_y_discrete(labels = c('isY' = 'Breeding island endemic', 'latitude' = 'Latitude', 'isY:latitude' = 'Latitude/island interaction',
                              'migratoryresident' = 'Non-migratory', 'modeprecocial' = 'Precociality', 'mass_geom' = 'Body mass')) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_text(family = 'serif', size = 16, color = 'black'),
        axis.title.x = element_text(family = 'serif', size = 20, color = 'black'),
        axis.text.x = element_text(family = 'serif', size = 16, color = 'black'),
        plot.margin = unit(c(1,0,0,0), 'cm')) + 
    geom_blank())

## FIGURE 1C - distributions of p-values.
(Fig1C <- ggplot(data = clutch_land_estimates %>% filter(predictor != '(Intercept)'), # removing intercept for clarity of scale.
                aes(y = predictor, x = p.value)) + 
  stat_pointinterval() + 
  geom_vline(xintercept = 0.05, col = 'red') + theme_bw() + xlab('P-value') +
  scale_y_discrete(labels = c('isY' = 'Breeding island endemic', 'latitude' = 'Latitude', 'isY:latitude' = 'Latitude/island interaction',
                              'migratoryresident' = 'Non-migratory', 'modeprecocial' = 'Precociality', 'mass_geom' = 'Body mass')) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_text(family = 'serif', size = 16, color = 'black'),
        axis.title.x = element_text(family = 'serif', size = 20, color = 'black'),
        axis.text.x = element_text(family = 'serif', size = 16, color = 'black'),
        plot.margin = unit(c(1,0,0,0), 'cm')) + 
    geom_blank())

Fig1_title <- ggdraw() + 
  draw_label("Landbird-only dataset (n = 4,288 species)", fontfamily = 'serif', size = 16,
  x = 0, hjust = 0)
Fig1 <- plot_grid(Fig1A, Fig1B, Fig1C, ncol = 1, labels = c('A)','B)','C)'), label_fontfamily = 'serif', label_size = 20)
Fig1 <- plot_grid(Fig1_title, Fig1, ncol = 1, rel_heights = c(0.05, 1))
ggsave(Fig1, filename = 'Figure1.png', device = 'png', width = 7, height = 9, units = 'in', dpi = 600)
ggsave(Fig1, filename = 'Figure1.svg', device = 'svg', width = 7, height = 9, units = 'in', dpi = 600)


##### RESULTS -- FIGURE 2. ====
# Use an example model for parameters of the trend lines in the figure:
m <- phylolm(clu_geom ~ is*latitude + migratory + mode + mass_geom, data = clutch_land, phy = trees[[1]], model = 'lambda')

(Fig2 <- ggplot(data = clutch_land %>% drop_na(), aes(y = clu_geom, x = latitude, shape = is)) + 
  geom_point(alpha = 0.7, size = 2.5, aes(col = is)) + 
  geom_smooth(method = 'lm', mapping = aes(y = predict(m, clutch_land), linetype = is), color = 'black', show.legend = F) + 
  ylim(0,8) + xlab('Latitude') +
  scale_colour_manual(name = '', values = c("#E69F00", "#0072B2"), labels = c('Continental', 'Breeding island endemic')) + 
  scale_shape_manual(name = '', values = c(19, 17), labels = c('Continental', 'Breeding island endemic'))+
  ylab('Clutch size') +
  scale_x_continuous(breaks = c(0,30,60,90), limits = c(0,90), 
                     labels = c('0' = '0°', '30' = '30°', '60' = '60°', '90' = '90°')) + 
  theme_bw() +
  theme(axis.title.y = element_text(family = 'serif', size = 20, color = 'black'),
        axis.text.y = element_text(family = 'serif', size = 16, color = 'black'),
        axis.title.x = element_text(family = 'serif', size = 20, color = 'black'),
        axis.text.x = element_text(family = 'serif', size = 16, color = 'black'),
        legend.text = element_text(family = 'serif', size = 16, color = 'black'),
        legend.position = 'top'))

Fig2_title <- ggdraw() + 
  draw_label("Landbird-only dataset (n = 4,288 species)", fontfamily = 'serif', size = 16,
             x = 0, hjust = 0)
Fig2 <- plot_grid(Fig2_title, Fig2, ncol = 1, rel_heights = c(0.05, 1))
ggsave(Fig2, filename = 'Figure2.png', device = 'png', width = 7, units = 'in', dpi = 600)
ggsave(Fig2, filename = 'Figure2.svg', device = 'svg', width = 7, units = 'in', dpi = 600)


##### MODELLING PREPARATION - AREA RELATIONSHIPS AMONG ISLAND ENDEMICS ====
## Investigate clutch size in relation to breeding range area among island endemics and prepare the dataset for Figure 3.
clutch_island <- clutch_land %>% filter(is == 'Y') # island endemics-only dataset.

## Investigate which evolutionary model is the best approximation of clutch size evolution in island-only dataset.
## Check VIF using basic lm() in R.
vif(lm(clu_geom ~ area + latitude + mode + migratory + mass_geom, data = clutch_island), type = 'predictor')

## As a reminder: area is set to log10 of area in km^2 in section MODELLING PREPARATION.
evo_fit <- function(x) {model.sel(A <- phylolm(formula = x, data = clutch_island, as.phylo(trees[[1]]), model = 'BM'),
                                  B <- phylolm(formula = x, data = clutch_island, as.phylo(trees[[1]]), model = 'OUrandomRoot'),
                                  C <- phylolm(formula = x, data = clutch_island, as.phylo(trees[[1]]), model = 'OUfixedRoot'),
                                  D <- phylolm(formula = x, data = clutch_island, as.phylo(trees[[1]]), model = 'lambda'),
                                  E <- phylolm(formula = x, data = clutch_island, as.phylo(trees[[1]]), model = 'kappa'),
                                  G <- phylolm(formula = x, data = clutch_island, as.phylo(trees[[1]]), model = 'delta'))
} # function that test 6 different evolutionary models (see METHODS)
evo_fit(clu_geom ~ area + latitude + mode + migratory + mass_geom) # for clumin lambda is the best.

## Check the health of the model and the linear model assumptions; NB: much smaller dataset so it runs faster.
cd <- comparative.data(clutch_island %>% rownames_to_column(var = 'treelabels'), phy = trees[[1]], names.col = treelabels) ## created phylogeny-ordered data.frame
nlmePagelclutch_is <- gls(clu_geom ~ area + latitude + mode + migratory + mass_geom, ## run gls with Pagel; NB: not so lengthy here, a few secs.
                          correlation = corPagel(1, phy = cd$phy),
                          data = cd$data, method = "ML")
plot(fitted(nlmePagelclutch_is), resid(nlmePagelclutch_is, type = 'normalized'))
abline(h = 0, col = "darkorange", lwd = 2)
qqnorm(resid(nlmePagelclutch, type = 'normalized'))
qqline(resid(nlmePagelclutch, type = 'normalized'))

##### MODELLING - AREA RELATIONSHIPS AMONG ISLAND ENDEMICS ====
## Fit the linear models using the same trees as before.
clutch_island_estimates <- data.frame()
for (i in 1:length(trees)) {
  m <- phylolm(clu_geom ~ area + latitude + mode + migratory + mass_geom, data = clutch_island, phy = trees[[i]], model = 'lambda')
  d <- as.data.frame(summary(m)$coefficients) %>% rownames_to_column(var = 'predictor') %>%
    mutate(lambda = m$optpar) %>% # extract Pagel's lambda for each model
    mutate(r2 = m$adj.r.squared) # extract R^squared for each model
  clutch_island_estimates <- rbind(clutch_island_estimates, d)
}
write_csv(clutch_island_estimates, 'clutch_island_estimates.csv')

##### RESULTS - AREA RELATIONSHIPS AMONG ISLAND ENDEMICS ====
## Median of beta coefficients - island landbird-only dataset.
round(median(clutch_island_estimates[clutch_island_estimates$predictor == 'area', 2]),3)
round(median(clutch_island_estimates[clutch_island_estimates$predictor == 'latitude', 2]),3) # NB: result significant for FIGURE 2.
round(median(clutch_island_estimates[clutch_island_estimates$predictor == 'modeprecocial', 2]),3)
round(median(clutch_island_estimates[clutch_island_estimates$predictor == 'migratoryresident', 2]),3)
round(median(clutch_island_estimates[clutch_island_estimates$predictor == 'mass_geom', 2]), 3)
### Quantiles of beta coefficients - island landbird-only dataset.
round(quantile(clutch_island_estimates[clutch_island_estimates$predictor == 'area', 2], c(0.05, 0.95)),3)
round(quantile(clutch_island_estimates[clutch_island_estimates$predictor == 'latitude', 2], c(0.05, 0.95)),3) # NB: result significant for FIGURE 2.
round(quantile(clutch_island_estimates[clutch_island_estimates$predictor == 'modeprecocial', 2], c(0.05, 0.95)),3)
round(quantile(clutch_island_estimates[clutch_island_estimates$predictor == 'migratoryresident', 2], c(0.05, 0.95)),3)
round(quantile(clutch_island_estimates[clutch_island_estimates$predictor == 'mass_geom', 2], c(0.05, 0.95)),3)
### Medians of standard errors - island landbird-only dataset
round(median(clutch_island_estimates[clutch_island_estimates$predictor == 'area', 3]),3)
round(median(clutch_island_estimates[clutch_island_estimates$predictor == 'latitude', 3]),3) # NB: result significant for FIGURE 2.
round(median(clutch_island_estimates[clutch_island_estimates$predictor == 'modeprecocial', 3]),3)
round(median(clutch_island_estimates[clutch_island_estimates$predictor == 'migratoryresident', 3]),3)
round(median(clutch_island_estimates[clutch_island_estimates$predictor == 'mass_geom', 3]), 3)
### Quantiles of standard errors - island landbird-only dataset.
round(quantile(clutch_island_estimates[clutch_island_estimates$predictor == 'area', 3], probs = c(0.05, 0.95)),3)
round(quantile(clutch_island_estimates[clutch_island_estimates$predictor == 'latitude', 3], c(0.05, 0.95)),3) # NB: result significant for FIGURE 2.
round(quantile(clutch_island_estimates[clutch_island_estimates$predictor == 'modeprecocial', 3], c(0.05, 0.95)),3)
round(quantile(clutch_island_estimates[clutch_island_estimates$predictor == 'migratoryresident', 3], c(0.05, 0.95)),3)
round(quantile(clutch_island_estimates[clutch_island_estimates$predictor == 'mass_geom', 3], c(0.05, 0.95)),3)
### Medians of p-values - island landbird-only dataset.
round(median(clutch_island_estimates[clutch_island_estimates$predictor == 'area', 5]),3)
round(median(clutch_island_estimates[clutch_island_estimates$predictor == 'latitude', 5]),3) # NB: result significant for FIGURE 2.
round(median(clutch_island_estimates[clutch_island_estimates$predictor == 'modeprecocial', 5]),3)
round(median(clutch_island_estimates[clutch_island_estimates$predictor == 'migratoryresident', 5]),3)
round(median(clutch_island_estimates[clutch_island_estimates$predictor == 'mass_geom', 5]), 3)
### Quantiles of p-values - island landbird-only dataset.
round(quantile(clutch_island_estimates[clutch_island_estimates$predictor == 'area', 5], probs = c(0.05, 0.95)),3)
round(quantile(clutch_island_estimates[clutch_island_estimates$predictor == 'latitude', 5], c(0.05, 0.95)),3) # NB: result significant for FIGURE 2.
round(quantile(clutch_island_estimates[clutch_island_estimates$predictor == 'modeprecocial', 5], c(0.05, 0.95)),3)
round(quantile(clutch_island_estimates[clutch_island_estimates$predictor == 'migratoryresident', 5], c(0.05, 0.95)),3)
round(quantile(clutch_island_estimates[clutch_island_estimates$predictor == 'mass_geom', 5], c(0.05, 0.95)),3)
### R-squared and lambda medians - island landbird-only dataset.
round(median(clutch_island_estimates[,6]),3) # lambda
round(median(clutch_island_estimates[,7]),3) # R-squared
### R-squared and lambda quantiles - island landbird-only dataset.
round(quantile(clutch_island_estimates[,6], c(0.05, 0.95)),3) # lambda
round(quantile(clutch_island_estimates[,7], c(0.05, 0.95)),3) # R-squared

##### RESULTS - FIGURE 3 ====
# Use an example model for parameters of the trend lines in the figure:
m <- phylolm(clu_geom ~ area + latitude + mode + migratory + mass_geom, data = clutch_island, phy = trees[[1]], model = 'lambda')

(Fig3 <- ggplot(data = clutch_island, aes(y = clu_geom, x = area)) + xlim(c(0,14)) + ylim(c(1,8)) +
  geom_jitter(alpha = 0.3, size = 2.5, col = "#0072B2") + xlab(expression(Log[10]~of~breeding~area~(m^{"2"}))) + ylab('Clutch size') +
  geom_smooth(method = 'lm', mapping = aes(y = predict(m, clutch_island)), col = 'black') + theme_bw() +
  theme(axis.title.x = element_text(family = 'serif', size = 20, color = 'black'),
        axis.text.x = element_text(family = 'serif', size = 16, color = 'black'),
        axis.title.y = element_text(family = 'serif', size = 20, color = 'black'),
        axis.text.y = element_text(family = 'serif', size = 16, color = 'black'),
        legend.position = 'none') + geom_blank())

Fig3_title <- ggdraw() + 
  draw_label("Breeding island endemic (landbirds) dataset (n = 470 species)", fontfamily = 'serif', size = 16,
             x = 0, hjust = 0)
Fig3 <- plot_grid(Fig3_title, Fig3, ncol = 1, rel_heights = c(0.05, 1))
ggsave(Fig3, filename = 'Figure3.png', device = 'png', width = 7, units = 'in', dpi = 600)
ggsave(Fig3, filename = 'Figure3.svg', device = 'svg', width = 7, units = 'in', dpi = 600)


##### DATASET PREPARATION - SEABIRDS AND THE ISLAND SYNDROME ====
## Investigate clutch size in relation to being a seabird and prepare the dataset for Figure 4.
clutch_sea <- clutch_sea %>% mutate(
  area = log10(area), # log area 
  latitude = abs(latitude)) %>% # force latitude to absolute
  column_to_rownames(var = 'treelabels')
clutch_land_and_sea <- clutch_land_and_sea %>% mutate( 
  area = log10(area), # log area 
  latitude = abs(latitude)) %>% # force latitude to absolute
  column_to_rownames(var = 'treelabels')

## Check VIF using basic lm() in R.
vif(lm(clu_geom ~ is*latitude + migratory + mode + mass_geom, data = clutch_sea), type = 'predictor')

## Investigate which evolutionary model is the best approximation of clutch size evolution in seabird-only dataset.
evo_fit <- function(x) {model.sel(A <- phylolm(formula = x, data = clutch_sea, as.phylo(trees[[1]]), model = 'BM'),
                                  B <- phylolm(formula = x, data = clutch_sea, as.phylo(trees[[1]]), model = 'OUrandomRoot'),
                                  C <- phylolm(formula = x, data = clutch_sea, as.phylo(trees[[1]]), model = 'OUfixedRoot'),
                                  D <- phylolm(formula = x, data = clutch_sea, as.phylo(trees[[1]]), model = 'lambda'),
                                  E <- phylolm(formula = x, data = clutch_sea, as.phylo(trees[[1]]), model = 'kappa'),
                                  G <- phylolm(formula = x, data = clutch_sea, as.phylo(trees[[1]]), model = 'delta'))
} # function that test 6 different evolutionary models (see METHODS)
evo_fit(clu_geom ~ is*latitude + migratory + mode + mass_geom) # lambda remains the best model.

## Check the health of the model and the linear model assumptions; NB: much smaller dataset so it runs faster.
cd <- comparative.data(clutch_sea %>% rownames_to_column(var = 'treelabels'), phy = trees[[1]], names.col = treelabels) ## created phylogeny-ordered data.frame
nlmePagelclutch_sea <- gls(clu_geom ~ is*latitude + migratory + mode + mass_geom, 
                       correlation = corPagel(1, phy = cd$phy),
                       data = cd$data, method = "ML")
plot(fitted(nlmePagelclutch_sea), resid(nlmePagelclutch_sea, type = 'normalized'))  # plots normalised residuals against fitted.
abline(h = 0, col = "darkorange", lwd = 2) # 0 line to compare.
qqnorm(resid(nlmePagelclutch, type = 'normalized')) # quantiles of residuals.
qqline(resid(nlmePagelclutch, type = 'normalized')) # expected trend line of residuals.

##### MODELLING - ISLAND SYNDROME IN SEABIRDS ====
clutch_sea_estimates <- data.frame()
for (i in 1:length(trees)) {
  m <- phylolm(clu_geom ~ is*latitude + mode + migratory + mass_geom, data = clutch_sea, phy = trees[[i]], model = 'lambda')
  d <- as.data.frame(summary(m)$coefficients) %>% rownames_to_column(var = 'predictor') %>%
    mutate(lambda = m$optpar) %>% # extract Pagel's lambda for each model
    mutate(r2 = m$adj.r.squared) # extract R^squared for each model
  clutch_sea_estimates <- rbind(clutch_sea_estimates, d)
}
write_csv(clutch_sea_estimates, 'clutch_sea_estimates.csv')

##### RESULTS - ISLAND SYNDROME IN SEABIRDS ====
## Median of beta coefficients - seabird-only dataset.
round(median(clutch_sea_estimates[clutch_sea_estimates$predictor == 'isY', 2]),3)
round(median(clutch_sea_estimates[clutch_sea_estimates$predictor == 'latitude', 2]),3) # NB: result significant for FIGURE 2.
round(median(clutch_sea_estimates[clutch_sea_estimates$predictor == 'modeprecocial', 2]),3)
round(median(clutch_sea_estimates[clutch_sea_estimates$predictor == 'migratoryresident', 2]),3)
round(median(clutch_sea_estimates[clutch_sea_estimates$predictor == 'mass_geom', 2]), 3)
round(median(clutch_sea_estimates[clutch_sea_estimates$predictor == 'seabirdseabird', 2]), 3) # NB: result significant for FIGURE 2.
### Quantiles of beta coefficients - seabird-only dataset.
round(quantile(clutch_sea_estimates[clutch_sea_estimates$predictor == 'isY', 2], c(0.05, 0.95)),3)
round(quantile(clutch_sea_estimates[clutch_sea_estimates$predictor == 'latitude', 2], c(0.05, 0.95)),3) # NB: result significant for FIGURE 2.
round(quantile(clutch_sea_estimates[clutch_sea_estimates$predictor == 'modeprecocial', 2], c(0.05, 0.95)),3)
round(quantile(clutch_sea_estimates[clutch_sea_estimates$predictor == 'migratoryresident', 2], c(0.05, 0.95)),3)
round(quantile(clutch_sea_estimates[clutch_sea_estimates$predictor == 'mass_geom', 2], c(0.05, 0.95)),3)
round(quantile(clutch_sea_estimates[clutch_sea_estimates$predictor == 'isY:latitude', 2], c(0.05, 0.95)) ,3)# NB: result significant for FIGURE 2.
### Medians of standard errors - seabird-only dataset
round(median(clutch_sea_estimates[clutch_sea_estimates$predictor == 'isY', 3]),3)
round(median(clutch_sea_estimates[clutch_sea_estimates$predictor == 'latitude', 3]),3) # NB: result significant for FIGURE 2.
round(median(clutch_sea_estimates[clutch_sea_estimates$predictor == 'modeprecocial', 3]),3)
round(median(clutch_sea_estimates[clutch_sea_estimates$predictor == 'migratoryresident', 3]),3)
round(median(clutch_sea_estimates[clutch_sea_estimates$predictor == 'mass_geom', 3]), 3)
round(median(clutch_sea_estimates[clutch_sea_estimates$predictor == 'isY:latitude', 3]), 3) # NB: result significant for FIGURE 2.
### Quantiles of standard errors - seabird-only dataset.
round(quantile(clutch_sea_estimates[clutch_sea_estimates$predictor == 'isY', 3], probs = c(0.05, 0.95)),3)
round(quantile(clutch_sea_estimates[clutch_sea_estimates$predictor == 'latitude', 3], c(0.05, 0.95)),3) # NB: result significant for FIGURE 2.
round(quantile(clutch_sea_estimates[clutch_sea_estimates$predictor == 'modeprecocial', 3], c(0.05, 0.95)),3)
round(quantile(clutch_sea_estimates[clutch_sea_estimates$predictor == 'migratoryresident', 3], c(0.05, 0.95)),3)
round(quantile(clutch_sea_estimates[clutch_sea_estimates$predictor == 'mass_geom', 3], c(0.05, 0.95)),3)
round(quantile(clutch_sea_estimates[clutch_sea_estimates$predictor == 'isY:latitude', 3], c(0.05, 0.95)),3) # NB: result significant for FIGURE 2.
### Medians of p-values - seabird-only dataset.
round(median(clutch_sea_estimates[clutch_sea_estimates$predictor == 'isY', 5]),3)
round(median(clutch_sea_estimates[clutch_sea_estimates$predictor == 'latitude', 5]),3) # NB: result significant for FIGURE 2.
round(median(clutch_sea_estimates[clutch_sea_estimates$predictor == 'modeprecocial', 5]),3)
round(median(clutch_sea_estimates[clutch_sea_estimates$predictor == 'migratoryresident', 5]),3)
round(median(clutch_sea_estimates[clutch_sea_estimates$predictor == 'mass_geom', 5]), 3)
round(median(clutch_sea_estimates[clutch_sea_estimates$predictor == 'isY:latitude', 5]), 3) # NB: result significant for FIGURE 2.
### Quantiles of p-values - seabird-only dataset.
round(quantile(clutch_sea_estimates[clutch_sea_estimates$predictor == 'isY', 5], probs = c(0.05, 0.95)),3)
round(quantile(clutch_sea_estimates[clutch_sea_estimates$predictor == 'latitude', 5], c(0.05, 0.95)),3) # NB: result significant for FIGURE 2.
round(quantile(clutch_sea_estimates[clutch_sea_estimates$predictor == 'modeprecocial', 5], c(0.05, 0.95)),3)
round(quantile(clutch_sea_estimates[clutch_sea_estimates$predictor == 'migratoryresident', 5], c(0.05, 0.95)),3)
round(quantile(clutch_sea_estimates[clutch_sea_estimates$predictor == 'mass_geom', 5], c(0.05, 0.95)),3)
round(quantile(clutch_sea_estimates[clutch_sea_estimates$predictor == 'isY:latitude', 5], c(0.05, 0.95)),3) # NB: result significant for FIGURE 2.
### R-squared and lambda medians - seabird-only dataset.
round(median(clutch_sea_estimates[,6]),3) # lambda
round(median(clutch_sea_estimates[,7]),3) # R-squared
### R-squared and lambda quantiles - seabird-only dataset.
round(quantile(clutch_sea_estimates[,6], c(0.05, 0.95)),3) # lambda
round(quantile(clutch_sea_estimates[,7], c(0.05, 0.95)),3) # R-squared

#### MODELLING PREPARATION - INTERACTION BETWEEN SEABIRD AND BREEDING ISLAND ENDEMICITY ====
## Check VIF using basic lm() in R - new variable added so best to check.
vif(lm(clu_geom ~ is*seabird + latitude + migratory + mode + mass_geom, data = clutch_land_and_sea), type = 'predictor')

## Investigate which evolutionary model is the best approximation of clutch size evolution in the complete dataset.
evo_fit <- function(x) {model.sel(A <- phylolm(formula = x, data = clutch_land_and_sea, as.phylo(trees[[1]]), model = 'BM'),
                                  B <- phylolm(formula = x, data = clutch_land_and_sea, as.phylo(trees[[1]]), model = 'OUrandomRoot'),
                                  C <- phylolm(formula = x, data = clutch_land_and_sea, as.phylo(trees[[1]]), model = 'OUfixedRoot'),
                                  D <- phylolm(formula = x, data = clutch_land_and_sea, as.phylo(trees[[1]]), model = 'lambda'),
                                  E <- phylolm(formula = x, data = clutch_land_and_sea, as.phylo(trees[[1]]), model = 'kappa'),
                                  G <- phylolm(formula = x, data = clutch_land_and_sea, as.phylo(trees[[1]]), model = 'delta'))
} # function that test 6 different evolutionary models (see METHODS)
evo_fit(clu_geom ~ is*latitude + migratory + mode + mass_geom) # lambda remains the best model.

## Check the health of the model and the linear model assumptions; NB: much smaller dataset so it runs faster.
cd <- comparative.data(clutch_land_and_sea %>% rownames_to_column(var = 'treelabels'), phy = trees[[1]], names.col = treelabels) ## created phylogeny-ordered data.frame
nlmePagelclutch_land_and_sea <- gls(clu_geom ~ is*seabird + latitude + migratory + mode + mass_geom, 
                           correlation = corPagel(1, phy = cd$phy),
                           data = cd$data, method = "ML")
plot(fitted(nlmePagelclutch_sea), resid(nlmePagelclutch_sea, type = 'normalized'))  # plots normalised residuals against fitted.
abline(h = 0, col = "darkorange", lwd = 2) # 0 line to compare.
qqnorm(resid(nlmePagelclutch, type = 'normalized')) # quantiles of residuals.
qqline(resid(nlmePagelclutch, type = 'normalized')) # expected trend line of residuals.

##### MODELLING - ISLAND SYNDROME IN SEABIRDS ====
clutch_landsea_estimates <- data.frame()
for (i in 1:length(trees)) {
  m <- phylolm(clu_geom ~ is*seabird + latitude + mode + migratory + mass_geom, data = clutch_land_and_sea, phy = trees[[i]], model = 'lambda')
  d <- as.data.frame(summary(m)$coefficients) %>% rownames_to_column(var = 'predictor') %>%
    mutate(lambda = m$optpar) %>% # extract Pagel's lambda for each model
    mutate(r2 = m$adj.r.squared) # extract R^squared for each model
  clutch_landsea_estimates <- rbind(clutch_landsea_estimates, d)
}
write_csv(clutch_landsea_estimates, 'clutch_land_and_sea_estimates.csv')

##### RESULTS - INTERACTION OF THE ISLAND SYNDROME AND BEING A SEABIRD ====
## Median of beta coefficients - complete dataset.
round(median(clutch_landsea_estimates[clutch_landsea_estimates$predictor == 'isY', 2]),3)
round(median(clutch_landsea_estimates[clutch_landsea_estimates$predictor == 'latitude', 2]),3) # NB: result significant for FIGURE 2.
round(median(clutch_landsea_estimates[clutch_landsea_estimates$predictor == 'modeprecocial', 2]),3)
round(median(clutch_landsea_estimates[clutch_landsea_estimates$predictor == 'migratoryresident', 2]),3)
round(median(clutch_landsea_estimates[clutch_landsea_estimates$predictor == 'mass_geom', 2]), 3)
round(median(clutch_landsea_estimates[clutch_landsea_estimates$predictor == 'seabirdseabird', 2]), 3) # NB: result significant for FIGURE 2.
round(median(clutch_landsea_estimates[clutch_landsea_estimates$predictor == 'isY:seabirdseabird', 2]), 3) # NB: result significant for FIGURE 2.
### Quantiles of beta coefficients - complete dataset.
round(quantile(clutch_landsea_estimates[clutch_landsea_estimates$predictor == 'isY', 2], c(0.05, 0.95)),3)
round(quantile(clutch_landsea_estimates[clutch_landsea_estimates$predictor == 'latitude', 2], c(0.05, 0.95)),3) # NB: result significant for FIGURE 2.
round(quantile(clutch_landsea_estimates[clutch_landsea_estimates$predictor == 'modeprecocial', 2], c(0.05, 0.95)),3)
round(quantile(clutch_landsea_estimates[clutch_landsea_estimates$predictor == 'migratoryresident', 2], c(0.05, 0.95)),3)
round(quantile(clutch_landsea_estimates[clutch_landsea_estimates$predictor == 'mass_geom', 2], c(0.05, 0.95)),3)
round(quantile(clutch_landsea_estimates[clutch_landsea_estimates$predictor == 'seabirdseabird', 2], c(0.05, 0.95)) ,3)# NB: result significant for FIGURE 2.
round(quantile(clutch_landsea_estimates[clutch_landsea_estimates$predictor == 'isY:seabirdseabird', 2], c(0.05, 0.95)) ,3)# NB: result significant for FIGURE 2.
### Medians of standard errors - complete dataset
round(median(clutch_landsea_estimates[clutch_landsea_estimates$predictor == 'isY', 3]),3)
round(median(clutch_landsea_estimates[clutch_landsea_estimates$predictor == 'latitude', 3]),3) # NB: result significant for FIGURE 2.
round(median(clutch_landsea_estimates[clutch_landsea_estimates$predictor == 'modeprecocial', 3]),3)
round(median(clutch_landsea_estimates[clutch_landsea_estimates$predictor == 'migratoryresident', 3]),3)
round(median(clutch_landsea_estimates[clutch_landsea_estimates$predictor == 'mass_geom', 3]), 3)
round(median(clutch_landsea_estimates[clutch_landsea_estimates$predictor == 'seabirdseabird', 3]), 3) # NB: result significant for FIGURE 2.
round(median(clutch_landsea_estimates[clutch_landsea_estimates$predictor == 'isY:seabirdseabird', 3]), 3) # NB: result significant for FIGURE 2.
### Quantiles of standard errors - complete dataset.
round(quantile(clutch_landsea_estimates[clutch_landsea_estimates$predictor == 'isY', 3], probs = c(0.05, 0.95)),3)
round(quantile(clutch_landsea_estimates[clutch_landsea_estimates$predictor == 'latitude', 3], c(0.05, 0.95)),3) # NB: result significant for FIGURE 2.
round(quantile(clutch_landsea_estimates[clutch_landsea_estimates$predictor == 'modeprecocial', 3], c(0.05, 0.95)),3)
round(quantile(clutch_landsea_estimates[clutch_landsea_estimates$predictor == 'migratoryresident', 3], c(0.05, 0.95)),3)
round(quantile(clutch_landsea_estimates[clutch_landsea_estimates$predictor == 'mass_geom', 3], c(0.05, 0.95)),3)
round(quantile(clutch_landsea_estimates[clutch_landsea_estimates$predictor == 'seabirdseabird', 3], c(0.05, 0.95)),3) # NB: result significant for FIGURE 2.
round(quantile(clutch_landsea_estimates[clutch_landsea_estimates$predictor == 'isY:seabirdseabird', 3], c(0.05, 0.95)),3) # NB: result significant for FIGURE 2.
### Medians of p-values - complete dataset.
round(median(clutch_landsea_estimates[clutch_landsea_estimates$predictor == 'isY', 5]),3)
round(median(clutch_landsea_estimates[clutch_landsea_estimates$predictor == 'latitude', 5]),3) # NB: result significant for FIGURE 2.
round(median(clutch_landsea_estimates[clutch_landsea_estimates$predictor == 'modeprecocial', 5]),3)
round(median(clutch_landsea_estimates[clutch_landsea_estimates$predictor == 'migratoryresident', 5]),3)
round(median(clutch_landsea_estimates[clutch_landsea_estimates$predictor == 'mass_geom', 5]), 3)
round(median(clutch_landsea_estimates[clutch_landsea_estimates$predictor == 'seabirdseabird', 5]), 3) # NB: result significant for FIGURE 2.
round(median(clutch_landsea_estimates[clutch_landsea_estimates$predictor == 'isY:seabirdseabird', 5]), 3) # NB: result significant for FIGURE 2.
### Quantiles of p-values - complete dataset.
round(quantile(clutch_landsea_estimates[clutch_landsea_estimates$predictor == 'isY', 5], probs = c(0.05, 0.95)),3)
round(quantile(clutch_landsea_estimates[clutch_landsea_estimates$predictor == 'latitude', 5], c(0.05, 0.95)),3) # NB: result significant for FIGURE 2.
round(quantile(clutch_landsea_estimates[clutch_landsea_estimates$predictor == 'modeprecocial', 5], c(0.05, 0.95)),3)
round(quantile(clutch_landsea_estimates[clutch_landsea_estimates$predictor == 'migratoryresident', 5], c(0.05, 0.95)),3)
round(quantile(clutch_landsea_estimates[clutch_landsea_estimates$predictor == 'mass_geom', 5], c(0.05, 0.95)),3)
round(quantile(clutch_landsea_estimates[clutch_landsea_estimates$predictor == 'seabirdseabird', 5], c(0.05, 0.95)),3) # NB: result significant for FIGURE 2.
round(quantile(clutch_landsea_estimates[clutch_landsea_estimates$predictor == 'isY:seabirdseabird', 5], c(0.05, 0.95)),3) # NB: result significant for FIGURE 2.
### R-squared and lambda medians - complete dataset.
round(median(clutch_landsea_estimates[,6]),3) # lambda
round(median(clutch_landsea_estimates[,7]),3) # R-squared
### R-squared and lambda quantiles - complete dataset.
round(median(clutch_landsea_estimates[,6], c(0.05, 0.95)),3) # lambda
round(median(clutch_landsea_estimates[,7], c(0.05, 0.95)),3) # R-squared
#### RESULTS - FIGURE 4 =====

#Modify the dataset for plotting purposes - create a variable to describe mix of continental + islands
clutch_land_and_sea <- clutch_land_and_sea %>% mutate(plot_var = ifelse(seabird == 'seabird' & is == 'Y', 'IslSea',
                                                                        ifelse(seabird == 'seabird' & is == 'N', 'ConSea',
                                                                               ifelse(seabird == 'landbird' & is == 'Y', 'IslLan', 'ConLan'))))

(Fig4 <- ggplot(data = clutch_land_and_sea, aes(x = plot_var, y = clu_geom, col = plot_var)) + geom_boxplot(width = 0.25) + 
    ylab('Clutch size\n') +
    stat_slab(width = 0.3, justification = -0.5, aes(fill = is), alpha = 0.3) + ylim(0,8) +
    scale_colour_manual(values = c("#E69F00", "#E69F00", "#0072B2", "#0072B2")) +
    scale_x_discrete(labels = c('ConLan' = 'Continental landbird', 'ConSea' = 'Continental seabird', 
                                'IslLan' = 'Island landbird', 'IslSea' = 'Island seabird')) +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(family = 'serif', size = 12, color = 'black'),
          axis.title.y = element_text(family = 'serif', size = 20, color = 'black'),
          axis.text.y = element_text(family = 'serif', size = 16, color = 'black'),
          legend.position = 'none') + geom_blank())

Fig4_title <- ggdraw() + 
  draw_label("Complete dataset (n = 4,530 species)", fontfamily = 'serif', size = 16,
             x = 0, hjust = 0)
Fig4 <- plot_grid(Fig4_title, Fig4, ncol = 1, rel_heights = c(0.05, 1))
ggsave(Fig4, filename = 'Figure4.png', device = 'png', width = 7, units = 'in', dpi = 600)
ggsave(Fig4, filename = 'Figure4.svg', device = 'svg', width = 7, units = 'in', dpi = 600)


