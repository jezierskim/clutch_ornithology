#### CLUTCH PAPER - ANALYSIS CODE
library(tidyverse)
library(ape)
library(phylolm)
library(nlme)
library(car)
library(MuMIn)
library(caper)
library(ggdist)
library(gghalves)

##### LOAD IN THE DATA REQUIRED: ====
## Load clutch dataset; create clutch-midpoint variable; remove outliers for all three values.
clutch <- read_csv('~/Desktop/clutch_dataset_geo.csv') # load in the clutch dataset.
clutch <- clutch %>% mutate(clumid = (clumin + clumax)/2) # calculate midpoint
Qmin <- quantile(clutch$clumin, probs=c(.25, .75), na.rm = FALSE) # calculate 25% and 75% quantiles.
Qmax <- quantile(clutch$clumax, probs=c(.25, .75), na.rm = FALSE)
Qmid <- quantile(clutch$clumid, probs=c(.25, .75), na.rm = FALSE)
IQRmin <- IQR(clutch$clumin) # calculate inter-quantile range
IQRmax <- IQR(clutch$clumax)
IQRmid <- IQR(clutch$clumid)
clutch <- clutch %>%
  filter(clumin < (Qmin[2]+1.5*IQRmin)) %>% # remove upper outliers
  filter(clumax < (Qmax[2]+1.5*IQRmax)) %>%
  filter(clumid < (Qmid[2]+1.5*IQRmid)) 
## END RESULT: datatframe of 4,522 species to be retained in the analysis.
rm(IQRmax, IQRmid, IQRmin, Qmax, Qmid, Qmin)

## Load phylogenetic trees, prune them to match the dataset, sample 100 random trees to be used in the analysis.
trees <- read.tree('~/Desktop/randomHackSix.tre') # load in the sample of 1,000 trees based on genetic data.
trees <- lapply(trees, keep.tip, tip = clutch$jetz) # prune all the trees to include species from 'clutch'.
trees <- sample(trees, 100) # sample 100 random trees without replacement. 


##### RESULTS -- GENERAL STATISTICS ====
clutch %>% group_by(is) %>% summarise(min = mean(clumin), max = mean(clumax), mid = mean(clumid))
### calculate means of each of clutch size metrics across the two states of island endemicity.

##### MODELLING PREPARATION ====
### Firstly, prepare variables: calculate natural logs of maximum body mass (bodmax) and area; check VIF; run
### different evolutionary models to see which one fits the data best; run a PGLS model with the best evolutionary model
### for each metric of clutch size to examine model diagnostics.

## Transform variables appropriately.
clutch <- clutch %>% mutate(bodmax = log10(bodmax), # log maximum body mass 
                            area = log10(area), # log area 
                            latitude = abs(latitude)) %>% # force latitude to absolute
  column_to_rownames(var = 'jetz') # move phylogeny names to rownames for downstream use with phylolm.

## Check VIF using basic lm() in R.
vif(lm(clumin ~ latitude*is + bodmax, data = clutch), type = 'predictor')

## Examine the best evolutionary model using different phylolm() run of the model.
evo_fit <- function(x) {model.sel(A <- phylolm(formula = x, data = clutch, as.phylo(trees[[1]]), model = 'BM'),
                       B <- phylolm(formula = x, data = clutch, phy = as.phylo(trees[[1]]), model = 'OUrandomRoot'),
                       C <- phylolm(formula = x, data = clutch, as.phylo(trees[[1]]), model = 'OUfixedRoot'),
                       D <- phylolm(formula = x, data = clutch, as.phylo(trees[[1]]), model = 'lambda'),
                       E <- phylolm(formula = x, data = clutch, as.phylo(trees[[1]]), model = 'kappa'),
                       G <- phylolm(formula = x, data = clutch, as.phylo(trees[[1]]), model = 'delta'))
} # function that test 6 different evolutionary models (see METHODS)
evo_fit(clumin ~ latitude*is + bodmax) # for clumin lambda is the best.
evo_fit(clumax ~ latitude*is + bodmax) # for clumax lambda is the best.
evo_fit(clumid ~ latitude*is + bodmax) # for clumid lambda is the best.

## Examine the model diagnostics using nlme() and Pagel's correlation matrix (Pagel's lambda).
clutch1 <- clutch %>% rownames_to_column(var = 'jetz') ## need to return rownames to a column for functionality.
cd <- comparative.data(clutch1, phy = trees[[1]], names.col = jetz) ## created phylogeny-ordered data.frame
nlmePagelclumin <- gls(clumin ~ latitude*is + bodmax, ## run gls with Pagel; NB: THIS STEP TAKES A WHILE.
                       correlation = corPagel(1, phy = cd$phy),
                 data = cd$data, method = "ML")
write_rds(nlmePagelclumin, '~/Desktop/nlmePagelclumin.rds')
plot(fitted(nlmePagelclumin), resid(nlmePagelclumin, type = 'normalized'))
abline(h = 0, col = "darkorange", lwd = 2)
qqnorm(resid(nlmePagelclumin, type = 'normalized'))
qqline(resid(nlmePagelclumin, type = 'normalized'))

nlmePagelclumax <- gls(clumax ~ latitude*is + bodmax, ## run gls with Pagel; NB: THIS STEP TAKES A WHILE.
                       correlation = corPagel(1, phy = cd$phy),
                       data = cd$data, method = "ML")
write_rds(nlmePagelclumax, '~/Desktop/nlmePagelclumax.rds')
plot(fitted(nlmePagelclumax), resid(nlmePagelclumin, type = 'normalized'))
abline(h = 0, col = "darkorange", lwd = 2)
qqnorm(resid(nlmePagelclumax, type = 'normalized'))
qqline(resid(nlmePagelclumax, type = 'normalized'))

nlmePagelclumid <- gls(clumid ~ latitude*is + bodmax, ## run gls with Pagel; NB: THIS STEP TAKES A WHILE.
                       correlation = corPagel(1, phy = cd$phy),
                       data = cd$data, method = "ML")
write_rds(nlmePagelclumid, '~/Desktop/nlmePagelclumid.rds')
plot(fitted(nlmePagelclumid), resid(nlmePagelclumin, type = 'normalized'))
abline(h = 0, col = "darkorange", lwd = 2)
qqnorm(resid(nlmePagelclumid, type = 'normalized'))
qqline(resid(nlmePagelclumid, type = 'normalized'))



##### MODELLING -- RUN ACROSS 100 TREES. ======
clumin_estimates <- data.frame()
for (i in 1:length(trees)) {
  m <- phylolm(clumin ~ latitude*is + bodmax, data = clutch, phy = trees[[i]], model = 'lambda')
  d <- as.data.frame(summary(m)$coefficients) %>% rownames_to_column(var = 'predictor') %>%
    mutate(lambda = m$optpar) %>% # extract Pagel's lambda for each model
    mutate(r2 = m$adj.r.squared) # extract R^squared for each model
  clumin_estimates <- rbind(clumin_estimates, d)
}
write_csv(clumin_estimates, '~/Desktop/clumin_estimates.csv')

clumax_estimates <- data.frame()
for (i in 1:length(trees)) {
  m <- phylolm(clumax ~ latitude*is + bodmax, data = clutch, phy = trees[[i]], model = 'lambda')
  d <- as.data.frame(summary(m)$coefficients) %>% rownames_to_column(var = 'predictor') %>%
    mutate(lambda = m$optpar) %>% # extract Pagel's lambda for each model
    mutate(r2 = m$adj.r.squared) # extract R^squared for each model
  clumax_estimates <- rbind(clumax_estimates, d)
}
write_csv(clumax_estimates, '~/Desktop/clumax_estimates.csv')

clumid_estimates <- data.frame()
for (i in 1:length(trees)) {
  m <- phylolm(clumid ~ latitude*is + bodmax, data = clutch, phy = trees[[i]], model = 'lambda')
  d <- as.data.frame(summary(m)$coefficients) %>% rownames_to_column(var = 'predictor') %>% 
    mutate(lambda = m$optpar) %>% # extract Pagel's lambda for each model
    mutate(r2 = m$adj.r.squared) # extract R^squared for each model
  clumid_estimates <- rbind(clumid_estimates, d)
}
write_csv(clumid_estimates, '~/Desktop/clumid_estimates.csv')

##### RESULTS -- MODEL ESTIMATES ASSOCIATED WITH FIGURE 1. ====
### For minimum clutch size:
summary(clumin_estimates %>% filter(predictor == 'isY'))
summary(clumin_estimates %>% filter(predictor == 'bodmax'))
summary(clumin_estimates %>% filter(predictor == 'lambda'))
summary(clumin_estimates %>% filter(predictor == 'r2'))

### For maximum clutch size:
summary(clumax_estimates %>% filter(predictor == 'isY'))
summary(clumax_estimates %>% filter(predictor == 'bodmax'))
summary(clumax_estimates %>% filter(predictor == 'lambda'))
summary(clumax_estimates %>% filter(predictor == 'r2'))

### For mid-point clutch size:
summary(clumid_estimates %>% filter(predictor == 'isY'))
summary(clumid_estimates %>% filter(predictor == 'bodmax'))
summary(clumid_estimates %>% filter(predictor == 'lambda'))
summary(clumid_estimates %>% filter(predictor == 'r2'))

##### RESULTS -- FIGURE 1. ====
### Raincloud plots forming part A) of Figure 1.
ggplot(data = clutch, aes(x = is, y = clumin)) + geom_boxplot(width = 0.25) +
  geom_half_point(side = 'l', range_scale = 0.2, aes(col = is)) +
  stat_slab(width = 0.3, justification = -0.5, fill = 'black') + ylim(0,10) + theme_minimal() +
  scale_colour_manual(values = c("#E69F00", "#0072B2"))

ggplot(data = clutch, aes(x = is, y = clumax)) + geom_boxplot(width = 0.25) +
  geom_half_point(side = 'l', range_scale = 0.2, aes(col = is)) +
  stat_slab(width = 0.3, justification = -0.5, fill = 'black') + ylim(0,10) + theme_minimal() +
  scale_colour_manual(values = c("#E69F00", "#0072B2"))

ggplot(data = clutch, aes(x = is, y = clumid)) + geom_boxplot(width = 0.25) +
  geom_half_point(side = 'l', range_scale = 0.2, aes(col = is)) +
  stat_slab(width = 0.3, justification = -0.5, fill = 'black') + ylim(0,10) + theme_minimal() +
  scale_colour_manual(values = c("#E69F00", "#0072B2"))

### Distribution plots of estimates and p-values of part B) i) and B) ii) of Figure 1.
clumin_estimates$predictor <- factor(clumin_estimates$predictor, levels=c("bodmax", "latitude:isY", "latitude", "isY"))
clumax_estimates$predictor <- factor(clumax_estimates$predictor, levels=c("bodmax", "latitude:isY", "latitude", "isY"))
clumid_estimates$predictor <- factor(clumid_estimates$predictor, levels=c("bodmax", "latitude:isY", "latitude", "isY"))
## Estimate distribution - part i).
ggplot(data = clumin_estimates %>% filter(predictor != '(Intercept)'), # removing intercept for clarity of scale.
       aes(y = predictor, x = Estimate)) + stat_pointinterval() + geom_vline(xintercept = 0, col = 'red') +
  theme_minimal() + xlim(-0.4, 0.1)

ggplot(data = clumax_estimates %>% filter(predictor != '(Intercept)'), # removing intercept for clarity of scale.
       aes(y = predictor, x = Estimate)) + stat_pointinterval() + geom_vline(xintercept = 0, col = 'red') +
  theme_minimal() + xlim(-0.4, 0.1)

ggplot(data = clumid_estimates %>% filter(predictor != '(Intercept)'), # removing intercept for clarity of scale.
       aes(y = predictor, x = Estimate)) + stat_pointinterval() + geom_vline(xintercept = 0, col = 'red') +
  theme_minimal() + xlim(-0.4, 0.1)

## P-value distributions - part ii)

ggplot(data = clumin_estimates %>% filter(predictor != '(Intercept)'), # removing intercept for clarity of scale.
       aes(y = predictor, x = p.value)) + stat_pointinterval() + geom_vline(xintercept = 0.05, col = 'red') +
  theme_minimal() + xlim(0,1)

ggplot(data = clumax_estimates %>% filter(predictor != '(Intercept)'), # removing intercept for clarity of scale.
       aes(y = predictor, x = p.value)) + stat_pointinterval() + geom_vline(xintercept = 0.05, col = 'red') +
  theme_minimal() + xlim(0,1)

ggplot(data = clumid_estimates %>% filter(predictor != '(Intercept)'), # removing intercept for clarity of scale.
       aes(y = predictor, x = p.value)) + stat_pointinterval() + geom_vline(xintercept = 0.05, col = 'red') +
  theme_minimal() + xlim(0,1)




##### RESULTS -- LATITUDE RELATED RESULT =====
### For minimum clutch size:
summary(clumin_estimates %>% filter(predictor == 'latitude'))
summary(clumin_estimates %>% filter(predictor == 'latitude:isY'))

### For maximum clutch size:
summary(clumax_estimates %>% filter(predictor == 'latitude'))
summary(clumax_estimates %>% filter(predictor == 'latitude:isY'))

### For mid-point clutch size:
summary(clumid_estimates %>% filter(predictor == 'latitude'))
summary(clumid_estimates %>% filter(predictor == 'latitude:isY'))

##### RESULTS -- FIGURE 2. ====
### Plot the figures from Figure 2. in the manuscript.
### Plot clutch sizes along latitude, and add two lines which assume body mass to be constant.
## Maximum clutch size
# Use an example model for parameters of the trend lines in the figure:
m <- phylolm(clumax ~ latitude*is + bodmax, data = clutch, phy = trees[[1]], model = 'lambda')
ggplot(data = clutch %>% drop_na(), aes(y = clumax, x = latitude, col = is)) + geom_jitter(alpha = 0.3, size = 2.5) + 
  geom_smooth(method = 'lm', mapping = aes(y = predict(m, clutch))) + 
  theme_minimal() + ylim(0,10) + scale_colour_manual(values = c("#E69F00", "#0072B2")) +
  scale_x_continuous(breaks = c(0,30,60,90), limits = c(0,90)) +
  theme(legend.position = 'none')

## Midpoint clutch size
# Use an example model for parameters of the trend lines in the figure:
m <- phylolm(clumid ~ latitude*is + bodmax, data = clutch, phy = trees[[1]], model = 'lambda')
ggplot(data = clutch %>% drop_na(), aes(y = clumid, x = latitude, col = is)) + geom_jitter(alpha = 0.3, size = 2.5) + 
  geom_smooth(method = 'lm', mapping = aes(y = predict(m, clutch))) + 
  theme_minimal() + ylim(0,10) + scale_colour_manual(values = c("#E69F00", "#0072B2")) +
  scale_x_continuous(breaks = c(0,30,60,90), limits = c(0,90)) +
  theme(legend.position = 'none')

