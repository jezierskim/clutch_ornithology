#### CLUTCH PAPER - ANALYSIS CODE
library(tidyverse)
library(ape)
library(phylolm)
library(nlme)
library(car)
library(MuMIn)
library(caper)

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


##### RESULTS - GENERAL STATISTICS ====
clutch %>% group_by(is) %>% summarise(min = mean(clumin), max = mean(clumax), mid = mean(clumid))
### calculate means of each of clutch size metrics across the two states of island endemicity.

##### MODELLING PREPARATION ====
### Firstly, prepare variables: calculate natural logs of maximum body mass (bodmax) and area; check VIF; run
### different evolutionary models to see which one fits the data best; run a PGLS model with the best evolutionary model
### for each metric of clutch size to examine model diagnostics.

## Transform variables appropriately.
clutch <- clutch %>% mutate(bodmax = log(bodmax), # log maximum body mass 
                            area = log(area), # log area 
                            latitude = abs(latitude)) %>% # force latitude to absolute
  column_to_rownames(var = 'jetz') %>% # move phylogeny names to rownames for downstream use with phylolm.

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





### FIGURE 1 - MODELS
clumin_estimates <- data.frame()

for (i in 1:length(trees)) {
  m <- phylolm(clumin ~ latitude*is + bodmax, phy = trees[[i]], model = 'lambda')
  
}


intercept            <-
  phylolm::summary.phylolm(mod)$coefficients[[1, 1]]
se.intercept         <-
  phylolm::summary.phylolm(mod)$coefficients[[1, 2]]
estimate                <-
  phylolm::summary.phylolm(mod)$coefficients[[2, 1]]
se.estimate             <-
  phylolm::summary.phylolm(mod)$coefficients[[2, 2]]
pval.intercept       <-
  phylolm::summary.phylolm(mod)$coefficients[[1, 4]]
pval.estimate           <-
  phylolm::summary.phylolm(mod)$coefficients[[2, 4]]
