# Multivariate generalized mixed effect model approach to estimating effect of year
# Jenilee Gobin
# Aug 8, 2022


# load libraries and preliminaries -----

library(MCMCglmm)
library(ggplot2)
library(reshape2)
library(tidyverse)
library(performance)
library(lme4)
library(nlme)
library(VCVglmm)
library(multilevel)
library(arm)
theme_set(theme_bw())


# load data -----

pred_dat <- read.csv("Table_DO_Lynx_Bobcat_2017_Oct1.csv", header=TRUE, sep=",") %>% 
  rename_with(tolower) %>% 
  rename(region = studyregion) %>% 
  mutate(region = factor(region), year = factor(year), species = factor(species)) %>% 
  dplyr::select (-envelope)

head(pred_dat)
str(pred_dat)
summary(pred_dat)

qplot(d13c, d15n, data = pred_dat) + geom_smooth(method = "lm")

set.seed(102)

# fit model -----

mglmm_mod <- MCMCglmm(cbind(d13c, d15n) ~ trait*species*region -1, 
                      random = ~year, 
                      rcov = ~us(trait):units, 
                      data = pred_dat, 
                      family = c("gaussian", "gaussian"), 
                      burnin = 3000, nitt = 100000, thin = 10, 
                      pr = TRUE, pl = TRUE, 
                      verbose = F)


(mod_sum <- summary.MCMCglmm(mglmm_mod)) # the value for year is very close to zero suggesting little to no effect (https://stat.ethz.ch/pipermail/r-sig-mixed-models/2014q2/022068.html, https://rdrr.io/github/Euphrasiologist/VCVglmm/f/README.md)
summary(mglmm_mod) # summary.MCMCglmm() is a wrapper for summary()

plot(mglmm_mod$Sol) 
plot(mglmm_mod$VCV) 

# estimate marginal and conditional R2 values

means_1 <- mod_sum$solutions
means_2 <- as.data.frame(mglmm_mod$Sol) %>% 
  colMeans()
# solutions are the fixed effects, which are posterior means

str(means_1)
str(means_2)
mFixed <- means_1[,1] * t(mglmm_mod$X)
str(as.matrix(mFixed))
head(as.matrix(mFixed))
tail(as.matrix(mFixed))

# omit the first two rows and the sum is the total of the variance explained by fixed effects
# do this separately for d13C and d15N (like Kelley et al. 2020) 

index_d15N <- grep('traitd15n', row.names(mFixed[-c(1:2),])) # get positions of rows that start with 'traitd15n' excluding rows 1 and 2
d15n_fixed <- mFixed[2+(index_d15N),] # take all rows starting with 'traitd15n'; we need to add 2 to the index because we omitted the first two rows in the index to get fixed effects for d15n
d13c_fixed <- mFixed[-c(1:2,2+(index_d15N)),] # take all of the rows that don't start with 'traitd15n' and omit first 2 rows to get fixed effects for d13c
sum_d15n_fixed <- colSums(d15n_fixed)
sum_d13c_fixed <- colSums(d13c_fixed)
var_d15n_fixed <- var(sum_d15n_fixed)
var_d13c_fixed <- var(sum_d13c_fixed)
total_fixed_var <- sum(var_d13c_fixed, var_d15n_fixed)

# variance attributed to year (random effects variance)
year_random_var <- mod_sum$Gcovariances[,'post.mean']

# variance attributed to residuals
rcov_resid_var <- sum(mod_sum$Rcovariances[,'post.mean'])

# percent fixed
fixed_percent <- total_fixed_var/(total_fixed_var + year_random_var + rcov_resid_var)*100

# percent random
random_percent <- year_random_var/(total_fixed_var + year_random_var + rcov_resid_var)*100

# check DIC for model without random effect
mglmm_mod2 <- MCMCglmm(cbind(d13c, d15n) ~ trait*species*region -1, 
                      random = NULL, 
                      rcov = ~us(trait):units, 
                      data = pred_dat, 
                      family = c("gaussian", "gaussian"), 
                      burnin = 3000, nitt = 100000, thin = 10, 
                      pr = TRUE, pl = TRUE, 
                      verbose = F)

mod_sum2 <- summary(mglmm_mod2)


 