# Clean up of Christa's code for lynx-bobcat paper ('Code_Lynx-Bobcat_Diet_Overlap.txt')
# Author: Jenilee Gobin
# Date: 15 June, 2022

# Objective #3 - run analysis with and without outliers (nicheROVER)

# 1 - CHECK/LOAD PACKAGES ----- 

# function to look for & install multiple packages
check.packages <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

#  vector of needed packages
analysis.packages <- c("tidyverse", 
                       "viridis", 
                       "data.table", 
                       "tibble",
                       "ggpubr", 
                       "cowplot", 
                       "ggplot2", 
                       "dplyr", 
                       "SIBER", 
                       "magrittr", 
                       "nicheROVER", 
                       "cowplot",
                       "ggpattern",
                       "grid")

# apply function to packages
check.packages(analysis.packages)

set.seed(2)

# 2 - nicheROVER DEMO -----

# analysis for fish data

# formatting data for use in nicheROVER
data(fish) # 4 fish, 3 isotopes
head(fish)
aggregate(fish[2:4], fish[1], mean) # isotope means calculated for each species

# generate the posterior distributions of μ (mean) and Σ (variance) for isotope values for each species with the default prior
# this step is not absolutely necessary in generating the niche region and overlap plots, but can be useful during exploratory data analyses

# fish data
data(fish)

# generate parameter draws from the "default" posteriors of each fish
nsamples <- 1e3
system.time({
  fish.par <- tapply(1:nrow(fish), fish$species,
                     function(ii) niw.post(nsamples = nsamples, X = fish[ii,2:4]))
})
##    user  system elapsed 
##   0.143   0.017   0.180

# various parameter plots
clrs <- c("black", "red", "blue", "orange") # colors for each species

# mu1 (del15N), mu2 (del13C), and Sigma12
par(mar = c(4, 4, .5, .1)+.1, mfrow = c(1,3))
niche.par.plot(fish.par, col = clrs, plot.index = 1)
niche.par.plot(fish.par, col = clrs, plot.index = 2)
niche.par.plot(fish.par, col = clrs, plot.index = 1:2)
legend("topleft", legend = names(fish.par), fill = clrs)

# all mu (del15N, del13C, del34S)
niche.par.plot(fish.par, col = clrs, plot.mu = TRUE, plot.Sigma = FALSE)
legend("topleft", legend = names(fish.par), fill = clrs)

# all mu and Sigma
par(mar = c(4.2, 4.2, 2, 1)+.1)
niche.par.plot(fish.par, col = clrs, plot.mu = TRUE, plot.Sigma = TRUE)
legend("topright", legend = names(fish.par), fill = clrs)

# Create 2-d projections of a subset of niche regions
# 2-d projections of 10 niche regions
clrs <- c("black", "red", "blue", "orange") # colors for each species
nsamples <- 10
fish.par <- tapply(1:nrow(fish), fish$species,
                   function(ii) niw.post(nsamples = nsamples, X = fish[ii,2:4]))

# format data for plotting function
fish.data <- tapply(1:nrow(fish), fish$species, function(ii) X = fish[ii,2:4])

niche.plot(niche.par = fish.par, niche.data = fish.data, pfrac = .05,
           iso.names = expression(delta^{15}*N, delta^{13}*C, delta^{34}*S),
           col = clrs, xlab = expression("Isotope Ratio (per mil)"))

#  Calculate and display niche overlap estimates
# We use the function overlap() to calculate overlap metric estimates from a specified number of Monte Carlo draws (nsamples) from the fish.par parameter list. 
# It is necessary to specify the α-level. In this case, we have decided to calculate the overlap metric at two niche regions sizes for comparison: alpha=0.95 and alpha=0.99, or 95% and 99%.
# Then, we calculate the mean overlap metric between each species. Remember that the overlap metric is directional, such that it represents the probability that an individual from Species A will be found in the niche of Species B.

# niche overlap plots for 95% niche region sizes
nsamples <- 1000
fish.par <- tapply(1:nrow(fish), fish$species,
                   function(ii) niw.post(nsamples = nsamples, X = fish[ii,2:4]))

# Overlap calculation.  use nsamples = nprob = 10000 (1e4) for higher accuracy.
# the variable over.stat can be supplied directly to the overlap.plot function

over.stat <- overlap(fish.par, nreps = nsamples, nprob = 1e3, alpha = c(.95, 0.99))

#The mean overlap metrics calculated across iteratations for both niche 
#region sizes (alpha = .95 and alpha = .99) can be calculated and displayed in an array.
over.mean <- apply(over.stat, c(1:2,4), mean)*100
round(over.mean, 2)

## , , alpha = 95%
## 
##          Species B
## Species A  ARCS  BDWF  LKWF  LSCS
##      ARCS    NA 11.25 66.21 81.63
##      BDWF  0.31    NA 25.96  4.52
##      LKWF  7.33 78.17    NA 52.86
##      LSCS 37.58 51.79 88.55    NA
## 
## , , alpha = 99%
## 
##          Species B
## Species A  ARCS  BDWF  LKWF  LSCS
##      ARCS    NA 33.32 87.29 92.08
##      BDWF  0.79    NA 41.46  8.94
##      LKWF 11.67 92.19    NA 70.23
##      LSCS 50.05 80.14 97.05    NA

over.cred <- apply(over.stat*100, c(1:2, 4), quantile, prob = c(.025, .975), na.rm = TRUE)
round(over.cred[,,,1]) # display alpha = .95 niche region

#In the returned plot, Species A is along the rows and Species B is along columns. 
# The plots represent the posterior probability that an individual from the species indicated by the row will be found within the niche of the species indicated by the column header. 
# Before you plot, you must decide upon your α-level, and make sure the variable over.stat reflects this choice of α.

# Overlap plot.Before you run this, make sure that you have chosen your alpha level.
clrs <- c("black", "red", "blue", "orange") # colors for each species
over.stat <- overlap(fish.par, nreps = nsamples, nprob = 1e3, alpha = .95)
overlap.plot(over.stat, col = clrs, mean.cred.col = "turquoise", equal.axis = TRUE,
             xlab = "Overlap Probability (%) -- Niche Region Size: 95%")

# Calculate and display niche size estimates

# See ?niche.size for exactly how niche size is defined as a function of the parameters μ and Σ. 
# In a Bayesian context, we calculate the posterior distribution of niche size by species. 
# This done by calculating the niche size for every posterior sample of μ and Σ.

# posterior distribution of (mu, Sigma) for each species
nsamples <- 1000
fish.par <- tapply(1:nrow(fish), fish$species,
                   function(ii) niw.post(nsamples = nsamples, X = fish[ii,2:4]))

# posterior distribution of niche size by species
fish.size <- sapply(fish.par, function(spec) {
  apply(spec$Sigma, 3, niche.size, alpha = .95)
})

# point estimate and standard error
rbind(est = colMeans(fish.size),
      se = apply(fish.size, 2, sd))

# boxplots
clrs <- c("black", "red", "blue", "orange") # colors for each species
boxplot(fish.size, col = clrs, pch = 16, cex = .5,
        ylab = "Niche Size", xlab = "Species")

# 3 - nichROVER ANALYSIS FOR LYNX AND BOBCAT -----

# formatting data for use in nicheROVER
pred_dat <- read.csv("Table_DO_Lynx_Bobcat_2017_Oct1_outliers.csv", header=TRUE, sep=",") %>% 
  rename_with(tolower) %>% 
  rename(region = studyregion) %>% 
  mutate(region = replace(region, region == 'NWON', 'ON_west'), 
         region = replace(region, region == 'SSM', 'ON_central'), 
         region = factor(region, levels = c('BC', 'ON_west', 'ON_central', 'QC')), year = factor(year), species = factor(species))

head(pred_dat)
str(pred_dat)
summary(pred_dat)

# 3A - PARAMETER ESTIMATES -----

# generate the posterior distributions of μ (mean) and Σ (variance) for isotope values for each species with the default prior
# this step is not absolutely necessary in generating the niche region and overlap plots, but can be useful during exploratory data analyses

# generate parameter draws from the "default" posteriors of each species and region
nsamples <- 1e3
system.time({
  pred.par <- tapply(1:nrow(pred_dat), list(pred_dat$species, pred_dat$region),
                     function(ii) niw.post(nsamples = nsamples, X = pred_dat[ii,5:6]))
})

# the above code is measuring the system time used to apply the function niw.post, which generate random draws from p(u, sigma|x) for the Normal-Inverse-Wishart (NIW) prior, for each group of data

# various parameter plots
clrs <- palette(viridis(8)) # colors for each species

# mu1 (del15N), mu2 (del13C), and Sigma12
par(mar = c(4, 4, .5, .1)+.1, mfrow = c(1,3))
niche.par.plot(pred.par, col = clrs, plot.index = 1)
niche.par.plot(pred.par, col = clrs, plot.index = 2)
niche.par.plot(pred.par, col = clrs, plot.index = 1:2)
legend("topleft", legend = names(pred.par), fill = clrs)

# to customize plotting with ggplot, need to convert pred.par to a dataframe

str(pred.par)
head(pred.par)

str(pred.par[[1]], list.len = 8)

pred.par[[1]]$mu

pred.par[[1]]$mu[1,1:2]

df <- data.frame(pred.par)
pred_df <- data.frame(t(df))

str(pred_df)
lynx_tib <- as_tibble(pred_df$Lynx_canadensis)
str(lynx_tib)

lynx_BC_tib <- as_tibble(lynx_tib$BC)

# will probably be easier to split the data before running niw.post and then recombine it...

# generate parameter draws from the "default" posteriors of each species and region
nsamples <- 1e3

# split data by region and then group analysis by species

pred_dat_split <- pred_dat %>% 
  select(-c(year, envelope)) %>% 
  group_by(region) %>% 
  group_split()
  
str(pred_dat_split)
  
pred_dat_BC <-  pred_dat_split[[1]] 
pred_dat_ON_west <-  pred_dat_split[[2]] 
pred_dat_ON_central<- pred_dat_split[[3]]
pred_dat_QC<- pred_dat_split[[4]]

BC.par <- tapply(1:nrow(pred_dat_BC), pred_dat_BC$species,
                     function(ii) niw.post(nsamples = nsamples, X = pred_dat_BC[ii,3:4]))

ON_west.par <- tapply(1:nrow(pred_dat_ON_west), pred_dat_ON_west$species,
                 function(ii) niw.post(nsamples = nsamples, X = pred_dat_ON_west[ii,3:4]))

ON_central.par <- tapply(1:nrow(pred_dat_ON_central), pred_dat_ON_central$species,
                      function(ii) niw.post(nsamples = nsamples, X = pred_dat_ON_central[ii,3:4]))

QC.par <- tapply(1:nrow(pred_dat_QC), pred_dat_QC$species,
                         function(ii) niw.post(nsamples = nsamples, X = pred_dat_QC[ii,3:4]))

# various parameter plots
clrs <- palette(viridis(2)) # colors for each species

# mu1 (del15N), mu2 (del13C), and Sigma12
par(mar = c(4, 4, .5, .1)+.1, mfrow = c(3,1))
niche.par.plot(BC.par, col = clrs, plot.index = 1, ylab = expression('P('*{theta}*'|X)'))
niche.par.plot(BC.par, col = clrs, plot.index = 2, ylab = expression('P('*{theta}*'|X)'))
niche.par.plot(BC.par, col = clrs, plot.index = 1:2, ylab = expression('P('*{theta}*'|X)'))
legend("topright", legend = names(BC.par), fill = clrs)

str(BC.par)
str(BC.par$Lynx_canadensis)
str(BC.par$Lynx_canadensis$mu)
str(BC.par$Lynx_canadensis$Sigma)

str(BC.par$Lynx_canadensis$Sigma[,,1]) # comprises 1000 covariance matrices
BC.par$Lynx_canadensis$Sigma[,,1] # this is the 1st one; I need one value from each matrix - d13c,d15n (or the other way around, which is the same value)
BC.par$Lynx_canadensis$Sigma[,,1][1,2] # this is the position of the value I need from each covariance matrix in the list
BC.par$Lynx_canadensis$Sigma[,,1][2,1] # this is the alternative position of the value I need from each covariance matrix in the list

# one way to do this would be to loop through each list element (i.e., each covariance matrix) and extract that value
# I wonder if there is an easier alternative though

BC_lynx_sigma <- as.data.frame(BC.par$Lynx_canadensis$Sigma)
str(BC_lynx_sigma)
# if I make it a dataframe then each row of each covariance matrix becomes a variable (2 per matrix * 1000 matrices = 2000 variables in total)

BC_lynx_sigma_t <- as.data.frame(t(BC_lynx_sigma))
str(BC_lynx_sigma_t) # if we transpose this we end up with 2000 observations of 2 variables


BC_lynx_sigma_reduced_d15n <- BC_lynx_sigma %>% 
  select(starts_with('d15n.')) # take only the second row of each covariance matrix, meaning the first value is the covariance and second value is variance of d15n
# tried pivoting to longer and then grouping by isotope before transposing but names complicate this approach


BC_lynx_sigma_reduced_d13c <- BC_lynx_sigma %>% 
  select(starts_with('d13c.')) # take only the second row of each covariance matrix, meaning the second value is the covariance and the first is the variance of d13c.
# tried pivoting to longer and then grouping by isotope before transposing but names complicate this approach


str(BC_lynx_sigma_reduced_d15n)
str(BC_lynx_sigma_reduced_d13c)

BC_lynx_sigma_reduced_t_d15n <- BC_lynx_sigma_reduced_d15n %>% 
  t() %>% # once transposed, we the first variable is the covariance and the second is the variance of d15n
  as.data.frame() %>% 
  rename(covariance = d13c, var_n = d15n)

str(BC_lynx_sigma_reduced_t_d15n)

BC_lynx_sigma_reduced_t_d13c <- BC_lynx_sigma_reduced_d13c %>% 
  t() %>% # once transposed, we the second variable is the covariance and the second is the variance of d13c
  as.data.frame() %>% 
  rename(covariance = d15n, var_c = d13c)

str(BC_lynx_sigma_reduced_t_d13c)

BC_lynx_vcv <- BC_lynx_sigma_reduced_t_d15n %>% 
  mutate(var_c = BC_lynx_sigma_reduced_t_d13c$var_c, 
         rep = as.integer(gsub('d15n.', '', row.names(.))))

# develop for loop...

str(BC.par$Lynx_canadensis)

sigma_df_0 <- as.data.frame(BC.par$Lynx_canadensis)
str(sigma_df_0)

# isolate first row of covariance matrices
sigma_df_1 <- as.data.frame(BC.par$Lynx_canadensis) %>% 
  select(starts_with('Sigma.d13c'))  # note rows end up just being repeated and we only need the first (var_c) and second (covariance)

sigma_df_1_long <- sigma_df_1[1:2,] %>% 
  t() %>% 
  as.data.frame() %>% 
  rename(var_c = `1`, covariance = `2`) %>% 
  mutate(rep = as.integer(gsub('Sigma.d13c.', '', row.names(.))))

sigma_df_2 <- as.data.frame(BC.par$Lynx_canadensis) %>% 
  select(starts_with('Sigma.d15n'))

sigma_df_2_long <- sigma_df_2[1:2,] %>% 
  t() %>% 
  as.data.frame() %>% 
  rename(covariance= `1`, var_n = `2`) %>% 
  mutate(rep = as.integer(gsub('Sigma.d15n.', '', row.names(.))))

sigma_df_long <- sigma_df_1_long %>% 
  merge(sigma_df_2_long)

str(sigma_df_long)

# generate a function to apply this across lists of covariance matrices for all species and regions

get_params <- function(region_species_list, species, region){
  
  # isolate carbon variance and covariance
  sigma_df1 <- as.data.frame(region_species_list) %>% 
    select(starts_with('Sigma.d13c')) 
  
  # convert to long and rename variables
  sigma_df_1long <- sigma_df1[1:2,] %>% 
    t() %>% 
    as.data.frame() %>% 
    rename(var_c = `1`, covariance = `2`) %>% 
    mutate(rep = as.integer(gsub('Sigma.d13c.', '', row.names(.))))
  
  # isolate nitrogen variance and covariance
  sigma_df2 <- as.data.frame(region_species_list) %>% 
    select(starts_with('Sigma.d15n')) 
  
  # convert to long and rename variables
  sigma_df_2long <- sigma_df2[1:2,] %>% 
    t() %>% 
    as.data.frame() %>% 
    rename(covariance = `1`, var_n = `2`) %>% 
    mutate(rep = as.integer(gsub('Sigma.d15n.', '', row.names(.))))
  
  # combine all variances/covariances
  sigma_df_long <- sigma_df_1long %>%
    merge(sigma_df_2long)
  
  # isolate mu values and combine with covariances
  df <- as.data.frame(region_species_list) %>%
    select(starts_with('mu')) %>%
    cbind(sigma_df_long) %>%
    mutate(species = species,
           region = region)
  
  return((df))
  
}

get_params(BC.par$Lynx_canadensis, 'Lynx_canadensis', 'BC') %>% arrange(rep)

str(get_params(BC.par$Lynx_canadensis, 'Lynx_canadensis', 'BC'))
head(get_params(BC.par$Lynx_canadensis, 'Lynx_canadensis', 'BC'),5)
tail(get_params(BC.par$Lynx_canadensis, 'Lynx_canadensis', 'BC'),25)

# apply get_dat function to each group of data (this portion could be improved with  automation but don't have time to develop that)
BC_lynx <- get_params(BC.par$Lynx_canadensis, 'Lynx_canadensis', 'BC')
ON_west_lynx <-get_params(ON_west.par$Lynx_canadensis, 'Lynx_canadensis', 'ON_west')
ON_central_lynx <-get_params(ON_central.par$Lynx_canadensis, 'Lynx_canadensis', 'ON_central')
QC_lynx <- get_params(QC.par$Lynx_canadensis, 'Lynx_canadensis', 'QC')
BC_bobcat <- get_params(BC.par$Lynx_rufus, 'Lynx_rufus', 'BC')
ON_west_bobcat <- get_params(ON_west.par$Lynx_rufus, 'Lynx_rufus', 'ON_west')
ON_central_bobcat <- get_params(ON_central.par$Lynx_rufus, 'Lynx_rufus', 'ON_central')
QC_bobcat <- get_params(QC.par$Lynx_rufus, 'Lynx_rufus', 'QC')

str(BC_lynx)
str(ON_west_lynx)
str(ON_central_lynx)
str(QC_lynx)
str(BC_bobcat)
str(ON_west_bobcat)
str(ON_central_bobcat)
str(QC_bobcat)

all_param <- rbind(BC_lynx, ON_west_lynx, ON_central_lynx, QC_lynx, 
                   BC_bobcat, ON_west_bobcat, ON_central_bobcat, QC_bobcat) %>% 
  dplyr::rename(d13c = mu.d13c, d15n = mu.d15n) %>% 
  pivot_longer(cols = c(d13c, d15n, covariance, var_c, var_n), names_to = 'parameter', values_to = 'value') %>% 
  mutate(parameter2 = factor(parameter, levels = c('d13c', 'd15n', 'var_c', 'var_n', 'covariance'), labels = c(expression({delta}^13*C~' (\u2030)'), 
                                                                                                               expression({delta}^15*N~' (\u2030)'), 
                                                                                                               expression({Sigma}[{delta}^13*C]~ ('\u2030' ^2)), 
                                                                                                               expression({Sigma}[{delta}^15*N]~ ('\u2030' ^2)), 
                                                                                                               expression({Sigma}[{delta}^13*C*{delta}^15*N]~ ('\u2030' ^2))
  )
  ), 
  region = factor(region, levels = c('BC', 'ON_west', 'ON_central', 'QC')))

str(all_param)

# parameter distribution plot

ann_text <- data.frame(value = Inf, 
                       density = c(3, 14, 9, 14.5), 
                       parameter = 'covariance', 
                       region = c('BC', 'ON_west', 'ON_central', 'QC'),
                       label = c('BC', 'ON_west', 'ON_central', 'QC'), 
                       species = 'Lynx_canadensis') %>% 
  mutate(parameter2 = factor(parameter, levels = c('d13c', 'd15n', 'var_c', 'var_n', 'covariance'), labels = c(expression({delta}^13*C~' (\u2030)'), 
                                                                                                               expression({delta}^15*N~' (\u2030)'), 
                                                                                                               expression({Sigma}[{delta}^13*C]~ ( '\u2030' ^2)), 
                                                                                                               expression({Sigma}[{delta}^15*N]~ ( '\u2030' ^2)), 
                                                                                                               expression({Sigma}[{delta}^13*C*{delta}^15*N]~ ('\u2030' ^2))
  )
  ), 
  region = factor(region, levels = c('BC', 'ON_west', 'ON_central', 'QC')))


(param_dist_plot <- ggplot(all_param, aes(value, color = species, fill = species))
  + geom_density_pattern(aes(pattern_angle = species, pattern_color = species, pattern_fill = species, pattern_density = species))
  + facet_grid(region~parameter2, scales = 'free', switch = 'x', labeller = label_parsed)
  + scale_color_viridis_d(direction = -1, labels = c('Canada lynx', 'bobcat'))
  + scale_fill_viridis_d(direction = -1, alpha = 0.1, labels = c('Canada lynx', 'bobcat'))
  + scale_pattern_colour_viridis_d(direction = -1, labels = c('Canada lynx', 'bobcat'))
  + scale_pattern_fill_viridis_d(direction = -1, labels = c('Canada lynx', 'bobcat'))
  + scale_pattern_density_discrete(range = c(0.5, 0.01), labels = c('Canada lynx', 'bobcat'))
  + scale_pattern_angle_discrete(range = c(-45, 45), labels = c('Canada lynx', 'bobcat'))
  + scale_x_continuous(name = '')
  + scale_y_continuous(name = expression('P('*{theta}*'|X)'))
  + theme_classic2()
  + theme(legend.position = 'top', legend.text = element_text(size = 8),  
          legend.box.margin = margin(0, 0, 0, 0), legend.margin = margin(0, 0, 0, 0), legend.title = element_blank())
  + theme(strip.placement = 'outside', strip.background = element_blank(), strip.text.y = element_blank(), strip.text.x = element_text(size = 12))
  + theme(axis.title = element_text(size = 12), axis.text = element_text(size = 10))
  + geom_text(aes(x = value ,y = density), data = ann_text, label = ann_text$label, col = 'black', fontface = 'bold', size = 4, hjust = 1)
)

ggsave('Appendix_Fig3_param_plot_wOutliers.tiff', param_dist_plot, dpi = 500, width = 7, height = 5, units = 'in')

# 3B - NICHE OVERLAP -----

# mean probability of overlap in 95% niche region
BC_over.stat <- overlap(BC.par, nreps = nsamples, nprob = 1e3, alpha = 0.95)
ON_west_over.stat <- overlap(ON_west.par, nreps = nsamples, nprob = 1e3, alpha = 0.95)
ON_central_over.stat <- overlap(ON_central.par, nreps = nsamples, nprob = 1e3, alpha = 0.95)
QC_over.stat <- overlap(QC.par, nreps = nsamples, nprob = 1e3, alpha = 0.95)

str(BC_over.stat)
str(over.stat)

#Remember that the overlap metric is directional, such that it represents the probability that an individual from Species A will be found in the niche of Species B.
BC_over.mean <- apply(BC_over.stat, c(1:1,2), mean)*100
round(BC_over.mean, 2)

ON_west_over.mean <- apply(ON_west_over.stat, c(1:1,2), mean)*100
round(ON_west_over.mean, 2)

ON_central_over.mean <- apply(ON_central_over.stat, c(1:1,2), mean)*100
round(ON_central_over.mean, 2)

QC_over.mean <- apply(QC_over.stat, c(1:1,2), mean)*100
round(QC_over.mean, 2)

# credible intervals
BC_over.cred <- apply(BC_over.stat*100, c(1:1, 2), quantile, prob = c(.025, .975), na.rm = TRUE)
round(BC_over.cred[,,1:2]) # display alpha = .95 niche region

ON_west_over.cred <- apply(ON_west_over.stat*100, c(1:1, 2), quantile, prob = c(.025, .975), na.rm = TRUE)
round(ON_west_over.cred[,,1:2]) # display alpha = .95 niche region

ON_central_over.cred <- apply(ON_central_over.stat*100, c(1:1, 2), quantile, prob = c(.025, .975), na.rm = TRUE)
round(ON_central_over.cred[,,1:2]) # display alpha = .95 niche region

QC_over.cred <- apply(QC_over.stat*100, c(1:1, 2), quantile, prob = c(.025, .975), na.rm = TRUE)
round(QC_over.cred[,,1:2]) # display alpha = .95 niche region

# Overlap plot.Before you run this, make sure that you have chosen your alpha level.
clrs <- c("black", "red") # colors for each species

BC_over.stat <- overlap(BC.par, nreps = nsamples, nprob = 1e3, alpha = .95)
overlap.plot(BC_over.stat, col = clrs, mean.cred.col = "turquoise", equal.axis = TRUE,
             xlab = "Overlap Probability (%) -- Niche Region Size: 95%")

ON_west_over.stat <- overlap(ON_west.par, nreps = nsamples, nprob = 1e3, alpha = .95)
overlap.plot(ON_west_over.stat, col = clrs, mean.cred.col = "turquoise", equal.axis = TRUE,
             xlab = "Overlap Probability (%) -- Niche Region Size: 95%")

ON_central_over.stat <- overlap(ON_central.par, nreps = nsamples, nprob = 1e3, alpha = .95)
overlap.plot(ON_central_over.stat, col = clrs, mean.cred.col = "turquoise", equal.axis = TRUE,
             xlab = "Overlap Probability (%) -- Niche Region Size: 95%")

QC_over.stat <- overlap(QC.par, nreps = nsamples, nprob = 1e3, alpha = .95)
overlap.plot(QC_over.stat, col = clrs, mean.cred.col = "turquoise", equal.axis = TRUE,
             xlab = "Overlap Probability (%) -- Niche Region Size: 95%")

# plot with ggplot ----------------------------------

#Remember that the overlap metric is directional, such that it represents the probability that an individual from Species A will be found in the niche of Species B.

# create function to get overlap metrics for plotting in ggplot

get_overlap <- function (mean, cred, stat, region){
  
  mean_df <- as.data.frame(mean) %>% 
    pivot_longer(everything(), names_to = 'species_B', values_to = 'mean') %>% 
    filter(!is.na(mean))
  
  cred_df <- as.data.frame(cred) %>%
     t() %>%
     as.data.frame() %>%
     mutate(species_A = str_split_fixed(row.names(.), '[.]', n=2)[,1],
            species_B = str_split_fixed(row.names(.), '[.]', n=2)[,2])

   stat_df <- as.data.frame(stat) %>%
     t() %>%
     as.data.frame() %>%
     mutate(species_B = gsub('[.,0-9]', "", row.names(.))) %>%
     pivot_longer(cols = 1:2, names_to = 'species_A', values_to = 'overlap') %>%
     filter(!is.na(overlap)) %>%
     merge(cred_df) %>%
     merge(mean_df) %>% 
     mutate(region = region)

  return(stat_df)
  
}

BC_overlap <- get_overlap(BC_over.mean, BC_over.cred, BC_over.stat, 'BC')
ON_west_overlap <- get_overlap(ON_west_over.mean, ON_west_over.cred, ON_west_over.stat, 'ON_west')
ON_central_overlap <- get_overlap(ON_central_over.mean, ON_central_over.cred, ON_central_over.stat, 'ON_central')
QC_overlap <- get_overlap(QC_over.mean, QC_over.cred, QC_over.stat, 'QC')

all_overlap <- rbind(BC_overlap, ON_west_overlap, ON_central_overlap, QC_overlap) %>% 
  mutate(region = factor(region, levels = c('BC', 'ON_west', 'ON_central', 'QC')))
str(all_overlap)

ann_text2 <- data.frame(overlap = 0, 
                       density = 0.125, 
                       parameter = 'covariance', 
                       region = c('BC', 'ON_west', 'ON_central', 'QC'),
                       label = c('BC', 'ON_west', 'ON_central', 'QC'), 
                       species_A = 'Lynx_canadensis') %>% 
  mutate(parameter = factor(parameter, levels = c('d13c', 'd15n', 'covariance'), labels = c(expression({delta}^13*C~' (\u2030)'), expression({delta}^15*N~' (\u2030)'), expression({Sigma}[{delta}^13*C*{delta}^15*N]))), 
         region = factor(region, levels = c('BC', 'ON_west', 'ON_central', 'QC')))

(overlap_plot1 <- ggplot(all_overlap, aes(overlap*100, color = species_A, fill = species_A))
  + geom_vline(aes(xintercept = mean, color = species_A), lwd = 1)
  + geom_vline(aes(xintercept = `2.5%`, color = species_A), lty = 2, lwd = 0.9)
  + geom_vline(aes(xintercept = `97.5%`, color = species_A), lty = 2, lwd = 0.9)
  + geom_density_pattern(aes(pattern_angle = species_A, pattern_color = species_A, pattern_fill = species_A, pattern_density = species_A))
  + facet_grid(region~species_A)
  + scale_color_viridis_d(direction = -1, labels = c('Canada lynx', 'bobcat'))
  + scale_fill_viridis_d(direction = -1, alpha = 0.1, labels = c('Canada lynx', 'bobcat'))
  + scale_pattern_colour_viridis_d(direction = -1, labels = c('Canada lynx', 'bobcat'))
  + scale_pattern_fill_viridis_d(direction = -1, labels = c('Canada lynx', 'bobcat'))
  + scale_pattern_density_discrete(range = c(0.5, 0.01), labels = c('Canada lynx', 'bobcat'))
  + scale_pattern_angle_discrete(range = c(-45, 45), labels = c('Canada lynx', 'bobcat'))
  + scale_x_continuous(name = 'niche overlap (%)')
  + scale_y_continuous(name = expression('P('*{theta}*'|X)'))
  + theme_classic2()
  + theme(legend.position = 'top', legend.text = element_text(size = 8),  
          legend.box.margin = margin(0, 0, 0, 0), legend.margin = margin(0, 0, 0, 0), legend.title = element_blank())
  + theme(strip.placement = 'outside', strip.background = element_blank(), strip.text = element_blank())
  + theme(axis.title = element_text(size = 12), axis.text = element_text(size = 10))
  + geom_text(aes(x = overlap ,y = density), data = ann_text2, label = ann_text2$label, col = 'black', fontface = 'bold', size = 4, hjust = 0)
)

ggsave('Appendix_Fig4_overlap_plot_v1_wOutliers.tiff', overlap_plot1, dpi = 500, width = 6, height = 6, units = 'in')

(overlap_plot2 <- ggplot(all_overlap, aes(overlap*100, color = species_A, fill = species_A))
  + geom_vline(aes(xintercept = mean, color = species_A), lwd = 1)
  + geom_vline(aes(xintercept = `2.5%`, color = species_A), lty = 2, lwd = 0.9)
  + geom_vline(aes(xintercept = `97.5%`, color = species_A), lty = 2, lwd = 0.9)
  + geom_density_pattern(aes(pattern_angle = species_A, pattern_color = species_A, pattern_fill = species_A, pattern_density = species_A))
  + facet_grid(species_A~region, scales = 'free')
  + scale_color_viridis_d(direction = -1, labels = c('Canada lynx', 'bobcat'))
  + scale_fill_viridis_d(direction = -1, alpha = 0.1, labels = c('Canada lynx', 'bobcat'))
  + scale_pattern_colour_viridis_d(direction = -1, labels = c('Canada lynx', 'bobcat'))
  + scale_pattern_fill_viridis_d(direction = -1, labels = c('Canada lynx', 'bobcat'))
  + scale_pattern_density_discrete(range = c(0.5, 0.01), labels = c('Canada lynx', 'bobcat'))
  + scale_pattern_angle_discrete(range = c(-45, 45), labels = c('Canada lynx', 'bobcat'))
  + scale_x_continuous(name = 'Probability of niche overlap with 95% niche region (%)')
  + scale_y_continuous(name = expression('P('*{theta}*'|X)'))
  + theme_classic2()
  + theme(legend.position = 'top', legend.text = element_text(size = 8),  
          legend.box.margin = margin(0, 0, 0, 0), legend.margin = margin(0, 0, 0, 0), legend.title = element_blank())
  + theme(strip.placement = 'outside', strip.background.y = element_blank(), strip.text.y = element_blank())
  + theme(axis.title = element_text(size = 12), axis.text = element_text(size = 10))
  #+ geom_text(aes(x = overlap ,y = density), data = ann_text2, label = ann_text2$label, col = 'black', fontface = 'bold', size = 4, hjust = 0)
)

ggsave('Appendix_Fig4_overlap_plot_v2_wOutliers.tiff', overlap_plot2, dpi = 500, width = 7, height = 4, units = 'in')


# Calculate and display niche size estimates ----------------------------------------
# posterior distribution of niche size by species
BC.size <- sapply(BC.par, function(spec) {
  apply(spec$Sigma, 3, niche.size, alpha = .95)
})

ON_west.size <- sapply(ON_west.par, function(spec) {
  apply(spec$Sigma, 3, niche.size, alpha = .95)
})


ON_central.size <- sapply(ON_central.par, function(spec) {
  apply(spec$Sigma, 3, niche.size, alpha = .95)
})


QC.size <- sapply(QC.par, function(spec) {
  apply(spec$Sigma, 3, niche.size, alpha = .95)
})


# point estimate and standard error
rbind(est = colMeans(BC.size),
      se = apply(BC.size, 2, sd))

rbind(est = colMeans(ON_west.size),
      se = apply(ON_west.size, 2, sd))

rbind(est = colMeans(ON_central.size),
      se = apply(ON_central.size, 2, sd))

rbind(est = colMeans(QC.size),
      se = apply(QC.size, 2, sd))
  
# boxplots
par(mar = c(4, 4, .5, .1)+.1, mfrow = c(2,2))

clrs <- c("black", "red") # colors for each species

boxplot(BC.size, col = clrs, pch = 16, cex = .5,
        ylab = "Niche Size", xlab = "Species")

boxplot(ON_west.size, col = clrs, pch = 16, cex = .5,
        ylab = "Niche Size", xlab = "Species")

boxplot(ON_central.size, col = clrs, pch = 16, cex = .5,
        ylab = "Niche Size", xlab = "Species")

boxplot(QC.size, col = clrs, pch = 16, cex = .5,
        ylab = "Niche Size", xlab = "Species")

# results appear consistent with those from SIBER  
  
  
  

















