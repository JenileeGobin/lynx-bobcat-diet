# Clean up of Christa's code for lynx-bobcat paper ('Code_Lynx-Bobcat_Diet_Overlap.txt')
# Author: Jenilee Gobin
# Date: 13 June, 2022

# Objective #3 - run analysis with and without outliers (SIBER)

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
                       "ggpattern")

# apply function to packages
check.packages(analysis.packages)

set.seed(2)

# 2 - SIBER DEMO -----
# https://cran.r-project.org/web/packages/SIBER/vignettes/Introduction-to-SIBER.html
# 2A - create SIBER object -----

fname <- system.file("extdata", "demo.siber.data.csv", package = "SIBER")
str(fname)

demo_dat <- read.csv(fname, header=T)
head(demo_dat) # this is how stable istope data should be organized
# Community 1 comprises 3 groups and drawn as black, red and green circles; community 2 comprises 3 groups and drawn as black, red and green triangles

siber_demo <- createSiberObject(demo_dat)
str(siber_demo)

# 2B - plot SIBER object -----

# iso.order is a vector of length 2 specifying which isotope should be plotted on the x and y axes.
# N.B. there is currently a problem with the addition of the group ellipses using if you deviate from the default of iso.order = c(1,2).
# I recommend you set up your original data with the chemical element you want plotted on the x-axis being the first column, and the y-axis in the second.
# ***Therefore have C in first column and N in second column***

# Ellipses are drawn for each group independently with ellipses = T. 
# These ellipses can be made to be maximum likelihood standard ellipses by setting p = NULL, or can be made to be prediction ellipses that contain approximately p proportion of data. 
# For example, p = 0.95 will draw an ellipse that encompasses approximately 95% of the data. 
# The parameter n determines how many points are used to make each ellipse and hence how smooth the curves are.

# Create lists of plotting arguments to be passed onwards to each of the three plotting functions.
community.hulls.args <- list(col = 1, lty = 1, lwd = 1)
group.ellipses.args  <- list(n = 100, p.interval = 0.95, lty = 1, lwd = 2)
group.hulls.args     <- list(lty = 2, col = "grey20")

# create viridis palette for the 3 groups
palette(viridis(3))

par(mfrow=c(1,1))
plotSiberObject(siber_demo,
                ax.pad = 2, #ax.pad determines the padding applied around the extremes of the data
                hulls = F, community.hulls.args = community.hulls.args, 
                ellipses = T, group.ellipses.args = group.ellipses.args,
                group.hulls = T, group.hulls.args = group.hulls.args,
                bty = "L",
                iso.order = c(1,2),
                xlab = expression({delta}^13*C~'\u2030'),
                ylab = expression({delta}^15*N~'\u2030')
)

#We create the bi-plot again here and this time add the additional ellipses overlayed on the basic plot that this time omits group hulls and group standard ellipses.

par(mfrow=c(1,1))

community.hulls.args <- list(col = 1, lty = 1, lwd = 1)
group.ellipses.args  <- list(n = 100, p.interval = 0.95, lty = 1, lwd = 2)
group.hull.args      <- list(lty = 2, col = "grey20")

# this time we will make the points a bit smaller by 
# cex = 0.5
plotSiberObject(siber_demo,
                ax.pad = 2, 
                hulls = F, community.hulls.args, 
                ellipses = F, group.ellipses.args,
                group.hulls = F, group.hull.args,
                bty = "L",
                iso.order = c(1,2),
                xlab=expression({delta}^13*C~'\u2030'),
                ylab=expression({delta}^15*N~'\u2030'),
                cex = 0.5
)



# Calculate summary statistics for each group: TA, SEA and SEAc
group_ML <- groupMetricsML(siber_demo)
print(group_ML)
#>            1.1       1.2       1.3       2.1       2.2       2.3
#> TA   21.924922 10.917715 17.945127 3.0714363 11.476354 1.4818061
#> SEA   5.783417  3.254484  5.131601 0.8623300  3.458824 0.4430053
#> SEAc  5.989967  3.370715  5.314872 0.8931275  3.582354 0.4588269


# You can add more ellipses by directly calling plot.group.ellipses()
# Add an additional p.interval % prediction ellilpse
plotGroupEllipses(siber_demo, n = 100, p.interval = 0.95,
                  lty = 1, lwd = 2)

# or you can add the XX% confidence interval around the bivariate means
# by specifying ci.mean = T along with whatever p.interval you want.
plotGroupEllipses(siber_demo, n = 100, p.interval = 0.95, ci.mean = T,
                  lty = 1, lwd = 2)

#Alternatively, we may wish to focus on comparing the two communities represented in these plots by the open circles and the open triangles. To illustrate these groupings, we might draw the convex hull between the means of each of the three groups comprising each community. Additionally, I have highlighted the location of each group by adding the 95% confidence interval of their bivariate mean.

# A second plot provides information more suitable to comparing
# the two communities based on the community-level Layman metrics

# this time we will make the points a bit smaller by 
# cex = 0.5
plotSiberObject(siber_demo,
                ax.pad = 2, 
                hulls = T, community.hulls.args, 
                ellipses = F, group.ellipses.args,
                group.hulls = F, group.hull.args,
                bty = "L",
                iso.order = c(1,2),
                xlab=expression({delta}^13*C~'\u2030'),
                ylab=expression({delta}^15*N~'\u2030'),
                cex = 0.5
)

# or you can add the XX% confidence interval around the bivariate means
# by specifying ci.mean = T along with whatever p.interval you want.
plotGroupEllipses(siber_demo, n = 100, p.interval = 0.95,
                  ci.mean = T, lty = 1, lwd = 2) 

# Calculate the various Layman metrics on each of the communities.
community_ML <- communityMetricsML(siber_demo) 
print(community_ML)


#>                  1         2
#> dY_range  6.526633  4.617564
#> dX_range  8.184069 10.184171
#> TA       13.079701 12.799232
#> CD        4.235518  4.809830
#> MNND      4.944052  5.111483
#> SDNND     3.526088  4.598012
# dY range (aka NR) - range of d15N 
# dx range (aka CR) - range of d13C
# CD - mean distance to centroid
# MNND - mean nearest neighbour distance
# SDNND - standard deviation of nearest neighbout distance
# TA - total area 
# SEA_B - standard Bayesian ellipse area
# SEA_C - size-corrected standard ellipse area

# 2C - fitting Bayesian models to the data -----

# Whether your intended analysis is to compare isotopic niche width among groups, or among communities, the initial step is to fit Bayesian multivariate normal distributions to each group in the dataset. 
# The decision as to whether you then want to compare the area of the ellipses among groups, or any / all of the 6 Layman metrics comes later.

# These multivariate normal distributions are fitted using the jags software run via the package rjags. 
# This method relies on an iterated Gibbs Sampling technique and some information on the length, number and iterations of sampling chains is required. 
# Additionally, the prior distributions for the parameters need to be specified. 
# In SIBER, these are bundled into two list objects: parms which holds the parameters defining how the sampling algorithm is to run; and priors which holds information on the prior distributions of the parameters to be estimated. 
# Typically, the priors are left vague and you should use these same values in your own analysis. 
# Since the data are z-scored internally before the models are fitted to the data, the expected means are inherently close to zero, and the marginal variances close to one. 
# This greatly aids the jags fitting process.

# After calling siberMVN() you will see output in the command window indicating that the jags models are being fitted, one block of output for each group in your dataset. 
# A subset of these blocks are shown below.

# options for running jags
parms <- list()
parms$n.iter <- 2 * 10^4   # number of iterations to run the model for
parms$n.burnin <- 1 * 10^3 # discard the first set of values
parms$n.thin <- 10     # thin the posterior by this many
parms$n.chains <- 2        # run this many chains

# define the priors
priors <- list()
priors$R <- 1 * diag(2)
priors$k <- 2
priors$tau.mu <- 1.0E-3

# fit the ellipses which uses an Inverse Wishart prior
# on the covariance matrix Sigma, and a vague normal prior on the 
# means. Fitting is via the JAGS method.
ellipses_posterior <- siberMVN(siber_demo, parms, priors)
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 30
#>    Unobserved stochastic nodes: 3
#>    Total graph size: 45
#> 
#> Initializing model
#> 
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 30
#>    Unobserved stochastic nodes: 3
#>    Total graph size: 45
#> 
#> Initializing model
#> 
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 30
#>    Unobserved stochastic nodes: 3
#>    Total graph size: 45
#> 
#> Initializing model
#> 
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 30
#>    Unobserved stochastic nodes: 3
#>    Total graph size: 45
#> 
#> Initializing model
#> 
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 30
#>    Unobserved stochastic nodes: 3
#>    Total graph size: 45
#> 
#> Initializing model
#> 
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 30
#>    Unobserved stochastic nodes: 3
#>    Total graph size: 45
#> 
#> Initializing model

# Comparing among groups using the Standard Ellipse Area
# When comparing individual groups with each other, be it within a single community, or groups among communities, the Standard Ellipse Area (SEA) is the recommended method. 
# Since the multivariate normal distributions have already been fitted to each group, it only remains to calculate the SEA on the posterior distribution of covariance matrix for each group, thereby yielding the Bayesian SEA or SEA-B. 
# We can also use the summary statistics we calculated earlier to add the maximum likelihood estimates of SEA-c to the Bayesian estimates.

# Plotting is via the function siberDensityPlot() which is essentially the same as siardensityplot() from the older version of SIAR. 
# Credible intervals can be extracted by calling the function hdr from the hdrcde package.

# The posterior estimates of the ellipses for each group can be used to
# calculate the SEA.B for each group.
SEA_B <- siberEllipses(ellipses_posterior)

siberDensityPlot(SEA_B, xticklabels = colnames(group_ML), 
                 xlab = c("Community | Group"),
                 ylab = expression("Standard Ellipse Area " ('\u2030' ^2) ),
                 bty = "L",
                 las = 1,
                 main = "SIBER ellipses on each group"
)

# Add red x's for the ML estimated SEA-c
points(1:ncol(SEA_B), group_ML[3,], col="red", pch = "x", lwd = 2)

# Calculate some credible intervals 
cr_p <- c(0.95, 0.99) # vector of quantiles

# call to hdrcde:hdr using lapply()
SEA_B_credibles <- lapply(
  as.data.frame(SEA_B), 
  function(x,...){tmp<-hdrcde::hdr(x)$hdr},
  prob = cr_p)

# do similar to get the modes, taking care to pick up multimodal posterior
# distributions if present
SEA_B_modes <- lapply(
  as.data.frame(SEA_B), 
  function(x,...){tmp<-hdrcde::hdr(x)$mode},
  prob = cr_p, all.modes=T)


# 2D - plot posterior ellipses with ggplot -----
# https://cran.r-project.org/web/packages/SIBER/vignettes/Plot-posterior-ellipses.html

library(SIBER)
library(tidyverse)
library(ellipse)

# load in the included demonstration dataset
data("demo.siber.data")
#
# create the siber object
siber.example <- createSiberObject(demo.siber.data)

# Calculate summary statistics for each group: TA, SEA and SEAc
group.ML <- groupMetricsML(siber.example)

# options for running jags
parms <- list()
parms$n.iter <- 2 * 10^4   # number of iterations to run the model for
parms$n.burnin <- 1 * 10^3 # discard the first set of values
parms$n.thin <- 10     # thin the posterior by this many
parms$n.chains <- 2        # run this many chains

# define the priors
priors <- list()
priors$R <- 1 * diag(2)
priors$k <- 2
priors$tau.mu <- 1.0E-3

# fit the ellipses which uses an Inverse Wishart prior
# on the covariance matrix Sigma, and a vague normal prior on the 
# means. Fitting is via the JAGS method.
ellipses.posterior <- siberMVN(siber.example, parms, priors)

# The posterior estimates of the ellipses for each group can be used to
# calculate the SEA.B for each group.
SEA.B <- siberEllipses(ellipses.posterior)

siberDensityPlot(SEA.B, xticklabels = colnames(group.ML), 
                 xlab = c("Community | Group"),
                 ylab = expression("Standard Ellipse Area " ('\u2030' ^2) ),
                 bty = "L",
                 las = 1,
                 main = "SIBER ellipses on each group"
)

# Now we want to create some plots of some sample ellipses from these distributions. We need to create a data.frame object of all the ellipses for each group. In this example we simply take the first 10 posterior draws assuming them to be independent of one another, but you could take a random sample if you prefer.

# how many of the posterior draws do you want?
n.posts <- 10

# decide how big an ellipse you want to draw
#p.ell <- 0.95

# for a standard ellipse use
 p.ell <- pchisq(1,2)

# a list to store the results
all_ellipses <- list()

# loop over groups
for (i in 1:length(ellipses.posterior)){
  
  # a dummy variable to build in the loop
  ell <- NULL
  post.id <- NULL
  
  for ( j in 1:n.posts){
    
    # covariance matrix
    Sigma  <- matrix(ellipses.posterior[[i]][j,1:4], 2, 2)
    
    # mean
    mu     <- ellipses.posterior[[i]][j,5:6]
    
    # ellipse points
    
    out <- ellipse::ellipse(Sigma, centre = mu , level = p.ell)
    
    
    ell <- rbind(ell, out)
    post.id <- c(post.id, rep(j, nrow(out)))
    
  }
  ell <- as.data.frame(ell)
  ell$rep <- post.id
  all_ellipses[[i]] <- ell
}

ellipse_df <- bind_rows(all_ellipses, .id = "id")

# now we need the group and community names

# extract them from the ellipses.posterior list
group_comm_names <- names(ellipses.posterior)[as.numeric(ellipse_df$id)]

# split them and convert to a matrix, NB byrow = T
split_group_comm <- matrix(unlist(strsplit(group_comm_names, "[.]")),
                           nrow(ellipse_df), 2, byrow = TRUE)

ellipse_df$community <- split_group_comm[,1]
ellipse_df$group     <- split_group_comm[,2]

ellipse_df <- dplyr::rename(ellipse_df, iso1 = x, iso2 = y)

# Now to create the plots. First plot all the raw data as we want.

first.plot <- ggplot(data = demo.siber.data, aes(iso1, iso2)) +
  geom_point(aes(color = factor(group):factor(community)), size = 2)+
  ylab(expression(paste(delta^{15}, "N (\u2030)")))+
  xlab(expression(paste(delta^{13}, "C (\u2030)"))) + 
  theme(text = element_text(size=15))
print(first.plot)

# Now we can try to add the posterior ellipses on top and facet by group

second.plot <- first.plot + facet_wrap(~factor(group):factor(community))
print(second.plot)

# rename columns of ellipse_df to match the aesthetics
third.plot <- second.plot + 
  geom_polygon(data = ellipse_df,
               mapping = aes(iso1, iso2,
                             group = rep,
                             color = factor(group):factor(community),
                             fill = NULL),
               fill = NA,
               alpha = 0.2)
print(third.plot)




# 3 - SIBER ANALYSIS FOR LYNX AND BOBCAT -----
# 3A - create SIBER object -----
pred_dat <- read.csv("Table_DO_Lynx_Bobcat_2017_Oct1_outliers.csv", header=TRUE, sep=",") %>% 
  rename_with(tolower) %>% 
  dplyr::rename(region = studyregion) %>% 
  mutate(region = replace(region, region == 'NWON', 'ON_west'), 
         region = replace(region, region == 'SSM', 'ON_central'), 
         region = factor(region, levels = c('BC', 'ON_west', 'ON_central', 'QC')), year = factor(year), species = factor(species))

head(pred_dat)
str(pred_dat)
summary(pred_dat)

siber_pred_dat <- pred_dat %>% 
  relocate(d13c, .before = region) %>% 
  relocate(d15n, .after = d13c) %>% 
  select(-c(year, envelope)) %>% 
  dplyr::rename(iso1 = d13c, 
         iso2 = d15n, 
         group = region, 
         community = species)
  

siber_pred <- createSiberObject(siber_pred_dat)
str(siber_pred)

# 3B - plot SIBER object -----

# iso.order is a vector of length 2 specifying which isotope should be plotted on the x and y axes.
# N.B. there is currently a problem with the addition of the group ellipses using if you deviate from the default of iso.order = c(1,2).
# I recommend you set up your original data with the chemical element you want plotted on the x-axis being the first column, and the y-axis in the second.
# ***Therefore have C in first column and N in second column***

# Ellipses are drawn for each group independently with ellipses = T. 
# These ellipses can be made to be maximum likelihood standard ellipses by setting p = NULL, or can be made to be prediction ellipses that contain approximately p proportion of data. 
# For example, p = 0.95 will draw an ellipse that encompasses approximately 95% of the data. 
# The parameter n determines how many points are used to make each ellipse and hence how smooth the curves are.

# Create lists of plotting arguments to be passed onwards to each of the three plotting functions.
community.hulls.args <- list(col = 1, lty = 1, lwd = 1)
group.ellipses.args  <- list(n = 100, p.interval = 0.95, lty = 1, lwd = 2)
group.hulls.args     <- list(lty = 2, col = "grey20")


par(mfrow=c(1,1))
plotSiberObject(siber_pred,
                ax.pad = 2, #ax.pad determines the padding applied around the extremes of the data
                hulls = F, community.hulls.args = community.hulls.args, 
                ellipses = T, group.ellipses.args = group.ellipses.args,
                group.hulls = T, group.hulls.args = group.hulls.args,
                bty = "L",
                iso.order = c(1,2),
                xlab = expression({delta}^13*C~'\u2030'),
                ylab = expression({delta}^15*N~'\u2030')
)

# to create customized plots, probably best to do with ggplot:

# faceted plot of isotope ratios by species and region

(SI_plot_1 <- ggplot(pred_dat, aes(x = d13c, y = d15n))
  + stat_ellipse(aes(col = species, lty = species), type = 'norm', lwd = 0.8, level = 0.95) #SIBER and nicheROVER assume multivariate normal distribution
  + geom_point(aes(fill = species, pch = species), cex = 2)
  + facet_grid(.~region, scales = 'free')
  + scale_x_continuous(name = expression({delta}^13*C~' (\u2030)'))
  + scale_y_continuous(name = expression({delta}^15*N~' (\u2030)'))
  + scale_fill_viridis_d(direction = -1, labels = c('Canada lynx', 'bobcat'))
  + scale_colour_viridis_d(direction = -1, labels = c('Canada lynx', 'bobcat'))
  + scale_linetype_manual(values = c(1,2), labels = c('Canada lynx', 'bobcat'))
  + scale_shape_manual(values = c(21, 24), labels = c('Canada lynx', 'bobcat'))
  + theme_classic2()
  + theme(legend.position = 'top', legend.text = element_text(size = 8), legend.key.size = unit(0.75, 'in'), 
          legend.box.margin = margin(--0, 0, -10, 0), legend.margin = margin(-10, 0, -10, 0), legend.title = element_blank())
  + theme(axis.title = element_text(size = 12), axis.text = element_text(size = 10))
)

ggsave('Appendix_Fig1_SI_plot_wOutliers.tiff', SI_plot_1, dpi = 500, width = 7, height = 4, units = 'in')
# ellipses drawn with stat_ellipse contain 95% of data; could also use geom_polygon with ellipse coordinate data from siber to plot prediction ellipse containing 95% of data
# https://cran.r-project.org/web/packages/SIBER/vignettes/Customising-Plots-Manually.html
# https://cran.r-project.org/web/packages/SIBER/vignettes/Plot-posterior-ellipses.html

# 3C - fit Bayesian models to data -----

# options for running jags
parms <- list()
parms$n.iter <- 1 * 10^4   # number of iterations to run the model for
parms$n.burnin <- 1 * 10^3 # discard the first set of values
parms$n.thin <- 10     # thin the posterior by this many
parms$n.chains <- 2        # run this many chains

# define the priors
priors <- list()
priors$R <- 1 * diag(2) # {R}{The scaling vector for the diagonal of Inverse Wishart distribution prior on the covariance matrix Sigma. Typically set to a 2x2 matrix matrix(c(1, 0, 0, 1), 2, 2).}
priors$k <- 2 # {k}{The degrees of freedom of the Inverse Wishart distribution for the covariance matrix Sigma. Typically set to the dimensionality of Sigma, which in this bivariate case is 2.}
priors$tau.mu <- 1.0E-3 # {tau}{The precision on the normal prior on the means mu.}

# fit the ellipses which uses an Inverse Wishart prior
# on the covariance matrix Sigma, and a vague normal prior on the 
# means. Fitting is via the JAGS method.
pred_ellipses_posterior <- siberMVN(siber_pred, parms, priors)
str(pred_ellipses_posterior)

# The posterior estimates of the ellipses for each group can be used to
# calculate the SEA.B for each group.
pred_SEA_B <- siberEllipses(pred_ellipses_posterior)
str(pred_SEA_B)

# Calculate summary statistics for each group: TA, SEA and SEAc
pred_group_ML <- groupMetricsML(siber_pred)
print(pred_group_ML)


pred_SEA_B_plot <- siberDensityPlot(pred_SEA_B, xticklabels = colnames(pred_group_ML), 
                 xlab = c("Species | Region"),
                 ylab = expression("Standard Ellipse Area " ('\u2030' ^2) ),
                 bty = "L",
                 las = 1,
                 main = "SIBER ellipses on each group"
)

# Add red x's for the ML estimated SEA-c
points(1:ncol(pred_SEA_B), pred_group_ML[3,], col="red", pch = "x", lwd = 2)

# recreate using ggplot - these are just boxplots ;)
# there are 8 columns of data because I have 8 groups; these are colnames(pred_group_ML); and there are 2000 datapoints for each based on the number of iterations, burn-in and thinning

pred_SEA_B_df <- as.data.frame(pred_SEA_B) %>% 
  set_colnames(colnames(pred_group_ML)) %>% 
  pivot_longer(everything(), names_to = c('species', 'region'), names_sep = '[.]') %>% 
  mutate(region = factor(region, levels = c('BC', 'ON_west', 'ON_central', 'QC')), 
         species = factor(species))

str(pred_SEA_B_df)

# add credible interval
cred_intervals <- pred_SEA_B_df %>% 
  group_by(species, region) %>% 
  summarise(mean = mean(value),
            ci_low = ci(value, method = 'HDI')[,2], 
            ci_high = ci(value, method = 'HDI')[,3]) 

(pred_SEA_B_ggplot <- ggplot(pred_SEA_B_df, aes(x = species, y = value, color = species, fill = species))
  + geom_violin_pattern(aes(pattern_angle = species, pattern_color = species, pattern_fill = species, pattern_density = species))
  + stat_summary(fun = "mean",
                 geom = "point",
                 color = "black")
  + geom_hline(data = cred_intervals, aes(yintercept = ci_low, col = species), lty = 2, lwd = 0.75)
  + geom_hline(data = cred_intervals, aes(yintercept = ci_high, col = species), lty = 2, lwd = 0.75)
  + facet_grid(.~region)
  + scale_y_continuous(name = expression("standard ellipse area " ('\u2030' ^2) ), limits = c(0,12), breaks = c(2*0:7))
  + scale_colour_viridis_d(direction = -1, labels = c('Canada lynx', 'bobcat'))
  + scale_fill_viridis_d(direction = -1, alpha = 0.1, labels = c('Canada lynx', 'bobcat'))
  + scale_pattern_colour_viridis_d(direction = -1, labels = c('Canada lynx', 'bobcat'))
  + scale_pattern_fill_viridis_d(direction = -1, labels = c('Canada lynx', 'bobcat'))
  + scale_pattern_density_discrete(range = c(0.5, 0.01), labels = c('Canada lynx', 'bobcat'))
  + scale_pattern_angle_discrete(range = c(-45, 45), labels = c('Canada lynx', 'bobcat'))
  + theme_classic2()
  + theme(legend.position = 'top', legend.text = element_text(size = 8),  
          legend.box.margin = margin(0, 0, 0, 0), legend.margin = margin(0, 0, 0, 0), legend.title = element_blank())
  + theme(axis.title = element_text(size = 12), axis.text = element_text(size = 10), axis.text.x = element_blank(), axis.ticks.x = element_blank())
)

ggsave('Appendix_Fig2_SEAb_plot_wOutliers.tiff', pred_SEA_B_ggplot, dpi = 500, width = 7, height = 4, units = 'in')
# differences with Christa's results must be coming from random aspects of the Bayesian model fitting if outliers were manually removed from the data
# however, repeating the analysis with different seeds yields similar answers, which makes me wonder why these small differences are occurring. 

SEA_B_foldchange_outliers <- pred_SEA_B_df %>% 
  group_by(species, region) %>% 
  summarize(mean = mean(value), sd = sd(value)) %>% 
  pivot_wider(names_from = species, values_from = c(mean, sd), names_sep = '.') %>% 
  mutate(mean_foldchange = mean.Lynx_rufus/mean.Lynx_canadensis)

summary(SEA_B_foldchange_outliers)
