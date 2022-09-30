
# Clean up of Christa's code for lynx-bobcat paper ('Code_Lynx-Bobcat_Diet_Overlap.txt')
# Author: Jenilee Gobin
# Date: 8 June, 2022

# Objective # 1 - investigate annual differences in isotope ratios to justify pooling of data across years
# Objective #2 - identify outliers
# Objective #3 - run analysis with and without outliers

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
                       "effectsize", 
                       "car", 
                       "broom", 
                       "MASS", 
                       "rstatix", 
                       "mvnormalTest", 
                       "mvtnorm", 
                       "dplyr", 
                       "purr", 
                       "psych",
                       "magrittr", 
                       "RVAideMemoire", 
                       "heplots",
                       "maps",
                       "mapdata", 
                       "maptools", 
                       "rnaturalearth", 
                       "sf", 
                       "dplyr", 
                       "scales", 
                       "RColorBrewer")

# apply function to packages
check.packages(analysis.packages)

# 2 - LOAD DATA ----- 

pred_dat <- read.csv("Table_DO_Lynx_Bobcat_2017_Oct1.csv", header=TRUE, sep=",") %>% 
  rename_with(tolower) %>% 
  rename(region = studyregion) %>% 
  mutate(region = factor(region), year = factor(year), species = factor(species))

head(pred_dat)
str(pred_dat)
summary(pred_dat)

# get summary stats
sum_stats1 <- pred_dat %>%
  group_by(species, region, year) %>%
  get_summary_stats(d13c, d15n, type = "mean_sd") %>% 
  view()

(sum_stats1_plot <- sum_stats1 %>% 
  ggplot()
  + geom_col(aes(x = as.factor(year), y = n, fill = region), position = 'dodge')
  + facet_grid(variable~species)
  + theme_bw()
  + geom_hline(yintercept = 3, lty = 2)
) # there is one instance of <3 data points per year for lynx, and 2 for bobcat


# pooled across years (see  1 for differences across years)
sum_stats2 <- pred_dat %>%
  group_by(species, region) %>%
  get_summary_stats(d13c, d15n, type = "mean_sd")

# get sample numbers byp species and region
sum_stats2 %>% filter(variable == 'd13c')

# 3 - VISUALIZE DATA -----

(plot1 <- ggplot(pred_dat, aes(x = d13c, y = d15n))
 + geom_point(aes(color = region, pch = region))
 + facet_grid(as.factor(year)~species)
 + theme_bw()
)

plot1b_dat <- pred_dat %>% 
  mutate(region = as.character(region), 
         region = replace(region, region == 'NWON', 'ON_west'), 
         region = replace(region, region == 'SSM', 'ON_central'), 
         region = factor(region, levels = c('BC', 'ON_west', 'ON_central', 'QC')), year = factor(year), 
         species = as.character(species), 
         species = replace(species, species == 'Lynx_canadensis', 'Canada lynx'), 
         species = replace(species, species == 'Lynx_rufus', 'bobcat'), 
         species = factor(species, levels = c('Canada lynx', 'bobcat')))
  

(plot1b <- ggplot(plot1b_dat, aes(x = d13c, y = d15n))
  + geom_point(aes(fill = year, pch = year), cex = 2)
  + facet_grid(species~region, scales = 'free')  
  #+ stat_ellipse(aes(col = year, lty = year), type = 'norm', lwd = 0.8, level = 0.95) # curious what annual clusters look like; looks like some differences in variance across years despite overlap
  + scale_fill_viridis_d(direction = -1, option = 'F')
  + scale_color_viridis_d(direction = -1, option = 'F')
  + scale_shape_manual(values = c(24, 21, 22, 25))
  + theme_classic2()
  + theme(legend.position = 'top', legend.text = element_text(size = 8),  
          legend.box.margin = margin(0, 0, 0, 0), legend.margin = margin(0, 0, 0, 0), legend.title = element_blank())
  + theme(axis.title = element_text(size = 12), axis.text = element_text(size = 10))
)

ggsave('Appendix_Fig0_databyyear.tiff', plot1b, dpi = 500, width = 7, height = 4, units = 'in')


# https://rpkgs.datanovia.com/ggpubr/reference/ggboxplot.html
(plot2 <- pred_dat %>% 
    ggboxplot(x = 'region', y = c('d13c', 'd15n'), color = 'year', facet.by = c('species'))
)
 
(plot2_lynx <- pred_dat %>% 
    filter(species == 'Lynx_canadensis') %>% 
    ggboxplot(x = 'year', y = c('d13c', 'd15n'), palette = 'grey', facet.by = c('region'))
)

(plot2_bobcat <- pred_dat %>% 
    filter(species == 'Lynx_rufus') %>% 
    ggboxplot(x = 'year', y = c('d13c', 'd15n'), palette = 'grey', facet.by = c('region'))
)

# 4A - OBJECTIVE #1 -----

# run manova

# cannot use car package because of missing data in certain years/regions:
pred_dat_manova_fitIII <- lm(cbind(d13c, d15n) ~ species * region * year, pred_dat,
             contrasts=list(species=contr.sum, region=contr.sum, year=contr.sum))
ManRes <- Manova(pred_dat_manova_fitIII, type="II")
summary(ManRes, multivariate=TRUE)

# note that manova() uses type I sum of squares and therefore is sensitive to the order in which predictor variables are specified
pred_dat_manova <- manova(cbind(d13c, d15n) ~ species * region * as.factor(year), pred_dat)
summary(pred_dat_manova) # default is Pillai, which is robust to violations of assumptions and should be just fine
summary(pred_dat_manova, test="Roy")
summary(pred_dat_manova, test="Wilks")
summary(pred_dat_manova, test="Hotelling-Lawley")
# outcomes are consistent across various test statistics

summary.aov(pred_dat_manova)
# both isotope ratios vary significantly  across species and regions
# annual differences are highly significant for C13 with significant species:year interaction term  
# annual differences are significant for d15n with significant interaction region:year interaction term and nearly significant species:year interaction

effectsize::eta_squared(pred_dat_manova) # large effect of species and region, and medium effect of year and species:year and region:year interactions

# check what happens if we change order of predictor variables: 
pred_dat_manova2 <- manova(cbind(d13c, d15n) ~ year * species * region, pred_dat)
summary(pred_dat_manova2) # default is Pillai, which is robust to violations of assumptions and should be just fine
summary(pred_dat_manova2, test="Roy")
summary(pred_dat_manova2, test="Wilks")
summary(pred_dat_manova2, test="Hotelling-Lawley")
# outcomes are consistent across various test statistics

summary.aov(pred_dat_manova2)

effectsize::eta_squared(pred_dat_manova2) # large effect of species and medium effect of year and species:year and region:year interactions
# if we specify year first in the model the amount of variance explained increases to 14%, with species and region explaining 35% and 29% of variability, respectively.
# proceed with initial manova

### SEE SMITH ET AL FOR USE OF MANOVA AND POSTHOC TESTS - AUG 5 2022 -----

# check univariate normality

# unable to test for years with <3 data points
# remove years with <3 data points (should also consider doing this before running the model at all)
univariate_shapiro <- pred_dat %>% 
  merge(sum_stats1 %>% filter(variable == 'd13c') %>% select(-c(variable))) %>% 
  filter(n>2) %>% 
  group_by(species, region, year) %>%
  shapiro_test(d13c, d15n) %>%
  arrange(species, variable, region, year) %>% 
  filter(p<0.05)
# in 2 of 40 instances p<0.05 but > 0.01; both cases are lynx (d13c BC 2011 and d15n NWON 2009)

# if assess at level of region
univariate_shapiro2 <- pred_dat %>%
  group_by(species, region) %>%
  shapiro_test(d13c, d15n) %>%
  arrange(species, variable, region)
# fail for bobcat d15n SSM (p<0.05)
# However, MANOVA is apparently fairly robust to deviations from normality so we will continue

# for sample sizer larger than 50, QQ-plots are recommended over Shapiro-Wilks, which becomes more sensitive to larger sample sizes

# QQ plot
(pred_dat_qqplot <- ggplot()
  + aes(sample = rstandard(pred_dat_manova))
  + stat_qq(distribution = qnorm)
  + stat_qq_line(line.p = c(0.25, 0.75))
  + theme_bw()
) # this works if pool data across years but not otherwise - likely due to numbers of values in each group

pred_dat_qqplot2 <- mqqnorm(pred_dat %>% select(d13c, d15n))
# departure from multivariate normality; 1 and 3 identified 

# Mardia's Skewness and Kurtois test
mvnormalTest::mardia(pred_dat[, c("d13c", "d15n")])$mv.test #fail

# manova fairly robust to deviations from normality, thus proceed. 

# check linearity
(plot3 <- ggplot(pred_dat, aes(x = d13c, y = d15n, color = region))
  + geom_point(aes(pch = region))
  + geom_smooth(method=lm, se=FALSE)
  + facet_grid(year~species)
  + theme_bw()
)

# multicollinearity assumption - correlation between d13c and d15c
multicor_all <- pred_dat %>% 
  cor_test(d13c, d15n) # overall correlation is ok (not too low or too high but significant p value obtained)

multicor_groups1 <- pred_dat %>% 
  group_by(species, region, year) %>% 
  cor_test(d13c, d15n)
# unable to test at level of year due to some years with too few data

# if remove years with <3 data points
multicor_groups1 <- pred_dat %>% 
  merge(sum_stats1 %>% filter(variable == 'd13c') %>% select(-c(variable))) %>% 
  filter(n>2) %>% 
  group_by(species, region, year) %>% 
  cor_test(d13c, d15n) 

summary(multicor_groups1)
filter(multicor_groups1, abs(cor) >0.9)
#one instance of a correlation coefficient >0.9 - bobcat SSM 2012 due to n = 3

multicor_groups2 <- pred_dat %>% 
  group_by(species, region) %>% 
  cor_test(d13c, d15n)
# no correlation coefficients >0.9

# homogeneity of variance

#Box's M-test - compares covariance matrices
# significance assessed at p<0.001
# this becomes problematic by group because do not have data for all years for all regions and species; also expects that n > p, where p is the number of variables which is not always the case
box_m <- boxM(cbind(d13c, d15n) ~ species * region * year, data = pred_dat)

# again, just try removing years with just one data point:
pred_dat_min3 <- pred_dat %>% 
  merge(sum_stats1 %>% filter(variable == 'd13c') %>% select(-c(variable))) %>% 
  filter(n>2) 

box_m <- boxM(cbind(d13c, d15n) ~ species * region * year, data = pred_dat_min3) # apparently this also includes zeros, which we can't address

# can assess at regional level ...

box_m_pred_dat <- boxM(cbind(d13c, d15n) ~ species * region , data = pred_dat) # fails

# if broken down by species
box_m_lynx_dat <- pred_dat %>% 
  filter(species == 'Lynx_canadensis') %>%
  select(-c(year, envelope, species))

box_m_lynx <- box_m(box_m_lynx_dat[,-1], box_m_lynx_dat[,1]) # highly significant for lynx suggesting covariance matrices are not equal - i.e., test failed
# sample sizes are also unbalanced across groups for lynx:
sum_stats2
# this test is also sensitive to large sample sizes (i.e., can report significant results when one does not exist) 

box_m_bobcat_dat <- pred_dat %>% 
  filter(species == 'Lynx_rufus') %>%
  select(-c(year, envelope, species))

box_m_bobcat <- box_m(box_m_bobcat_dat[,-1], box_m_bobcat_dat[,1]) # non-significant for bobcat

# levene's test
levene_test_pred_dat <- pred_dat %>% 
  gather(key = "variable", value = "value", d13c, d15n) %>%
  group_by(variable) %>%
  levene_test(value ~ species*region*year)
# fails for both d13c and d15c

# post-hoc tests

# linear discriminant analysis (LDA)
# not sure if possible to do this with multiple grouping variables? Is MDA an option here???????

# split data into groups
pred_dat_split <- pred_dat %>% 
  group_by(species, region) %>% 
  group_split()

# lynx BC
lynx_BC_lda <- lda(year~cbind(d13c, d15n), data = pred_dat_split[[1]]) %>% 
  predict() %>% 
  data.frame(pred_dat_split[[1]][, "year"]) 
  
summary(lynx_BC_lda$class == lynx_BC_lda$year) 
as.integer(summary(lynx_BC_lda$class == lynx_BC_lda$year)["TRUE"])/length(lynx_BC_lda$year)*100
    
(lynx_BC_lda_plot <- lynx_BC_lda %>% ggplot()
  + geom_point(aes(x = x.LD1, y = x.LD2, colour = year), size = 3)
  + theme_bw()
)
# don't see super obvious groupings by year and can predict year correctly at rate of 57%

# lynx NWON
lynx_NWON_lda <- lda(year~cbind(d13c, d15n), data = pred_dat_split[[2]]) %>% 
    predict() %>% 
    data.frame(pred_dat_split[[2]][, "year"]) 

summary(lynx_NWON_lda$class == lynx_NWON_lda$year)
as.integer(summary(lynx_NWON_lda$class == lynx_NWON_lda$year)["TRUE"])/length(lynx_NWON_lda$year)*100

(lynx_NWON_lda_plot <- lynx_NWON_lda %>% ggplot()
  + geom_density(aes(LD1, fill = year), alpha = .2)
  + theme_bw()
)

# fairly substantial overlap and predict year correctly at rate of 58%

# lynx QC
lynx_QC_lda <- lda(year~cbind(d13c, d15n), data = pred_dat_split[[3]]) %>% 
    predict() %>% 
    data.frame(pred_dat_split[[3]][, "year"]) 
  
summary(lynx_QC_lda$class == lynx_QC_lda$year)
as.integer(summary(lynx_QC_lda$class == lynx_QC_lda$year)["TRUE"])/length(lynx_QC_lda$year)*100
  
(lynx_QC_lda_plot <- lynx_QC_lda %>% ggplot()
  + geom_density(aes(LD1, fill = year), alpha = .2)
  + theme_bw()
) 
# fairly substantial overlap but predict year correctly at rate of 84%

# lynx SSM
lynx_SSM_lda <- lda(year~cbind(d13c, d15n), data = pred_dat_split[[4]]) %>% 
    predict() %>% 
    data.frame(pred_dat_split[[4]][, "year"]) 

summary(lynx_SSM_lda$class == lynx_SSM_lda$year)
as.integer(summary(lynx_SSM_lda$class == lynx_SSM_lda$year)["TRUE"])/length(lynx_SSM_lda$year)*100
  
(lynx_SSM_lda_plot <- lynx_SSM_lda  %>% ggplot()
  + geom_point(aes(x = x.LD1, y = x.LD2, colour = year), size = 3)
  + theme_bw()
) 
# no obvious groupings and predict year correctly at rate of 72%; only one datapoint for 2012

# bobcat BC
bobcat_BC_lda <- lda(year~cbind(d13c, d15n), data = pred_dat_split[[5]]) %>% 
    predict() %>% 
    data.frame(pred_dat_split[[5]][, "year"]) 

summary(bobcat_BC_lda$class == bobcat_BC_lda$year)
as.integer(summary(bobcat_BC_lda$class == bobcat_BC_lda$year)["TRUE"])/length(bobcat_BC_lda$year)*100
  
(bobcat_BC_lda_plot <- bobcat_BC_lda %>% ggplot()
  + geom_point(aes(x = x.LD1, y = x.LD2, colour = year), size = 3)
  + theme_bw()
)
# no obvious groupings and predict year correctly at rate of 56%, only one datapoint for 2012

# bobcat NWON
bobcat_NWON_lda <- lda(year~cbind(d13c, d15n), data = pred_dat_split[[6]]) %>% 
    predict() %>% 
    data.frame(pred_dat_split[[6]][, "year"]) 
  
summary(bobcat_NWON_lda$class == bobcat_NWON_lda$year)
as.integer(summary(bobcat_NWON_lda$class == bobcat_NWON_lda$year)["TRUE"])/length(bobcat_NWON_lda$year)*100  
  
(bobcat_NWON_lda <- bobcat_NWON_lda  %>% ggplot()
  + geom_point(aes(x = x.LD1, y = x.LD2, colour = year), size = 3)
  + theme_bw()
)
# no super obvious groupings and predict year correctly at rate or 58%, only one datapoint for 2011 

# bobcat QC
bobcat_QC_lda <- lda(year~cbind(d13c, d15n), data = pred_dat_split[[7]]) %>% 
    predict() %>% 
    data.frame(pred_dat_split[[7]][, "year"]) 

summary(bobcat_QC_lda$class == bobcat_QC_lda$year)
as.integer(summary(bobcat_QC_lda$class == bobcat_QC_lda$year)["TRUE"])/length(bobcat_QC_lda$year)*100 

(bobcat_QC_lda <- bobcat_QC_lda %>% ggplot()
  + geom_point(aes(x = x.LD1, y = x.LD2, colour = year), size = 3)
  + theme_bw()
)
# 2012 separated on LD1 an predict year correctly at rate of 78%

# bobcat SSM
bobcat_SSM_lda <- lda(year~cbind(d13c, d15n), data = pred_dat_split[[8]]) %>% 
    predict() %>% 
    data.frame(pred_dat_split[[8]][, "year"]) 
  
summary(bobcat_SSM_lda$class == bobcat_SSM_lda$year)
as.integer(summary(bobcat_SSM_lda$class == bobcat_SSM_lda$year)["TRUE"])/length(bobcat_SSM_lda$year)*100 
    
(bobcat_SSM_lda <- bobcat_SSM_lda  %>% ggplot()
  + geom_density(aes(LD1, fill = year), alpha = .2)
  + theme_bw()
)
# no obvious groupings and predict year correctly at rate of 63% (note there are only 8 data points in total), 2009 distribution really  flat though

# 4B - OBJECTIVE #2 -----

# multivariate outliers

md <- mahalanobis_distance(data = pred_dat[, c("d13c", "d15n")])$is.outlier

pred_dat_md <- pred_dat %>% 
  cbind(md) %>% 
  filter(md == 'TRUE')
# identified 4 bobcats in BC: 1 from 2009, 1 from 2010, 2 from 2011

# assess by group
pred_dat_md_grouped <- pred_dat %>%
  group_by(species, region, year) %>%
  select(-envelope) %>% 
  mahalanobis_distance() %>%
  filter(is.outlier == TRUE) %>%
  as.data.frame()
# outcome is the same as above

# double check if this is still the case if split groups because they are not the same as the outliers identified by Christa

pred_dat_split <- pred_dat %>% 
  group_by(species, region) %>% 
  select(-envelope, year) %>% 
  group_split()

lynx_BC_md <- pred_dat_split[[1]][,4:5] %>% 
  mahalanobis_distance() %>%
  filter(is.outlier == TRUE)

lynx_NWON_md <- pred_dat_split[[2]][,4:5] %>% 
  mahalanobis_distance() %>%
  filter(is.outlier == TRUE)

lynx_QC_md <- pred_dat_split[[3]][,4:5] %>% 
  mahalanobis_distance() %>%
  filter(is.outlier == TRUE)

lynx_SSM_md <- pred_dat_split[[4]][,4:5] %>% 
  mahalanobis_distance() %>%
  filter(is.outlier == TRUE)

bobcat_BC_md <- pred_dat_split[[5]][,4:5] %>% 
  mahalanobis_distance() %>%
  filter(is.outlier == TRUE)

bobcat_NWON_md <- pred_dat_split[[6]][,4:5] %>% 
  mahalanobis_distance() %>%
  filter(is.outlier == TRUE)

bobcat_QC_md <- pred_dat_split[[7]][,4:5] %>% 
  mahalanobis_distance() %>%
  filter(is.outlier == TRUE)

bobcat_SSM_md <- pred_dat_split[[8]][,4:5] %>% 
  mahalanobis_distance() %>%
  filter(is.outlier == TRUE)

# outliers have been removed from the dataset at a whole because we don't have the datapoints she cited in the text in the data; 
# note that mahalanobis_distance() does not appear to work on grouped data

# NOTE THAT CHRISTA HAS HARE STABLE ISOTOPE DATA BUT WAS NOT INCLUDED IN THE FILE TRANSFER - 
# COULD BE LOCATED IN "C:/Users/Christa/Google Drive/UofM/RFiles" OF LAPTOP OR "C:/Users/Schmorb/Google Drive/UofM/RFiles" OF DESKTOP
# SHOUDLD ALSO HAVE ISOTOPE DATA FOR RED SQUIRRELS


# 4C - OBJECTIVE #3 -----
# SEE 2_SIBER_ANALYSIS_JL.R

# 5 - MAP SAMPLES -----

# maps (https://r-graph-gallery.com/278-the-maps-library.html)

usa <- map("usa")
canada <- map("canada.cities")

# rnaturalearth and sf (https://community.rstudio.com/t/combining-usa-canada-specific-province-with-plotly-plot-geo-or-ggplot2/12375/2)

# use crs 4326 
canada <- ne_states(country = "canada", returnclass = "sf") %>%
  st_transform(crs = 4326)

usa <- ne_states(country = "united states of america", returnclass = "sf") %>%
  st_transform(crs = 4326)

# read in lynx and bobcat location data and combine with pred_dat

lynx_locations <- read.csv('TotalLynxLocations.csv')
bobcat_locations <- read.csv('TotalBobcatLocations.csv')

str(lynx_locations)
str(bobcat_locations)
str(pred_dat)

all_location_dat <- bobcat_locations %>% 
  mutate(Species = 'Lynx_rufus') %>% 
  select(c(ID, Species, Lat, Long)) %>% 
  rename_with(~ tolower(.)) %>% 
  rbind(lynx_locations %>% 
          select(c(ID, Species, X, Y)) %>% 
          rename(lat = X, long = Y) %>% 
          rename_with(~tolower(.)) %>% 
          mutate(species = 'Lynx_canadensis')) %>% 
  rename(envelope = id) %>% 
  mutate(envelope = tolower(envelope))

tolower(pred_dat$envelope) %in% all_location_dat$envelope

location_matches <- pred_dat %>% 
  filter(tolower(envelope) %in% all_location_dat$envelope)
str(location_matches)

# no match for one envelope name
(no_matches <- pred_dat %>% 
  filter(!(tolower(envelope) %in% all_location_dat$envelope))
  )

# go in the opposite direction to get summary of lat/long dat available owing to NAs being introduced upon merge
location_matches2 <- all_location_dat %>% 
  filter(tolower(envelope) %in% tolower(pred_dat$envelope)) %>% 
  filter(is.na(lat) | is.na(long))

summary(location_matches2) # of the location data that matches envelopes/ids in pred_dat, 27 have NA values for lat or long

# this is because these bobcat locations are planar x y...
location_matches2b <- bobcat_locations %>% 
  filter(tolower(ID) %in% tolower(pred_dat$envelope)) %>% 
  filter(is.na(Lat)|is.na(Long))

summary(location_matches2b)

# merge for now and then re-run when location data has been added for bobcat 11a774; 
# also map existing points to determine UTM zone for bobcat that need conversion to lat/long

pred_dat_locations <- pred_dat %>% 
  mutate(envelope = tolower(envelope)) %>% 
  merge(all_location_dat, by = c('envelope', 'species'), all.x = TRUE) 

str(pred_dat_locations)
summary(pred_dat_locations) # this is from the 27 NAs above and the one bobcat id with no match (11a774)

na_locations <- pred_dat_locations %>% filter(is.na(lat)|is.na(long)) %>% filter(envelope != '11a774')
na_locations %>% select(c(species, region)) %>% table()
# those with planar coordinates are bobcats from BC and NWON regions...

# go back to bobcat location data and convert
convert_dat <- bobcat_locations %>% 
  filter(tolower(ID) %in% na_locations$envelope) %>% 
  select(c(ID, Planar_X, Planar_Y)) %>% 
  rename(envelope = ID) %>% 
  rename_with(~tolower(.)) %>% 
  mutate(envelope = tolower(envelope)) %>% 
  merge(pred_dat_locations)

# projections are NAD 1983 Lambert North America
crs_string <- '+proj=lcc +lat_0=40 +lon_0=-96 +lat_1=20 +lat_2=60 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +type=crs'

# incorporate location data for sample 11a774 that is not found in Arthur's location data but is found Christa's 'Geolocations.xls' file
(dat_11a774 <- pred_dat %>% 
  filter(envelope == '11a774') %>% 
  mutate(planar_x = -1531995.80766, planar_y = 1256984.57956, lat = NA, long = NA) %>% 
  relocate(envelope, planar_x, planar_y, species, region, year, .before = d13c)
)
  
  
(convert_dat_sf_bc <- convert_dat %>% 
  filter(region == 'BC') %>% 
  rbind(dat_11a774) %>% 
  st_as_sf(coords = c("planar_x", "planar_y"), crs = crs_string) %>% 
  st_transform(crs = 4326) %>% 
  select(-c(lat, long))
)

(convert_dat_sf_on <- convert_dat %>% 
  filter(region == 'NWON') %>% 
  st_as_sf(coords = c("planar_x", "planar_y"), crs = crs_string) %>% 
  st_transform(crs = 4326) %>% 
  select(-c(lat, long))
)

str(convert_dat)

# make sf object
pred_dat_locations_sf <- pred_dat %>% 
  merge(all_location_dat, by = c('envelope', 'species'), all.x = TRUE) %>% 
  filter(!(is.na(lat)|is.na(long))) %>% 
  st_as_sf(coords = c("long", "lat"), crs = 4326) %>% 
  rbind(convert_dat_sf_bc, convert_dat_sf_on)


land_cols <- hcl.colors(n=4, palette = 'Earth')

# Plot the map
(map <- ggplot() 
  + geom_sf(data = canada, fill = land_cols[3])
  + geom_sf(data = usa, fill = land_cols[2])
  + geom_text(data = canada, aes(label = postal, x = longitude, y = latitude), size = 3, check_overlap = T, vjust = "inward")
  + geom_text(data = usa, aes(label = postal, x = longitude, y = latitude), size = 3, check_overlap = T, vjust = "inward")
  + geom_sf(data = pred_dat_locations_sf, aes(fill = species, pch = species), cex = 2)
  #+ geom_sf(data = convert_dat_sf_bc, col = 'red')
  #+ geom_sf(data = convert_dat_sf_on, col = 'blue')
  + coord_sf(xlim = c(-130,-66), ylim = c(43,57))
  + scale_fill_viridis_d(direction = -1, labels = c('Canada lynx', 'bobcat'))
  + scale_shape_manual(values = c(21,24), labels = c('Canada lynx', 'bobcat'))
  + theme_void()
  + theme(legend.position = c(0.999,0.995), legend.justification = c('right', 'top'), legend.title = element_blank(), legend.background = element_rect(fill = alpha('white', 0.9)), 
          legend.text = element_text(size = 8), legend.margin = margin(r=0.4, unit='cm'), legend.direction = 'horizontal')
  + theme(panel.background = element_rect(fill = 'lightblue', colour = 'black'), panel.border = element_rect(color = "black", size = 0.8))
  #+ theme(axis.title = element_blank(), panel.grid = element_blank())
  + theme(axis.text = element_text(size = 8, margin = margin(2,2,2,2)), axis.ticks = element_line(), axis.ticks.length = unit(2, "pt"))
)

ggsave('Fig0_map.tiff', map, dpi = 500, width = 7, height = 3, units = 'in')
ggsave('Fig0_map_large.tiff', map, dpi = 500, width = 14, height = 4, units = 'in')


# centroids in this context are the coordinates of the center of polygon defined by the outermost points for each region
# (presumably including both bobcat and lynx for each region)





















