

library(data.table)
library(sf)
library(ggplot2)
library(maps) 

setwd('C:/Users/rburstein/Dropbox (IDM)/IHME/data')

# load and trim data
d <- fread('U5MR_LMICS_admin2_data_full.csv')

# Do we have Indonesia data comparable to the article which will be published ie: % of district/cities that have achieved SDG 3,2 targets  and performance/ progress 2000 - 2017 in neonatal, infant and U5 mortality?
d1 <- d[year == 2017 & ADM0_NAME == 'Indonesia']

mean(d1$u5mr_mean<.025)
sum(d1$u5mr_mean<.025)

mean(d1$u5mr_upper<.025)
sum(d1$u5mr_upper<.025)

mean(d1$u5mr_lower<.025)
sum(d1$u5mr_lower<.025)

mean(d1$nnmr_mean<.012)
sum(d1$nnmr_mean<.012)

mean(d1$nnmr_upper<.012)
sum(d1$nnmr_upper<.012)

mean(d1$nnmr_lower<.012)
sum(d1$nnmr_lower<.012)

