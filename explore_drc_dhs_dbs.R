## ###########################################################
## Author:  Roy Burstein (rburstein@idmod.org)
## Date:    July 2019
## Purpose: Explore DRC 2014 DHS with Linked DBS Polio data
##          Just understand and visualize the data. Make Maps
##          Try to model a surface
## ###########################################################


# libs
library(data.table)
library(ggplot2)
library(ggridges)
library(lubridate)
library(raster)


# locations - relative to user, but in AMUG Dropbox
#  Note: the MikeF path corresponds to data used for this post:
#      https://wiki.idmod.org/pages/viewpage.action?pageId=77041385&focusedCommentId=96077085
user      <- Sys.info()[['user']]
datpath   <- sprintf('C:/Users/%s/Dropbox (IDM)/AMUG/roy/scratch',user)
datpath2  <- sprintf('C:/Users/%s/Dropbox (IDM)/AMUG/MikeF/DRC 2013-14 DHS Serosurvey',user)

# load in data
d   <- fread(sprintf('%s/drc_dhs_polio/kids_master_dataset_v6.csv',datpath))
d2  <- fread(sprintf('%s/kids_and_adults.csv',datpath2))
shp <- shapefile(sprintf('%s/shapefile/cod_admbnda_adm1_rgc_20170711.shp',datpath2))
  
  
## ###########################################################
## Clean

# some are missing ages in the adult file but do have dob and doi, so we can calc their age 
d[, age_y := as.numeric(round((ymd(doi) - ymd(dob)) / 365,1))]
table(d$age, round(d$age_y), useNA = 'ifany')

# make one cleaned file with all ages sero data
setnames(d2,c('age','Sabin_1','Sabin_2','Sabin_3'),c('age_y','sabin.1','sabin.2','sabin.3'))
vnames <- c('cluster','hh','sabin.1','sabin.2','sabin.3','age_y','province','lat','lon','urban','urbanDetail')
dd <- rbind(d[,vnames,with=F],d2[adult==1,vnames,with=F])

# make age band ids
agebreaks <- c(0,1,2,3,4,5,10,15,20,25,30,35,40,45,50,55,60)
agelabels <- c("0-1","1-2","2-3","3-4","4-5","5-9","10-14","15-19","20-24","25-29","30-34",
               "35-39","40-44","45-49","50-54","55-59")
dd[ , agegroup := cut(age_y, breaks = agebreaks, right = FALSE, labels = agelabels)]
dd$age_min  <- as.numeric(do.call('rbind',strsplit(as.character(dd$agegroup),'-'))[,1])
dd$age_max  <- as.numeric(do.call('rbind',strsplit(as.character(dd$agegroup),'-'))[,2])


## ###########################################################
## Visualize and Map

# immunity across ages
ggplot(d, aes(x = sabin.2, y = factor(round(age_y)))) + geom_density_ridges() +
  theme_bw() + xlab('Log2 Type 2 Titer') + ylab('age') + 
  geom_vline(xintercept=3,color='red')

ggplot(d, aes(x = sabin.2, y = age_y)) + geom_hex(bins=20) +
  theme_bw() + xlab('Log2 Type 2 Titer') + ylab('age') + 
  geom_vline(xintercept=3,color='red')

# replicate the plot from Mike's blog post
ggplot(dd, aes(x=age_y,y=sabin.2,color=agegroup)) + theme_minimal() + 
  geom_point(position = position_jitter(w = 0.5, h = 0.15), alpha = .2, stroke = 0, shape = 16) + 
  geom_hline(yintercept=3,color='red') +
  geom_segment(data=dd[, .(mean_sabin.2 = mean(sabin.2,na.rm=TRUE)), by = .(age_min,age_max)],
               aes(x=age_min,xend=age_max,y=mean_sabin.2,yend=mean_sabin.2),color='black') +
  ylab('Log2 Type-2 Titer') + xlab('age') +
  theme(legend.position="none")

# replicate maps from Mike's blogpost
kids_agg <- dd[age_y <= 5, 
          .(sabin.1 = mean(sabin.1>3,na.rm=TRUE),
              sabin.2 = mean(sabin.2>3,na.rm=TRUE),
              sabin.3 = mean(sabin.3>3,na.rm=TRUE),
              N       = sum(!is.na(sabin.1))),
            by = .(province,lat,lon,cluster,urban,urbanDetail)]
adult_agg <- dd[age_y > 5, 
               .(sabin.1 = mean(sabin.1>3,na.rm=TRUE),
                 sabin.2 = mean(sabin.2>3,na.rm=TRUE),
                 sabin.3 = mean(sabin.3>3,na.rm=TRUE),
                 N       = sum(!is.na(sabin.1))),
               by = .(province,lat,lon,cluster,urban,urbanDetail)]

 # is there a correlation within clusters between kids and adults? within HHs?


## ###########################################################
## Check out correlation across Types




## ###########################################################
## Check out correlation between potential geospatial covariates and immunity




