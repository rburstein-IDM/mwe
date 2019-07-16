## ###########################################################
## Author:  Roy Burstein (rburstein@idmod.org)
## Date:    July 2019
## Purpose: Explore DRC 2014 DHS with Linked DBS Polio data
##          Just understand and visualize the data. Make Maps
##          Try to model a surface
## ###########################################################


# libs
libs <- c('data.table','ggplot2','ggridges','lubridate','raster','sf','mgcv','gridExtra')
for(l in libs) library(l,character.only=TRUE)

# locations - relative to user, but in AMUG Dropbox
#  Note: the MikeF path corresponds to data used for this post:
#      https://wiki.idmod.org/pages/viewpage.action?pageId=77041385&focusedCommentId=96077085
user      <- Sys.info()[['user']]
datpath   <- sprintf('C:/Users/%s/Dropbox (IDM)/AMUG/roy/scratch',user)
datpath2  <- sprintf('C:/Users/%s/Dropbox (IDM)/AMUG/MikeF/DRC 2013-14 DHS Serosurvey',user)

# load in data
d    <- fread(sprintf('%s/drc_dhs_polio/kids_master_dataset_v6.csv',datpath))
d2   <- fread(sprintf('%s/kids_and_adults.csv',datpath2))
shp  <- st_read(sprintf('%s/shapefile/cod_admbnda_adm1_rgc_20170711.shp',datpath2))
vacc1 <- raster(sprintf('%s/rasters/vaccine/dpt1_coverage/mean/IHME_AFRICA_DPT_2000_2016_DPT1_COVERAGE_PREV_MEAN_2014_Y2019M04D01.TIF',datpath))
vacc3 <- raster(sprintf('%s/rasters/vaccine/dpt3_coverage/mean/IHME_AFRICA_DPT_2000_2016_DPT3_COVERAGE_PREV_MEAN_2014_Y2019M04D01.TIF',datpath))
pop <- raster(sprintf('%s/rasters/population/under5/worldpop_raked_a0004t_1y_2014_00_00.tif',datpath))
edu <- raster(sprintf('%s/rasters/education/years_f_15_49/mean/IHME_AFRICA_EDU_2000_2015_YEARS_FEMALE_15_49_MEAN_2014_Y2018M02D28.TIF',datpath))


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

# crop rasters
vacc1 <- mask(crop(vacc1,shp), shp)
vacc3 <- mask(crop(vacc3,shp), shp)
pop   <- mask(crop(pop  ,shp), shp)
edu   <- mask(crop(edu  ,shp), shp)


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
dd[, child := paste0('Child = ',age_y < 6)]
agg <- dd[lat!=0,.(sabin.1 = mean(sabin.1>3,na.rm=TRUE),
                   sabin.2 = mean(sabin.2>3,na.rm=TRUE),
                   sabin.3 = mean(sabin.3>3,na.rm=TRUE),
                   N       = sum(!is.na(sabin.1))),
          by = .(province,lat,lon,cluster,urban,urbanDetail,child)]
agg <- melt(agg, id.vars = c('province','lat','lon','cluster','urban','urbanDetail','N','child'))
setnames(agg, 'value', 'pct_seropositive')

ggplot(agg) + geom_sf(data = shp, fill = 'white', colour = 'grey') + theme_minimal() + coord_sf(datum = NA) + 
  geom_point(aes(lon,lat,size=N,color=pct_seropositive),alpha=0.75, stroke = 0, shape = 16) + 
  facet_grid(child ~ variable) +
  scale_color_gradientn(values = c(0,0.8,1.0), colours = c("#a6172d","#EFDC05","#4f953b") ) +
  ylab('') + xlab('')




 # is there a correlation within clusters between kids and adults? within HHs?


## ###########################################################
## Check out correlation across Types



## ###########################################################
## Check out correlation across space




## ###########################################################
## Check out correlation between potential geospatial covariates and immunity (univariate plots for now)

# DPT1
agg$spatial_dpt1  <- raster::extract(vacc, SpatialPoints(cbind(agg$lon,agg$lat)))
d$spatial_dpt1   <- raster::extract(vacc1, SpatialPoints(cbind(d$lon,d$lat)))
d$spatial_dpt3   <- raster::extract(vacc3, SpatialPoints(cbind(d$lon,d$lat)))
d$spatial_edu    <- raster::extract(edu, SpatialPoints(cbind(d$lon,d$lat)))
d$spatial_logpop <- log(raster::extract(pop, SpatialPoints(cbind(d$lon,d$lat))))
d[, spatial_dpt_dropout := spatial_dpt1-spatial_dpt3] 

# predict prevalence of seropositivity
univ_mod_plot <- function(var, types = 1:3){
  
  varstr  <- var
  out <- data.table()
  
  for(t in types){
    tmp <- data.table(var = d[[varstr]], svp = d[[paste0('sabin.',t)]]>3)
  
    m    <- gam(svp~ s(var), family = 'binomial', data = tmp)
    pred <- predict.gam(m,newdata=tmp,se.fit=TRUE)
    tmp[, m := plogis(pred[[1]])]
    tmp[, l := plogis(pred[[1]]-pred[[2]]*1.96)]
    tmp[, u := plogis(pred[[1]]+pred[[2]]*1.96)]
    tmp[, type := t]
    out <- rbind(out,tmp)
  }  

  colz <- c('#30A9DE', '#E53A40', '#EFDC05'); leg <- 'Type'
  g1 <- ggplot(out,aes(x=var)) + geom_ribbon(aes(ymin=l,ymax=u,fill=factor(type)), alpha = 0.5) +
    geom_line(aes(y=m,color=factor(type))) + theme_minimal() + ylab('Predicted Seropositive Prevalence') +
    ylim(min(out$l),1) + xlab('') + scale_color_manual(values=colz, name=leg) + scale_fill_manual(values=colz, name=leg) +
    theme(legend.position = c(0.85, 0.1), legend.direction = "horizontal")
  g2 <- ggplot(out) + geom_density(aes(var),colour=NA,fill='grey') + theme_minimal() + xlab(varstr) + ylab("Data Density")
  grid.arrange(g1,g2, heights=c(2,1))
  
}  


univ_mod_plot(var='spatial_dpt1')
univ_mod_plot(var='spatial_dpt3')
univ_mod_plot(var='spatial_dpt_dropout')
univ_mod_plot(var='spatial_edu')
univ_mod_plot(var='spatial_logpop')

univ_mod_plot(var='age_y')

univ_mod_plot(var='tOPV')

univ_mod_plot(var='age_mother')




# predict mean log2 titer and get prediction intervals
m2 <- gam(sabin.2 ~ s(dpt1), data = d)

