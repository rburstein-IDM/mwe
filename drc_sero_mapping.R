## ###########################################################
## Author:  Roy Burstein (rburstein@idmod.org)
## Date:    July 2019
## Purpose: Make a map of serology data DRC for immunity to Type 2
## ###########################################################


# libs
libs <- c('data.table','ggplot2','ggridges','lubridate','fasterize','gbm','dismo','rgeos',
          'raster','sf','mgcv','gridExtra','gstat','polycor','INLA')
for(l in libs) require(l,character.only=TRUE)

# locations - relative to user, but in AMUG Dropbox
#  Note: the MikeF path corresponds to data used for this post:
#      https://wiki.idmod.org/pages/viewpage.action?pageId=77041385&focusedCommentId=96077085
user      <- Sys.info()[['user']]
datpath   <- sprintf('C:/Users/%s/Dropbox (IDM)/AMUG/roy/scratch',user)
datpath2  <- sprintf('C:/Users/%s/Dropbox (IDM)/AMUG/MikeF/DRC 2013-14 DHS Serosurvey',user)

# load in data
d        <- fread(sprintf('%s/drc_dhs_polio/kids_master_dataset_v6.csv',datpath))
d2       <- fread(sprintf('%s/kids_and_adults.csv',datpath2))
shp      <- st_read(sprintf('%s/shapefile/cod_admbnda_adm1_rgc_20170711.shp',datpath2))
dpt1     <- raster(sprintf('%s/rasters/vaccine/dpt1_coverage/mean/IHME_AFRICA_DPT_2000_2016_DPT1_COVERAGE_PREV_MEAN_2014_Y2019M04D01.TIF',datpath))
dpt3     <- raster(sprintf('%s/rasters/vaccine/dpt3_coverage/mean/IHME_AFRICA_DPT_2000_2016_DPT3_COVERAGE_PREV_MEAN_2014_Y2019M04D01.TIF',datpath))
u5pop    <- raster(sprintf('%s/rasters/population/under5/worldpop_raked_a0004t_1y_2014_00_00.tif',datpath))
mat_edu  <- raster(sprintf('%s/rasters/education/years_f_15_49/mean/IHME_AFRICA_EDU_2000_2015_YEARS_FEMALE_15_49_MEAN_2014_Y2018M02D28.TIF',datpath))
imp_sani  <- brick(sprintf('%s/rasters/sanitation/s_imp_mean_ras_cod_0.tif',datpath))[[15]]


# helper functions stolen from seegMBG to use later on
insertRaster <- function (raster, new_vals, idx = NULL) {
  cellIdx <- function (x) 1:ncell(x)
  if (is.null(idx)) idx <- cellIdx(raster)
  stopifnot(length(idx) == nrow(new_vals))
  stopifnot(max(idx) <= ncell(raster))
  n <- ncol(new_vals)
  raster_new <- raster::brick(replicate(n,
                                        raster[[1]],
                                        simplify = FALSE))
  names(raster_new) <- colnames(new_vals)

  for(i in 1:n) {
    raster_new[[i]][idx] <- new_vals[, i]
  }
  return (raster_new)
}


## ###########################################################
## Clean up and prep data

# cut down the columns a bit
d <- d[, c('hh','psu','weight','province','sabin.1','sabin.2','sabin.3','lat','lon'), with = FALSE]

# na omit for now (note that in Arie's paper they did some raking for missingness)
d <- na.omit(d)
d <- d[lat!=0]

# aggregate to the cluster level (naive to weight for now)
dagg <- d[, .(type1 = sum(sabin.1>3),
              type2 = sum(sabin.2>3),
              type3 = sum(sabin.3>3),
              N     = .N), 
          by = .(psu,lat,lon,province)]

# make alternative variables
dpt_dropout <- dpt1 - dpt3
u5pop[u5pop==0] <- 0.01
log_u5pop   <- log(u5pop)

# make a raster version of the DRC shape
ext_raster <- trim(fasterize(shp, dpt1))

# clean up rasters and extract values at data locations
rlist <- c('dpt1','dpt3','dpt_dropout','u5pop','log_u5pop','mat_edu','imp_sani')
csfr  <- data.table() # to track center and scaling
for(r in rlist){
  
  message(sprintf(' .. crop it, mask it, extract it, center it, scale it: %s',r))

  # mask and extract
  assign(r, mask(crop(get(r), ext_raster), ext_raster))
  dagg[[r]] <- raster::extract(get(r), SpatialPoints(cbind(dagg$lon,dagg$lat)))
  
  # center and scale
  tmp  <- data.table(var=r,mean=mean(dagg[[r]], na.rm=T),sd=sd(dagg[[r]], na.rm=T))
  csfr <- rbind(csfr, tmp)
  dagg[[paste0(r,'_cs')]] <- ( dagg[[r]]-tmp$mean[1] ) / tmp$sd[1]
  
}

# make a prediction frame representing the full raster
predfr <- data.table(idx = 1:ncell(get(rlist[1])))
for(r in c(rlist)){ 
  predfr[[r]] <- as.vector(get(r))
  predfr[[paste0(r,'_cs')]] <- (predfr[[r]]-csfr[var==r]$mean) / csfr[var==r]$sd
}

# add lat and lon to prediction frame in case we use it
predfr[, lon := coordinates(ext_raster)[,1]]
predfr[, lat := coordinates(ext_raster)[,2]]




## ###########################################################
## Super simple sanity check model

# simple logistic regression
m1 <- glm( cbind(type2,N-type2) ~ dpt1+dpt_dropout+mat_edu+imp_sani+log_u5pop+lat*lon,
              data = dagg, family='binomial') 
summary(m1)

# logistic gam
m2 <- gam( cbind(type2,N-type2) ~ s(dpt1) + s(dpt_dropout) + s(mat_edu) + s(imp_sani) + s(log_u5pop) +
             s(lat) + s(lon),
             data = dagg, family='binomial', predict=TRUE )
summary(m2) # looks like non-linear terms are more predictive.


# how much explanatory power can we get out of a gbm?
incl <- c('type2','dpt1','dpt_dropout','mat_edu','imp_sani','log_u5pop','lat','lon')
m3 <- gbm.step(data    = dagg[,incl,with=F],
               gbm.y   = 1,
               gbm.x   = incl[-1],
               bag.fraction = 0.5,
               tree.complexity = 4,
               n.trees = 120,
               learning.rate = .001,
               offset  = log(dagg$N),
               family  = 'poisson')
               
# how did the various predictors do?
predtest <- data.table(type2 = dagg$type2/dagg$N,
                       pred1 = plogis(predict(m1,newdata=dagg)),
                       pred2 = plogis(predict(m2,newdata=dagg)),
                       pred3 = predict.gbm(m3,n.trees=m3$gbm.call$best.trees,type='response'))
cor(na.omit(predtest))

# map them
pred <- insertRaster(ext_raster,cbind(plogis(predict(m1,newdata=predfr)),
                                      plogis(predict.gam(m2,newdata=predfr)),
                                      plogis(predict(m3,n.trees=m3$gbm.call$best.trees,
                                                     type='response',newdata=predfr))))
plot(pred) # three simple predictions
                      


## #########################################################
## TODO: Incorporate SIAs as a spatial covariate.




## ###########################################################
## Basic geostats model - spatial covariates (enter linearly) with a spatial Matern GP spde
##   start with Sabin.2 in under-5s



# make spde object
spde <- inla.spde2.pcmatern(
          mesh=mesh, alpha=2, ### mesh and smoothness parameter
          prior.range=c(0.05, 0.01), ### P(practic.range<0.05)=0.01
          prior.sigma=c(1, 0.01)) ### P(sigma>1)=0.01
space <- inla.spde.make.index("space",
                             n.spde = spde$n.spde)

# make design matrix
narmdagg <- na.omit(dagg)
design_matrix <- data.frame(int = 1,
                            narmdagg)
outcome       <- narmdagg$type2
N             <- narmdagg$N               

# make a mesh
mesh <- inla.mesh.2d(loc      = cbind(narmdagg$lon, narmdagg$lat),
                     max.edge = c(.5,5),
                     offset   = 1,
                     cutoff   = .5)
plot(mesh); points(dagg$lon,dagg$lat,col='red')

# make projector matrix
A <- inla.spde.make.A(mesh, loc=cbind(narmdagg$lon, narmdagg$lat))


# make a data stack
stack.obs <- inla.stack(data    = list(outcome = outcome),
                        A       = list(A, 1),
                        effects = list(space,
                                       design_matrix),
                        tag     = 'est')

# 
formla <- outcome ~ 0 + int + dpt1 + dpt_dropout + 
  mat_edu + imp_sani + log_u5pop + f(space, model=spde) # + f(psu,model='iid')

# fit inla model
res_fit <- inla(formla,
                data = inla.stack.data(stack.obs),
                control.predictor = list(A       = inla.stack.A(stack.obs),
                                         link    = 1,
                                         compute = FALSE),
               # control.fixed = list(expand.factor.strategy = 'inla',
              #                       prec.intercept = int_prior_prec,
               #                      mean.intercept = int_prior_mn,
              #                       mean = fe_mean_prior,
              #                       prec = fe_sd_prior),
              #  control.inla = list(int.strategy = intstrat, h = 1e-3, tolerance = 1e-6),
              # weights = wgts,
                control.compute = list(dic = TRUE,
                                       cpo = TRUE,
                                       config = TRUE),
                family      = 'binomial',
                num.threads = 1,
                Ntrials     = N,
                verbose     = TRUE,
                keep        = TRUE)




# uncertainy

# population at risk






















