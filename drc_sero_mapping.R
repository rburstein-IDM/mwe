## ###########################################################
## Author:  Roy Burstein (rburstein@idmod.org)
## Date:    July 2019
## Purpose: Make a map of serology data DRC for immunity to Type 2
## ###########################################################


# libs
libs <- c('data.table','ggplot2','ggridges','lubridate','fasterize','gbm','dismo','rgeos',
          'raster','sf','mgcv','gridExtra','gstat','polycor','INLA','matrixStats')
for(l in libs) require(l,character.only=TRUE)

# locations - relative to user, but in AMUG Dropbox
#  Note: the MikeF path corresponds to data used for this post:
#      https://wiki.idmod.org/pages/viewpage.action?pageId=77041385&focusedCommentId=96077085
user      <- Sys.info()[['user']]
datpath   <- sprintf('C:/Users/%s/Dropbox (IDM)/AMUG/roy/scratch',user)
datpath2  <- sprintf('C:/Users/%s/Dropbox (IDM)/AMUG/MikeF/DRC 2013-14 DHS Serosurvey',user)
temptrash <- sprintf('%s/temptrash',datpath)


# load in data
d1       <- fread(sprintf('%s/drc_dhs_polio/kids_master_dataset_v6.csv',datpath))
d2       <- fread(sprintf('%s/kids_and_adults.csv',datpath2))
load(sprintf('%s/sia_draft_20190711.rdata',datpath))
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
d <- d2[, c('hh','cluster','adult','province','Sabin_1','Sabin_2','Sabin_3','lat','lon'), with = FALSE]

# na omit for now (note that in Arie's paper they did some raking for missingness)
d <- na.omit(d)
d <- d[lat!=0]

# aggregate to the cluster level (naive to weight for now)
dagg <- d[, .(type1 = sum(Sabin_1>3),
              type2 = sum(Sabin_2>3),
              type3 = sum(Sabin_3>3),
              N     = .N), 
          by = .(cluster,lat,lon,province,adult)]

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


## #########################################################
## TODO: Incorporate SIAs as a spatial covariate.

# bring in sia data
sia <- data.table(df_sia_final)
df_sia_final <- NULL

# keep to the 5 years preceding this survey and DRC only
sia <- sia[ADM0_NAME == 'DEMOCRATIC REPUBLIC OF THE CONGO']
sia[, year := year(ymd(start_date))]
sia <- sia[year >= 2009 & year < 2015]

# keep only 0 to 5 and faction >= .5
sia <- sia[agegroup == '0 to 5 years']
sia <- sia[fraction >= 0.5]

# collapse by ADM_2
sia <- sia[, .(type1_sias = sum(vaccinetpye %in% c('bOPV','mOPV1','tOPV')),
               type2_sias = sum(vaccinetpye %in% c('tOPV')),
               type3_sias = sum(vaccinetpye %in% c('tOPV'))) ]


## ###########################################################
## Super simple sanity check model

# simple logistic regression
m1ch <- glm( cbind(type2,N-type2) ~ dpt1+dpt_dropout+mat_edu+imp_sani+log_u5pop, # + lat*lon,
              data = dagg[adult==0], family='binomial')
m1ad <- glm( cbind(type2,N-type2) ~ dpt1+dpt_dropout+mat_edu+imp_sani+log_u5pop, # + lat*lon,
              data = dagg[adult==1], family='binomial')

summary(m1ch)
summary(m1ad)

# logistic gam
m2ch <- gam( cbind(type2,N-type2) ~ s(dpt1) + s(dpt_dropout) + s(mat_edu) + s(imp_sani) + s(log_u5pop), # + s(lat) + s(lon), 
             data = dagg[adult==0], family='binomial', predict=TRUE )
m2ad <- gam( cbind(type2,N-type2) ~ s(dpt1) + s(dpt_dropout) + s(mat_edu) + s(imp_sani) + s(log_u5pop), # + s(lat) + s(lon), 
             data = dagg[adult==1], family='binomial', predict=TRUE )
summary(m2ch) # looks like non-linear terms are more predictive.
summary(m2ad)

# how much explanatory power can we get out of a gbm?
incl <- c('type2','dpt1','dpt_dropout','mat_edu','imp_sani','log_u5pop') #,'lat','lon')
run_brt <- function(ddd)
  gbm.step(data            = ddd[,incl,with=F],
           gbm.y           = 1,
           gbm.x           = incl[-1],
           bag.fraction    = 0.5,
           tree.complexity = 4,
           n.trees         = 120,
           learning.rate   = .001,
           offset          = log(ddd$N),
           family          = 'poisson')
m3ch <- run_brt(ddd = dagg[adult==0])
m3ad <- run_brt(ddd = dagg[adult==1]) 

summary(m3ch)
summary(m3ad); plot.gbm(m3ad,4)

# how did the various predictors do?
predtest_ch <- data.table(type2 = dagg[adult==0]$type2/dagg[adult==0]$N,
                          pred1 = plogis(predict(m1ch,newdata=dagg[adult==0])),
                          pred2 = plogis(predict(m2ch,newdata=dagg[adult==0])),
                          pred3 = predict.gbm(m3ch,n.trees=m3ch$gbm.call$best.trees,type='response'))
predtest_ad <- data.table(type2 = dagg[adult==1]$type2/dagg[adult==1]$N,
                          pred1 = plogis(predict(m1ad,newdata=dagg[adult==1])),
                          pred2 = plogis(predict(m2ad,newdata=dagg[adult==1])),
                          pred3 = predict.gbm(m3ad,n.trees=m3ad$gbm.call$best.trees,type='response'))
cor(na.omit(predtest_ch))
cor(na.omit(predtest_ad))
# the covariates are a bit more predictive of child rather than adult immunity.
# the more complex brt model does about the same for both



# map them
pred_ch <- insertRaster(ext_raster,cbind(plogis(predict(m1ch,newdata=predfr)),
                                         plogis(predict.gam(m2ch,newdata=predfr)),
                                         predict(m3ch,n.trees=m3ch$gbm.call$best.trees,
                                                 type='response',newdata=predfr[,incl[-1],with=FALSE])))
pred_ad <- insertRaster(ext_raster,cbind(plogis(predict(m1ad,newdata=predfr)),
                                         plogis(predict.gam(m2ad,newdata=predfr)),
                                         predict(m3ad,n.trees=m3ad$gbm.call$best.trees,
                                                  type='response',newdata=predfr[,incl[-1],with=FALSE])))
plot(pred_ch[[3]]) # three simple predictions
plot(pred_ad[[3]])  







## ###########################################################
## Basic geostats model - spatial covariates (enter linearly) with a spatial Matern GP spde
##   start with Sabin.2 in under-5s


# make a wrapper function to do the model  heavy lifting

fitinlamodel <- function(data = dagg, fevars, outcometype){
  
  # make design matrix
  message(' .. data munge')
  narmdagg      <- na.omit(data)
  outcome       <- narmdagg[[outcometype]]
  N             <- narmdagg$N      
  locs          <- cbind(narmdagg$lon, narmdagg$lat)
  narmdagg      <- narmdagg[, fevars, with=FALSE]
  design_matrix <- data.frame(int = 1, narmdagg)
  
  # make a mesh
  mesh <- inla.mesh.2d(loc      = locs,
                       max.edge = c(.5,5),
                       offset   = 1,
                       cutoff   = 0.5)
  plot(mesh);  points(locs, col = 'red')
  
  
  # make spde object
  spde  <- inla.spde2.matern(mesh = mesh,  alpha = 2)
  space <- inla.spde.make.index("space", n.spde = spde$n.spde)
  
  
  # make projector matrix
  A <- inla.spde.make.A(mesh, loc = locs)
  
  
  # make a data stack
  stack.obs <- inla.stack(data    = list(o = outcome),
                          A       = list(A, 1),
                          effects = list(space, design_matrix),
                          tag     = 'est')
  
  
  # inla formula
  formla <- as.formula( sprintf ('o ~ -1 + int + %s + f(space, model = spde)', paste(fevars,collapse='+')) )
  
  
  # fit inla model
  message(' .. model fit')
  res_fit <- inla(formla,
                  data = inla.stack.data(stack.obs),
                  control.predictor = list(A       = inla.stack.A(stack.obs),
                                           link    = 1,
                                           compute = FALSE),
                  control.compute   = list(dic     = TRUE,
                                           cpo     = TRUE,
                                           config  = TRUE),
                  control.fixed     = list(expand.factor.strategy = 'inla'),
                  family            = 'binomial', #'gaussian',
                  num.threads       = 1,
                  Ntrials           = N,
                  verbose           = FALSE,
                  keep              = FALSE,
                  working.directory = temptrash)
  
  summary(res_fit)
  
  #########
  ##  make a prediction surface
  message(' .. prediction')
  draws <- inla.posterior.sample(1000, res_fit)
  
  ## get samples as matrices
  par_names <- rownames(draws[[1]]$latent)
  l_idx <- match(res_fit$names.fixed, par_names) # main effects
  s_idx <- grep('^space.*', par_names) # spatial effects
  pred_s <- sapply(draws, function (x) x$latent[s_idx])
  pred_l <- sapply(draws, function (x) x$latent[l_idx])
  rownames(pred_l) <- res_fit$names.fixed
  
  ## if we fit with a nugget, we also need to take draws of the nugget precision
  #if(length(grep('^IID.ID.*', par_names)) > 0){
  #  pred_n <- sapply(draws, function(x) {
  #    nug.idx <- which(grepl('IID.ID', names(draws[[1]]$hyper)))
  #    x$hyperpar[[nug.idx]]}) ## this gets the precision for the nugget
  #}else{
  #  pred_n <- NULL
  #}
  
  coords <- cbind(predfr$lon,predfr$lat)
  
  # Projector matrix
  A.pred <- inla.spde.make.A(mesh = mesh, loc = coords)
  
  predfr$int <- 1
  vals <- predfr[, c('int',fevars), with=FALSE]
  
  # get 1000 obs at each predframe location
  cell_l <- unname(as.matrix(as(data.matrix(vals), "dgeMatrix") %*% pred_l))
  cell_s <- as.matrix(crossprod(t(A.pred), pred_s))
  
  ## add on nugget effect if applicable
  #if(!is.null(pred_n)){
  #  cell_n <- sapply(pred_n, function(x){
  #    rnorm(n = nrow(cell_l), sd = 1 / sqrt(x), mean = 0)
  #  })
  #} else
  
  # make a cell pred object predframe by 1000 draws
  message(' .. summarize and save predictions')
  cell_pred <- plogis(cell_s+cell_l)
    
  # make summaries for plotting
  pred_q     <- rowQuantiles(cell_pred, probs = c(0.025, 0.975))
  inlapredfr <- data.frame(mean=rowMeans(cell_pred),lower=pred_q[,1],upper=pred_q[,2])
  predsumr   <- insertRaster(ext_raster,inlapredfr) 
  
  
  return(list(model_object = res_fit,
              cell_pred    = cell_pred,
              pred_summary = inlapredfr,
              raster       = predsumr))
  
}






# fit a separate child and adult model:
vars <- c('dpt1', 'dpt_dropout', 'mat_edu', 'imp_sani', 'log_u5pop')
type <- 'type2'

ch_inla_mod <- fitinlamodel(data = dagg[adult==0], fevars = vars, outcometype = type)
ad_inla_mod <- fitinlamodel(data = dagg[adult==1], fevars = vars, outcometype = type)






## plot it
hist(ch_inla_mod$inlapredfr[,1]) 
plot(ch_inla_mod$raster) #, zlim = c(0.5, 1.0))
        
hist(ad_inla_mod$inlapredfr[,1]) 
plot(ad_inla_mod$raster) #, zlim = c(0.2, 1.0)) #, main = 'Predicted adult seropositive prevalence')


# children minus adult
plot(ch_inla_mod$raster[[1]]-ad_inla_mod$raster[[1]], main = 'Additional immunity in children over adults')
         
diff_cp    <- ch_inla_mod$cell_pred - ad_inla_mod$cell_pred
diffq      <- rowQuantiles(diff_cp, probs = c(0.025, 0.975))
diffpredfr <- data.frame(mean=rowMeans(diff_cp),lower=diffq[,1],upper=diffq[,2])
diff       <- insertRaster(ext_raster,diffpredfr) 

plot(diff[[2:3]]) #, main = 'Predicted adult seropositive prevalence')


## TODO aggregate these results at the district level

## Do rank districts











######################################
# how do estimates compare to cluster level data?
res <- data.table(y = dagg[[outcometype]], N = dagg$N,
              raster::extract(pred, cbind(dagg$lon,dagg$lat)))
res[, prev := y/N +rnorm(nrow(res),0,0.01)] # jitter for plotting

ggplot(res, aes(x=prev,y=mean,ymin=lower,ymax=upper)) + theme_minimal() + 
  geom_abline(intercept=0,slope=1,color='red') +
  xlab('Empirical Estimate') + ylab('Model Estimate') +
  geom_errorbar(alpha=0.25) + geom_point(alpha=0.75,size=0.5)




######################################
# population at risk

pred_popatrisk <- inlapredfr*cbind(predfr$u5pop)

# about 1 million children at risk, type 2
sum(pred_popatrisk[,1], na.rm = TRUE)
sum(pred_popatrisk[,2], na.rm = TRUE)
sum(pred_popatrisk[,3], na.rm = TRUE)

plot(insertRaster(ext_raster,cbind(log(pred_popatrisk)[,1])), main='log population at risk')





## TODOS: 
# add a campaign covariate
# quality, residual analysis
# adults map -> joint fit with kids, and with other types















