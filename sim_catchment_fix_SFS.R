########################################################################################
# Author:   Roy Burstein
# Date:     August 2019
# Purpose:  Simulate a simple version of Mike's catchment model to answer the question:
#           is the remaining 'latent field' proportional to incidence?
########################################################################################

library(data.table)
library(magrittr)
library(ggplot2)
library(INLA)


########################################################################################
## Simulate data (admin, site, pathogen cases by week)

### set params

# number of admin areas
n_admins <- 5

# number of pathogens tracked
n_pathogens <- 3

# number of site types
n_sitetypes <- 5



# set up param table
params <- expand.grid(admin    = paste0('AD_',1:n_admins),
                      pathogen = paste0('PATH_',letters[1:n_pathogens])) %>%
              data.table()


# for now, we assume cases are independent (the current method of using total other pathogens depends on this assumption, I think)
# later we can change this
# make them zero inflated but if have cases some random number of cases
params[, totalcases := rbinom(.N, 1, 0.9) * rpois(.N, exp(rnorm(.N, 7, 1)))]


# Simulate overall incidence over a 20 week period for each admin area and pathogen (again independently)
d <- data.table()
for(r in 1:nrow(params)){
  
  Nc <- params$totalcases[r] 
  
  if(Nc != 0){
  
    # an epi curve ( just do a simple gaussian distribution, moving mean and sd a bit for each admin-pathogen combo)
    peak <- abs(round(rnorm(1, 10, 3)))
    stdv <- abs(round(rnorm(1, 4,  1)))
    week <- ceiling(abs(round(rnorm(Nc, mean = peak, sd = stdv))))
    
    # for each case, draw a sitetype it was detected at (for now assume equal probability detection at each site)
    # this is pretty unrealistic in the real data so we'll need to think aboiut the effect of this assumption (what is most cases are from one site one admin (ie UW))
    styp <- paste0('SITE_',ceiling(runif(Nc,0,n_sitetypes)))
    
    # append to d
    out <- data.table(admin = params$admin[r], pathogen = params$pathogen[r], week= week, sitetype = styp)
    d   <- rbind(d, out)
  } 
}

d <- d[order(admin,pathogen,week)]

# plot to see whats happening
ggplot(d, aes(week,fill=sitetype)) + geom_histogram() + facet_grid(admin~pathogen)


# get catchment as currently defined (all pathogens by Admin, site other than the pathogen in question)
# late can do a smooth model like mike did for this and then wont get any zeroes 
catch <- data.table()
for(p in unique(d$pathogen)){
  tmp <- d[pathogen != p, .(catchment = .N+0.001), by = .(admin, sitetype)]
  tmp$pathogen <- p
  catch <- rbind(catch, tmp)
}
d <- merge(d, catch, by = c('admin','sitetype','pathogen'), all.x = T)
d[is.na(catchment), catchment := 0.1] # in cases with no other pathogens present

# aggregate date by week admin site pathogen for the mode
dagg <- d[, .(cases = .N), by = .(admin,sitetype,pathogen,week,catchment)]

# exapnd the aggregated data for all possible combos
expanded <-  expand.grid(admin    = paste0('AD_',1:n_admins),
                         pathogen = paste0('PATH_',letters[1:n_pathogens]),
                         sitetype = paste0('SITE_',1:n_sitetypes),
                         week     = 1:max(d$week)) %>% data.table()
dagg <- merge(expanded, dagg, by = c('admin', 'pathogen', 'sitetype', 'week'), all.x = TRUE)
dagg[is.na(catchment) , catchment := 0.1]
dagg[is.na(cases),     cases := 0]

dagg[, catchment := log(catchment)-mean(log(catchment))]

########################################################################################
## Run model, try to do similar to mike's iid model as possible

## model one pathogen at a time.
path_to_model <- 'PATH_c'

# subset
inputData <- dagg[pathogen == path_to_model]

# priors
hyper <- list()
hyper$global <- list(prec = list( prior = "pc.prec", param = 1/10, alpha = 0.01))
hyper$local  <- list(prec = list( prior = "pc.prec", param = 1/200, alpha = 0.01))
hyper$age    <- list(prec = list( prior = "pc.prec", param = 1, alpha = 0.01))
hyper$time   <- list(prec = list( prior = "pc.prec", param = 1/50, alpha = 0.01))

family <- 'poisson'

outcome <- inputData$cases

# initialize formula  
formula <- as.formula('cases ~ 1 + catchment + sitetype')

# time as a latent fied
inputData$time_row_rw2 <- inputData$week
inputData$time_row_IID <- inputData$week

inputData$admin_row <- match(inputData$admin,unique(inputData$admin))
inputData$time_row_admin <- inputData$week
  
formula <- update(formula,  ~ . + f(admin_row, model='iid', hyper=hyper$local, # constr = TRUE, 
                                     group = time_row_admin, control.group=list(model="rw2")))
formula <- update(formula,  ~ . + f(time_row_rw2, model='rw2', hyper= hyper$time) ) # +
                 # f(time_row_IID, model='iid', hyper=hyper$local,  constr = TRUE) )

validLatentFieldColumns <- c('admin_row','time_row_admin','time_row_rw2') #,'time_row_IID')
  
# linear combination of pathogen and latent fields

# find unique rows after discarding factors that are being averaged over
lc.data <- data.frame(inputData[,names(inputData) %in% validLatentFieldColumns, with = FALSE])
lc.rowIdx <- !duplicated(lc.data)
lc.data <- lc.data[lc.rowIdx,]

# generate list of desired linear combinations # https://groups.google.com/forum/#!topic/r-inla-discussion-group/_e2C2L7Wc30
lcIdx=c()
spentColumn<-rep(FALSE,length(validLatentFieldColumns))
for(COLUMN in validLatentFieldColumns){
  if(!spentColumn[validLatentFieldColumns %in% COLUMN]) {
    lcIdx[[COLUMN]] <- inla.idx(lc.data[[COLUMN]])          
  }
  
  spentColumn[validLatentFieldColumns %in% COLUMN]<-TRUE
}

# generate list of desired linear combinations # https://groups.google.com/forum/#!topic/r-inla-discussion-group/_e2C2L7Wc30
lc.latentField <- vector("list", nrow(lc.data))

w<-vector("list", length(names(lcIdx))+1)
w[[length(names(lcIdx))+1]]<-1 #pathogen

for(k in 1:nrow(lc.data)){
  
  for(n in 1:length(names(lcIdx))){
    w[[n]]<-rep(0,nrow(lc.data))
    w[[n]][lcIdx[[n]][k]]<-1
  }
  names(w) <- c(names(lcIdx),'(Intercept)')
  
  lc <- inla.make.lincomb(w)
  names(lc)<- paste0('latent_field',k)
  lc.latentField[k]<-lc
  lc.data$latentField[k]<-names(lc)

}


# get original values for linear combination categories
lc.colIdx <- (names(inputData) %in% c('admin','week'))
lc.data <-inputData[which(lc.rowIdx),c('admin','week'),with=FALSE]

#fit
model <- INLA::inla(formula           = formula,
                    family            = family, 
                    data              = inputData, 
                    lincomb           = lc.latentField,
                   # Ntrials           = inputData$cases,
                    control.predictor = list(compute=TRUE,link=1),
                    control.compute   = list(config=TRUE,dic=TRUE),
                    verbose           = TRUE,
                    control.inla      = list(int.strategy="auto", strategy = "gaussian"))

