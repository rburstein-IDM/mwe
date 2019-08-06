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
n_sitetypes <- 3



# set up param table
params <- expand.grid(admin    = paste0('AD_',1:n_admins),
                      pathogen = paste0('PATH_',letters[1:n_pathogens])) %>%
              data.table()


# for now, we assume cases are independent (the current method of using total other pathogens depends on this assumption, I think)
# later we can change this
# make them zero inflated but if have cases some random number of cases
params[, totalcases := rbinom(.N, 1, 0.9) * rpois(.N, exp(rnorm(.N, 2, 1)))]


# Simulate overall incidence over a 20 week period for each admin area and pathogen (again independently)
d <- data.table()
for(r in 1:nrow(params)){
  
  Nc <- params$totalcases[r] 
  
  if(Nc != 0){
  
    # an epi curve ( just do a simple gaussian distribution, moving mean and sd a bit for each admin-pathogen combo)
    peak <- abs(round(rnorm(1, 10, 3)))
    stdv <- abs(round(rnorm(1, 4,  1)))
    week <- floor(abs(round(rnorm(Nc, mean = peak, sd = stdv))))
    
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
  tmp <- d[pathogen != p, .(catchment = .N), by = .(admin, sitetype)]
  tmp$pathogen <- p
  catch <- rbind(catch, tmp)
}
d <- merge(d, catch, by = c('admin','sitetype','pathogen'), all.x = T)
d[is.na(catchment), catchment := 0] # in cases with no other pathogens present

# aggregate date by week admin site pathogen for the mode
dagg <- d[, .(cases = .N), by = .(admin,sitetype,pathogen,week,catchment)]

# exapnd the aggregated data for all possible combos
expanded <-  expand.grid(admin    = paste0('AD_',1:n_admins),
                         pathogen = paste0('PATH_',letters[1:n_pathogens]),
                         sitetype = paste0('SITE_',1:n_sitetypes),
                         week     = 0:max(d$week)) %>% data.table()
dagg <- merge(expanded, dagg, by = c('admin', 'pathogen', 'sitetype', 'week'), all.x = TRUE)
dagg[is.na(catchment), catchment := 0]
dagg[is.na(cases),     cases := 0]


########################################################################################
## Run model

formula <- 

  outcome ~ catchment + site_type + f(residence_neighborhood_district_nameRow, 
                                      model = "iid", hyper = modelDefinition$local, constr = TRUE, 
                                      replicate = replicateIdx, group = time_row_residence_neighborhood_district_name, 
                                      control.group = list(model = "rw2")) + f(time_row_rw2, model = "rw2", 
                                                                               hyper = modelDefinition$hyper$time, replicate = replicateIdx) + 
  f(time_row_IID, model = "iid", hyper = modelDefinition$hyper$local, 
    replicate = replicateIdx, constr = TRUE)
<environment: 0x1bb39078>
  
c("INLA::inla(formula = modelDefinition$formula, family = modelDefinition$family, ",  "    data = modelDefinition$inputData, Ntrials = modelDefinition$inputData$n, ",  "    verbose = TRUE, lincomb = modelDefinition$lincomb, control.compute = list(config = TRUE, ",  "        dic = TRUE), control.predictor = list(compute = TRUE, ",  "        link = 1), control.inla = list(int.strategy = \"auto\", ",  "        strategy = \"gaussian\"))")


