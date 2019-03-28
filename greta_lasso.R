## Roy Burstein
## March 2019
## Example of Bayesian Hierachical Lasso Models in R



## INSTRUCTIONS:
## Make a folder called test on your desktop. Put the data in it. We will link this volume to a docker container
## Run the following command in powershell:
#    docker run -e PASSWORD=GRETAFORSTEVE -v C:\Users\<<YOUR USERNAME>>\Desktop\test:/home/rstudio/test -d -p 8787:8787 earthlab/r-greta
## Once the above line has successfully run. Type the following in your Google Chrome address bar:
#    http://localhost:8787/
#    (When prompted, use rstudio as a username and GRETAFORSTEVE as a password)
#    This should open up an RStudio environment in chome!


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load  packages

library(greta)



####################################################################################
### FIRST THING: SIMULATE DATA AND RUN MODEL TO DEMONSTRATE  
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Simulate some data
set.seed(54321)

# number of groups
ng <- 20

# set parameters
a  <- 1              # global intercept
re <- rnorm(ng,0,5) # random intercepts

# set fixed effects. some are zeroes, well see if lasso picks them out
bs <- c(20,0,0,30,10)

# group ID (uneven number)
group_id <- rep(1:ng,rpois(ng,250)) # uneven number of groups like in steves data
N        <- length(group_id)

# indep variables. Simulate a matrix
X       <- mapply(rep(N,length(bs)), FUN = function(x) {rnorm(x,0,1)} )

# response variable
Y <- a + X%*%bs + rnorm(N,0,3) +  re[group_id] #

# plot data
hist(Y)
plot(X[,4],Y,col=group_id)



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Fit with greta

# priors
int          <- normal(0, 10)
sd           <- lognormal(0,1)
lambda       <- gamma(1, 1)
coefs        <- laplace(0, 1/lambda, dim = ncol(X))
group_sd     <- lognormal(0,1)
group_effect <- normal(0, group_sd, dim = length(unique(group_id)))

# operation
mu    <- int + X %*% coefs + group_effect[group_id]

# likelihood
distribution(Y) <- normal(mu, sd)

# define the model
m <- model(int, coefs, sd, group_effect, group_sd, lambda)

# plot out the model 
plot(m)

# sampling
draws <- mcmc(m, n_samples = 1000)

# model summary
summary(draws)




####################################################################################
### NOW WITH STEVE'S REAL DATA

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Grab and set up the data

# note if running on docker, assuming you have linked a volume with the data in it to ./test
# steve's dataset, contains object ri_corr_select_centered_scaled
load('./test/RI covariates - selected for low missingness - NAs omitted - 20190327.RData')

# rename ri_corr_select_centered_scaled for ease of use
d <- ri_corr_select_centered_scaled

# make a design matrix with the covariates
X <- as.matrix(d[c("adolescent_fert", "controlofcorruption", "fertility",              
                   "governmenteffectiveness", "ln_gdp", "ln_healthexp", "politicalstability",     
                   "regulatoryquality", "ruleoflaw", "rural_pop", "voiceandaccountability")])

# treat year seperately (may or may not choose to include)
year <- d$year

# grab a vector of country IDs
country_id_tab <- data.frame(country_id=1:length(unique(d$iso)), iso=unique(d$iso))
d              <- merge(d, country_id_tab, by='iso')
country_id     <- d$country_id

# transform outcome
Y <- qlogis(d$wuenic_dtp3)

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Fit the model 
# priors
int          <- normal(0, 10)
sd           <- lognormal(0,1)
lambda       <- gamma(1, 1)
coefs        <- laplace(0, 1/lambda, dim = ncol(X))
group_sd     <- lognormal(0,1)
group_effect <- normal(0, group_sd, dim = length(unique(country_id)))

# operation
mu    <- int + X %*% coefs + group_effect[country_id]

# likelihood
distribution(Y) <- normal(mu, sd)

# define the model
m <- model(int, coefs, sd, group_effect, group_sd, lambda)

# plot out the model 
plot(m)

# sampling
draws <- mcmc(m, n_samples = 1000)

# model summary
summary(draws)$statistics








