## Roy Burstein
## Interval censored data simulation and model for adam

library(data.table)

# set params
nsims <- 10000


#############
# simulate underlying true failure times from a weibull distribution
truth <- data.table(sim           = 1:nsims,
                    true_inf_date = rbinom(nsims, 1, 0.4) * rweibull(nsims, shape = 3, scale = 30) )
truth[true_inf_date == 0, true_inf_date := 999] # set never infected to 999

# plot out the true values
hist(truth$true_inf_date[truth$true_inf_date!=999], breaks = 50)


# assume simply each simulant gets tested once every two years on the first day of the year ( can complicate this later)
d <- data.table(expand.grid(sim           = 1:sims,
                            test_year     = seq(15,60, by = 2)))
d <- merge(d, truth, by = 'sim')

# make a variable for whether or not they were infected as of that test date
d[, status := as.numeric(true_inf_date < test_year)]




#############
# define an appropriate likelihood
# note, we must account for zero inflation, left truncation, interval censoring, and right censoring 

