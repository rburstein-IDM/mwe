# MWE random slopes same 


# setup
library(lme4)
set.seed(54321)

# number of groups
ng <- 75 

# set parameters
b0  <- 1 # global intercept
b0g <- rnorm(ng,0,20) # random intercepts
b1  <- 4 # global slope
b1g <- rnorm(ng,0,12) # random slopes


## Sim data
# group ID (uneven number)
groupid <- rep(1:ng,rpois(ng,200)) # uneven number of groups like in steves data
N       <- length(groupid)

# indep variable
X     <- runif(N,-10,10)

# response
Y     <- rep(NA,N)
for(i in 1:N) {
  g    <- groupid[i]
  Y[i] <- b0 + b0g[g] + X[i]*(b1+b1g[g]) + rnorm(1,0,3)
}

# plot data
plot(X,Y,col=groupid)

# fit models (w/ and w/o random effects (intercepts and slopes))
m1 <- lmer(Y ~ 1 + X + (1 + X | groupid) )
m2 <- lm(Y ~ 1 + X)

# summarize, see different overall effects of X
summary(m1) 
summary(m2) # both the overall intercept and slope differ

# confirm those random slopes are actually mean zero
mean(coef(m1)$groupid$X) # nope
sum(coef(m1)$groupid$X) # just checking: its not a sum to zero constraint either

# simple diagnostics
hist(coef(m1)$groupid$X, breaks=20) # hist of random slopes
plot(b1g,coef(m1)$groupid$X) # actual versus predicted random slopes






### Another example, 3 very distinct groups to easily visualize
x1 <- runif(100,0,10)
y1 <- 100 + x1*5 + rnorm(100,0,1)

x2 <- runif(100,5,15)
y2 <- 5 + x2*-5 + rnorm(100,0,1)

x3 <- runif(100,12,20)
y3 <- 30 + x3*-5 + rnorm(100,0,1)

x <- c(x1,x2,x3)
y <- c(y1,y2,y3)
g <- rep(1:3,each=100)

# mods
m1 <- lm(y~1+x)
m2 <- lmer(y~1+x+(1+x|g))

#plot
plot(x,y,col=g)

# add lm line to plot
lines(x,summary(m1)$coefficients[1,1]+summary(m1)$coefficients[2,1]*x)

# add global lmer line to plot
lines(x,summary(m2)$coefficients[1,1]+summary(m2)$coefficients[2,1]*x,col='blue')










