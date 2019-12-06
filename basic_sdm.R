# File:    Customers.R
# Author:  Jim Duggan
# Summary: One-stock example of customers , with fraction increases/descreases

library(deSolve)
library(ggplot2)

# Setup simulation times and time step
START<-2015; FINISH<-2030; STEP<-0.5

# Create time vector
simtime <- seq(START, FINISH, by=STEP)

# Create stock and auxs
stocks  <- c(sCustomers=10000)
auxs    <- c(aGrowthFraction=0.08, aDeclineFraction=0.03)


# The Model function, takes 3 arguments from ode()
model <- function(time, stocks, auxs){
  with(as.list(c(stocks, auxs)),{
    
    # The inflow equation
    fRecruits <- sCustomers*aGrowthFraction
    
    # The outflow equation
    fLosses <- sCustomers*aDeclineFraction
    
    # The stock equation
    dC_dt <- fRecruits - fLosses
    
    # All the results for the time step
    ans <- list(c(dC_dt),
                Recruits=fRecruits, 
                Losses=fLosses,
                NetFlow=dC_dt,
                GF=aGrowthFraction,
                DF=aDeclineFraction)
  })
}


# Run simulation
o<-data.frame(ode(y=stocks, times=simtime, func = model, 
                  parms=auxs, method='euler'))

# Plot results
qplot(x=time,y=sCustomers,data=o) + geom_line()





######### 
## HF model from alonge

# Setup simulation times and time step
START<-0; FINISH<-75; STEP<-1

# Create time vector
simtime <- seq(START, FINISH, by=STEP)

stocks <- c(R = 0, Q = 0, V = -.25)
auxs   <- c(rR = 0.3, rRQ = 0.2, rQ = 0.3, rV = 0.1)

L <- function(x) { (1-exp(-2*x))/(1+exp(02*x))  }

# The Model function, takes 3 arguments from ode()
model2 <- function(time, stocks, auxs){
  with(as.list(c(stocks, auxs)),{

    dR_dt <- rR*L(V) - rRQ*R
    dQ_dt <- rRQ*R - rQ*L(V)
    dV_dt <- rV*L(Q) # need to delay quality somehow! (look up dede)
    
    # All the results for the time step
    ans <- list(c(dR_dt,dQ_dt,dV_dt))
  })
}
o2<-data.frame(ode(y=stocks, times=simtime, func = model2, 
                  parms=auxs))

o2 <- data.table(melt(o2, id = 'time'))

ggplot(o2, aes(x=time,y=value,color=variable)) + geom_line()



####
# LAGS
## the derivative function
derivs <- function(t, y, parms) {
  if (t < 1)
    dy <- -1
  else
    dy <- - lagvalue(t - 1)
  list(c(dy))
}

## initial values and times
yinit <- 1
times <- seq(0, 30, 0.1)

## solve the model
yout <- dede(y = yinit, times = times, func = derivs, parms = NULL)


## display, plot results
plot(yout, type = "l", lwd = 2, main = "dy/dt = -y(t-1)")




######### 
## HF model from alonge, but add lag delay in quality

# Setup simulation times and time step
START<-0; FINISH<-200; STEP<-1

# Create time vector
simtime <- seq(START, FINISH, by=STEP)

stocks <- c(R = 0, Q  = 0, V = -0.5)
auxs   <- c(rR = 0.3, rRQ = 0.2, rQ = 0.3, rV = 0.1, delay = 3)

L <- function(x) { (1-exp(-2*x))/(1+exp(02*x))  }

# The Model function, takes 3 arguments from ode()
model3 <- function(time, stocks, auxs){
  with(as.list(c(stocks, auxs)),{
    
    if (time >= delay)
      QQ <- lagvalue(time - delay, 2)
    else
      QQ <- lagvalue(1,2)
    
    dR_dt <- rR*L(V) - rRQ*R
    dQ_dt <- rRQ*R - rQ*L(V)
    dV_dt <- rV*L(QQ) # need to delay quality somehow! (look up dede)
    
    # All the results for the time step
    ans <- list(c(dR_dt,dQ_dt,dV_dt))
  })
}
o3<-data.frame(dede(y=stocks, times=simtime, func = model3, 
                   parms=auxs))

o3 <- data.table(melt(o3, id = 'time'))

ggplot(o3, aes(x=time,y=value,color=variable)) + geom_line()






######### 
## HF model from alonge, full with gaming for complete replication. 

# Setup simulation times and time step
START<-0; FINISH<-200; STEP<-1

# Create time vector
simtime <- seq(START, FINISH, by=STEP)

stocks <- c(R = 0, Q  = 0, V = -0.5)
auxs   <- c(rR        = 0.3, 
            rRQ       = 0.2,
            rQ        = 0.3, 
            rV        = 0.1, 
            Qdelay    = 3,
            Vdelay    = 1,
            P4Pfactor = 0,
            G_0       = 0,
            M_0       = 0)

L <- function(x) { (1-exp(-2*x))/(1+exp(-2*x))  }


# The Model function, takes 3 arguments from ode()
model4 <- function(time, stocks, auxs){
  with(as.list(c(stocks, auxs)),{
   
    if (time >= Qdelay)
      Qdelay <- lagvalue(time - Qdelay, 2)
    else
      Qdelay <- lagvalue(1,2)
    
    if (time >= Vdelay)
      Vdelay <- lagvalue(time - Vdelay, 3)
    else
      Vdelay <- lagvalue(1,2)
    
    # dependent variables
    P4P <- P4Pfactor * L(Vdelay) # add delay here
    G   <- G_0 + 1.1^P4P - 1
    M   <- M_0 + P4P - G

    # differential equations
    dR_dt <- rR*L(V) - rRQ*R*1.3^M + P4P
    
    if(V>0)
      dQ_dt <- rRQ*R*1.3^M - rQ*L(V)*(1.1^G)/(1.3^M)
    else 
      dQ_dt <- rRQ*R*1.3^M - rQ*L(V)*(1.3^G)/(1.1^M)
    
    if(Qdelay>0)
      dV_dt <- rV*L(Qdelay)*1.3^M
    else 
      dV_dt <- rV*L(Qdelay)/1.3^M
    
    # All the results for the time step
    ans <- list(c(dR_dt,dQ_dt,dV_dt))
  })
}

o4 <- data.frame(dede(y=stocks, times=simtime, func = model4,  parms=auxs))
o4 <- data.table(melt(o4, id = 'time'))

ggplot(o4[variable!='R'], aes(x=time,y=value,color=variable)) + geom_line(lwd=2) +
  ylim(-1,1) + xlim(0,100) + theme_bw() + 
  scale_color_manual(values=c("#3465eb", "#eb6134"), labels = c('Quality', "Volume"))

