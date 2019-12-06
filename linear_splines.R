# Example code to send to Sam Dolan

library(splines)

# simulate super simple data 
x <- 1:99
b <- rep(c(-1,2,0), each = 33) # betas, assuming different linear functions from 0-33, 34-66, 67-99
y <- 2 + x*b + rnorm(99,0,4)

# notice the way I set it up there are potential discontinuities in the function we want to capture in the model
plot(x,y)

# bs just makes bases functions, what would these look like? 
xbs <- bs(x, degree = 1, knots = c(33,66))

# Note that the splines are additive where they overlap, this means two things:
# 1. the function is assumed to be continuous
# 2. the coefficients will not be very interpretable
# BONUS 3. Note that bs() rescales everything from 0 to 1, also making interepretation of effect sizes different. 
plot(1:99,xbs[,1], type = 'l', col = 'red')
lines(1:99,xbs[,2], col = 'blue')
lines(1:99,xbs[,3], col = 'green')

# we can fit a model, and predict, and show this model doesnt work well
m <- lm(y ~ bs(x, degree = 1, knots = c(33,66)))
summary(m)
p <- predict(m)
plot(x,y); lines(x,p, col='red')


# Solution: I think for the features you want, it may be easier to make your own bespoke linear functions. See below
x1 <- x * x %in%  1:33
x2 <- x * x %in% 34:66
x3 <- x * x %in% 67:99

# coefficents now match what we simulated, and the predictions are correct
m2 <- lm(y ~ x1 + x2 + x3)
summary(m2)
p2 <- predict(m2)
plot(x,y); lines(x,p2, col='red')

