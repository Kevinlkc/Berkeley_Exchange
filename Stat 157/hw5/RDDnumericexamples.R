library(rdd)
library(rdrobust)
library(rddtools)

## RDD numerical examples
n   = 500
x   = rnorm(n)
y0  = x + rnorm(n, 0, 0.5)
y1  = y0 + 5
z   = (x>=0)
y   = z*y1 + (1-z)*y0
plot(y0 ~ x, col = "grey", pch = 19, cex = 0.1,
     ylim = c(min(y), max(y)),
     xlab = "X", ylab = "Y")
points(y1 ~ x, col = "grey", pch = 19, cex = 0.1)
points(y ~ x, col = "black", pch = 19, cex = 0.1)
abline(v = 0, lty = 2)

plot(rdd_data(x=x, y=y,cutpoint=0), 
     xlab = "X", ylab = "Y", cex = 0.3)

RDDest = rdrobust(y, x)
cbind(RDDest$coef, RDDest$ci)

Greg = lm(y ~ z + x + z*x)
cbind(coef(Greg)[2], confint(Greg, 'zTRUE'))


y0  = x + rnorm(n, 0, 0.5)
y1  = y0 + 1
z   = (x>=0)
y   = z*y1 + (1-z)*y0
plot(y0 ~ x, col = "grey", pch = 19, cex = 0.1,
     ylim = c(min(y), max(y)))
points(y1 ~ x, col = "grey", pch = 19, cex = 0.1)
points(y ~ x, col = "black", pch = 19, cex = 0.1)
abline(v = 0, lty = 2)

plot(rdd_data(x=x, y=y,cutpoint=0), 
     xlab = "X", ylab = "Y", cex = 0.3)

RDDest = rdrobust(y, x)
cbind(RDDest$coef, RDDest$ci)

Greg = lm(y ~ z + x + z*x)
cbind(coef(Greg)[2], confint(Greg, 'zTRUE'))



y0  = 2*x + rnorm(n, 0, 0.5)
y1  = 5 + 0.5*x + rnorm(n, 0, 0.5)
z   = (x>=0)
y   = z*y1 + (1-z)*y0
plot(y0 ~ x, col = "grey", pch = 19, cex = 0.1,
     ylim = c(min(y), max(y)))
points(y1 ~ x, col = "grey", pch = 19, cex = 0.1)
points(y ~ x, col = "black", pch = 19, cex = 0.1)
abline(v = 0, lty = 2)

plot(rdd_data(x=x, y=y,cutpoint=0), 
     xlab = "X", ylab = "Y", cex = 0.3)

RDDest = rdrobust(y, x)
cbind(RDDest$coef, RDDest$ci)

Greg = lm(y ~ z + x + z*x)
cbind(coef(Greg)[2], confint(Greg, 'zTRUE'))


y0  = 2*x + rnorm(n, 0, 0.5)
y1  = 1 + 0.5*x + rnorm(n, 0, 0.5)
z   = (x>=0)
y   = z*y1 + (1-z)*y0
plot(y0 ~ x, col = "grey", pch = 19, cex = 0.1,
     ylim = c(min(y), max(y)))
points(y1 ~ x, col = "grey", pch = 19, cex = 0.1)
points(y ~ x, col = "black", pch = 19, cex = 0.1)
abline(v = 0, lty = 2)

plot(rdd_data(x=x, y=y,cutpoint=0), 
     xlab = "X", ylab = "Y", cex = 0.3)

RDDest = rdrobust(y, x)
cbind(RDDest$coef, RDDest$ci)

Greg = lm(y ~ z + x + z*x)
cbind(coef(Greg)[2], confint(Greg, 'zTRUE'))


## simulation under nonlinear models - hw

## data analysis
data(house)
head(house)

plot(y ~ x, data = house, pch = 19, cex = 0.3)
rdd_house = rdd_data(x=x, y=y, data=house, cutpoint=0)
summary(rdd_house)
plot(rdd_house, cex = 0.3,
     xlab = "X", ylab = "Y", main = "data: house")

RDDest = rdrobust(house$y, house$x)
cbind(RDDest$coef, RDDest$ci)

house$z = (house$x >= 0)
Greg = lm(y ~ z + x + z*x, data = house)
cbind(coef(Greg)[2], confint(Greg, 'zTRUE'))

x0.subset = seq(0.05, 1, 0.05)
local.lm = sapply(x0.subset,
                  function(x0){
                    Greg = lm(y ~ z + x + z*x, data = house,
                              subset = (abs(x)<=x0))
                    cbind(coef(Greg)[2], confint(Greg, 'zTRUE'))
                  })

plot(local.lm[1, ] ~ x0.subset, type = "b",
     pch = 19, cex = 0.5,
     ylim = range(local.lm),
     xlab = "x0",
     ylab = "point and interval estimates",
     main = "subset linear regression: |X|<x0")
lines(local.lm[2, ] ~ x0.subset, type = "b",
      pch = 19, cex = 0.5) 
lines(local.lm[3, ] ~ x0.subset, type = "b",
      pch = 19, cex = 0.5)