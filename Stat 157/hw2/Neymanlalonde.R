library(Matching)
data(lalonde)
head(lalonde)

z = lalonde$treat
y = lalonde$re78

## Neymanian inference
n1= sum(z)
n0= length(z) - n1
tauhat = mean(y[z==1]) - mean(y[z==0])
vhat   = var(y[z==1])/n1 + var(y[z==0])/n0
sehat  = sqrt(vhat)
tauhat
sehat

## OLS
olsfit = lm(y ~ z)
summary(olsfit)$coef[2, 1: 2]
library(car)
sqrt(hccm(olsfit, type = "hc2")[2, 2])
sqrt(hccm(olsfit, type = "hc3")[2, 2])

