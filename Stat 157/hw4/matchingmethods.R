library("car")
library("Matching")

## experimental data
data("lalonde")
head(lalonde)

y = lalonde$re78
z = lalonde$treat
x = as.matrix(lalonde[, c("age", "educ", "black",
                          "hisp", "married", "nodegr",
                          "re74", "re75")])

## analysis the randomized experiment
neymanols = lm(y ~ z)
round(summary(neymanols)$coef[2, ], 4)
sqrt(hccm(neymanols)[2, 2])

xc = scale(x)
linols = lm(y ~ z + xc + z*xc)
round(summary(linols)$coef[2, ], 4)
sqrt(hccm(linols)[2, 2])


## analysis as if it is an observational study
## Without adjustment 
matchest = Match(Y = y, Tr = z, X = x)
summary(matchest)

matchest.adj = Match(Y = y, Tr = z, X = x, BiasAdjust = TRUE)
summary(matchest.adj)

 


## observational data
dat <- read.table("cps1re74.csv",header=T)
# unemployed
dat$u74 <- as.numeric(dat$re74==0)
dat$u75 <- as.numeric(dat$re75==0)

head(dat)

y = dat$re78
z = dat$treat
x = as.matrix(dat[, c("age", "educ", "black",
                          "hispan", "married", "nodegree",
                          "re74", "re75", "u74", "u75")])

## analyze as if it is from a randomized experiment
neymanols = lm(y ~ z)
round(summary(neymanols)$coef[2, ], 4)
sqrt(hccm(neymanols)[2, 2])

xc = scale(x)
linols = lm(y ~ z + xc + z*xc)
round(summary(linols)$coef[2, ], 4)
sqrt(hccm(linols)[2, 2])

## analyze the observational study
## Without adjustment 
matchest = Match(Y = y, Tr = z, X = x)
summary(matchest)

matchest.adj = Match(Y = y, Tr = z, X = x, BiasAdjust = TRUE)
summary(matchest.adj)

matchest.adj$index.treated
matchest.adj$index.control

