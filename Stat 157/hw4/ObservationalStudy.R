ObsCausal.est = function(z, y, x, out.family = gaussian, 
                         truncpscore = c(0, 1))
{
     ## fitted propensity score
     pscore   = glm(z ~ x, family = binomial)$fitted.values
     pscore   = pmax(truncpscore[1], pmin(truncpscore[2], pscore))
     
     ## fitted potential outcomes
     outcome1 = glm(y ~ x, weights = z, 
                    family = out.family)$fitted.values
     outcome0 = glm(y ~ x, weights = (1 - z), 
                    family = out.family)$fitted.values
     
     ## regression imputation estimator
     ace.reg  = mean(outcome1 - outcome0) 
     ## propensity score weighting estimator
     ace.ipw0 = mean(z*y/pscore - (1 - z)*y/(1 - pscore))
     ace.ipw  = mean(z*y/pscore)/mean(z/pscore) - 
                   mean((1 - z)*y/(1 - pscore))/mean((1 - z)/(1 - pscore))
     ## doubly robust estimator
     res1     = y - outcome1
     res0     = y - outcome0
     ace.dr   = ace.reg + mean(z*res1/pscore - (1 - z)*res0/(1 - pscore))

     return(c(ace.reg, ace.ipw0, ace.ipw, ace.dr))     
}


ObsCausal = function(z, y, x, n.boot = 10^2,
                     out.family = gaussian, truncpscore = c(0, 1))
{
     point.est  = ObsCausal.est(z, y, x, out.family, truncpscore)
     
     ## nonparametric bootstrap
     n.sample   = length(z)
     x          = as.matrix(x)
     boot.est   = replicate(n.boot, 
                  {id.boot = sample(1:n.sample, n.sample, replace = TRUE)
                   ObsCausal.est(z[id.boot], y[id.boot], x[id.boot, ], 
                                 out.family, truncpscore)})

     boot.se    = apply(boot.est, 1, sd)
     
     res        = rbind(point.est, boot.se)
     rownames(res) = c("est", "se")
     colnames(res) = c("reg", "HT", "Hajek", "DR")
     
     return(res)
}


n.sim   = 100
ACE0    = rep(0, n.sim)
ACE     = matrix(0, 4, n.sim)
SEboot  = matrix(0, 4, n.sim)
n       = 200

## simulation with correct models
for(r in 1:n.sim)
{
  x       = matrix(rnorm(n*2), n, 2)
  x1      = cbind(1, x)
  beta.z  = c(0, 1, 1)
  pscore  = 1/(1 + exp(- as.vector(x1%*%beta.z)))
  z       = rbinom(n, 1, pscore)
  beta.y1 = c(1, 2, 1)
  beta.y0 = c(1, 2, 1)
  y1      = rnorm(n, x1%*%beta.y1)
  y0      = rnorm(n, x1%*%beta.y0)
  y       = z*y1 + (1 - z)*y0
  ACE0[r] = mean(y1 - y0)
  
  causaleffect = ObsCausal(z, y, x)
  ACE[, r]     = causaleffect[1, ]
  SEboot[, r]  = causaleffect[2, ]
}

apply(ACE, 1, mean) - mean(ACE0)
apply(ACE, 1, sd)
apply(SEboot, 1, mean)



## simulation with an incorrect propensity score model
for(r in 1:n.sim)
{
  x       = matrix(rnorm(n*2), n, 2)
  x1      = cbind(1, x, exp(x))
  beta.z  = c(-1, 0, 0, 1, -1)
  pscore  = 1/(1 + exp(- as.vector(x1%*%beta.z)))
  z       = rbinom(n, 1, pscore)
  beta.y1 = c(1, 2, 1, 0, 0)
  beta.y0 = c(1, 1, 1, 0, 0)
  y1      = rnorm(n, x1%*%beta.y1)
  y0      = rnorm(n, x1%*%beta.y0)
  y       = z*y1 + (1 - z)*y0
  ACE0[r] = mean(y1 - y0)
  
  causaleffect = ObsCausal(z, y, x)
  ACE[, r]     = causaleffect[1, ]
  SEboot[, r]  = causaleffect[2, ]
}

apply(ACE, 1, mean) - mean(ACE0)
apply(ACE, 1, sd)
apply(SEboot, 1, mean)



## simulation with an incorrect outcome model
for(r in 1:n.sim)
{
  x       = matrix(rnorm(n*2), n, 2)
  x1      = cbind(1, x, exp(x))
  beta.z  = c(0, 1, 1, 0, 0)
  pscore  = 1/(1 + exp(- as.vector(x1%*%beta.z)))
  z       = rbinom(n, 1, pscore)
  beta.y1 = c(1, 0, 0, 0.2, -0.1)
  beta.y0 = c(1, 0, 0, -0.2, 0.1)
  y1      = rnorm(n, x1%*%beta.y1)
  y0      = rnorm(n, x1%*%beta.y0)
  y       = z*y1 + (1 - z)*y0
  ACE0[r] = mean(y1 - y0)
  
  causaleffect = ObsCausal(z, y, x)
  ACE[, r]     = causaleffect[1, ]
  SEboot[, r]  = causaleffect[2, ]
}

apply(ACE, 1, mean) - mean(ACE0)
apply(ACE, 1, sd)
apply(SEboot, 1, mean)



## simulation without correct models
for(r in 1:n.sim)
{
  x       = matrix(rnorm(n*2), n, 2)
  x1      = cbind(1, x, exp(x))
  beta.z  = c(-1, 0, 0, 1, -1)
  pscore  = 1/(1 + exp(- as.vector(x1%*%beta.z)))
  z       = rbinom(n, 1, pscore)
  beta.y1 = c(1, 0, 0, 2, -1)
  beta.y0 = c(1, 0, 0, -2, 1)
  y1      = rnorm(n, x1%*%beta.y1)
  y0      = rnorm(n, x1%*%beta.y0)
  y       = z*y1 + (1 - z)*y0
  ACE0[r] = mean(y1 - y0)
  
  causaleffect = ObsCausal(z, y, x)
  ACE[, r]     = causaleffect[1, ]
  SEboot[, r]  = causaleffect[2, ]
}

apply(ACE, 1, mean) - mean(ACE0)
apply(ACE, 1, sd)
apply(SEboot, 1, mean)


## an application: the karolinska dataset
karolinska = read.table("karolinska.txt", header = TRUE)
head(karolinska)

z = karolinska$HighVolDiagHosp
y = 1 - (karolinska$YearsSurvivingAfterDiagnosis == 1)
x = as.matrix(karolinska[, c(3, 4, 5)])

causaleffects = ObsCausal(z, y, x, 
                          out.family = binomial, n.boot = 10^3)
causaleffects
rbind(causaleffects[1, ] - 1.96*causaleffects[2, ],
      causaleffects[1, ] + 1.96*causaleffects[2, ])

## checking the data
pscore   = glm(z ~ x, family = binomial)$fitted.values
hist(pscore[z==1], col="grey", border = NA, freq = FALSE,
     ylim = c(0, 4), main = "", breaks = 10,
     xlab = expression(hat(e)(X)), ylab = "")
hist(pscore[z==0], add= T, freq = FALSE, breaks = 10)



## an application: the NHANES BMI dataset
library(ATE)
data(nhanes_bmi)
z = nhanes_bmi$School_meal
y = nhanes_bmi$BMI
x = as.matrix(nhanes_bmi[, -c(1, 2)])

causaleffects = ObsCausal(z, y, x, n.boot = 10^3)
causaleffects
rbind(causaleffects[1, ] - 1.96*causaleffects[2, ],
      causaleffects[1, ] + 1.96*causaleffects[2, ])

## checking the data
pscore   = glm(z ~ x, family = binomial)$fitted.values
hist(pscore[z==1], col="grey", border = NA, freq = FALSE,
     ylim = c(0, 4.5), breaks = 30, main = "",
     xlab = expression(hat(e)(X)), ylab = "")
hist(pscore[z==0], add= T, freq = FALSE, breaks = 30)

## truncated propensity score
causaleffects = ObsCausal(z, y, x, n.boot = 10^3,
                          truncpscore = c(0.1, 0.9))
causaleffects
rbind(causaleffects[1, ] - 1.96*causaleffects[2, ],
      causaleffects[1, ] + 1.96*causaleffects[2, ])

