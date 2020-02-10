ATT.est = function(z, y, x, out.family = gaussian, Utruncpscore = 1)
{
  ## sample size
  nn  = length(z)
  nn1 = sum(z)
  
  ## fitted propensity score
  pscore   = glm(z ~ x, family = binomial)$fitted.values
  pscore   = pmin(Utruncpscore, pscore)
  odds.pscore = pscore/(1 - pscore)
  
  ## fitted potential outcomes
  outcome0 = glm(y ~ x, weights = (1 - z), 
                 family = out.family)$fitted.values
  
  ## regression imputation estimator
  ace.reg0 = lm(y ~ z + x)$coef[2]
  ace.reg  = mean(y[z==1]) - mean(outcome0[z==1]) 
  ## propensity score weighting estimator
  ace.ipw0 = mean(y[z==1]) - 
                mean(odds.pscore*(1 - z)*y)*nn/nn1
  ace.ipw  = mean(y[z==1]) - 
                mean(odds.pscore*(1 - z)*y)/mean(odds.pscore*(1 - z))
  ## doubly robust estimator
  res0     = y - outcome0
  ace.dr   = ace.reg - mean(odds.pscore*(1 - z)*res0)*nn/nn1
  
  return(c(ace.reg0, ace.reg, ace.ipw0, ace.ipw, ace.dr))     
}


ObsCausal.ATT = function(z, y, x, n.boot = 10^2,
                         out.family = gaussian, Utruncpscore = 1)
{
  point.est  = ATT.est(z, y, x, out.family, Utruncpscore)
  
  ## nonparametric bootstrap
  n.sample   = length(z)
  x          = as.matrix(x)
  boot.est   = replicate(n.boot, 
                         {id.boot = sample(1:n.sample, n.sample, replace = TRUE)
                          ATT.est(z[id.boot], y[id.boot], x[id.boot, ], 
                                  out.family, Utruncpscore)})
  
  boot.se    = apply(boot.est, 1, sd)
  
  res        = rbind(point.est, boot.se)
  rownames(res) = c("est", "se")
  colnames(res) = c("reg0", "reg", "HT", "Hajek", "DR")
  
  return(res)
}



n.sim   = 100
ATT0    = rep(0, n.sim)
ATT     = matrix(0, 5, n.sim)
SEboot  = matrix(0, 5, n.sim)
n       = 200

## simulation with correct models
## parallel y1 and y0
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
  ATT0[r] = mean(y1[z==1]) - mean(y0[z==1])
  
  causaleffect = ObsCausal.ATT(z, y, x)
  ATT[, r]     = causaleffect[1, ]
  SEboot[, r]  = causaleffect[2, ]
}

apply(ATT, 1, mean) - mean(ATT0)
apply(ATT, 1, sd)
apply(SEboot, 1, mean)


## nonparallel y1 and y0
for(r in 1:n.sim)
{
  x       = matrix(rnorm(n*2), n, 2)
  x1      = cbind(1, x)
  beta.z  = c(0, 1, 1)
  pscore  = 1/(1 + exp(- as.vector(x1%*%beta.z)))
  z       = rbinom(n, 1, pscore)
  beta.y1 = c(1, 2, 1)
  beta.y0 = c(1, 1, 1)
  y1      = rnorm(n, x1%*%beta.y1)
  y0      = rnorm(n, x1%*%beta.y0)
  y       = z*y1 + (1 - z)*y0
  ATT0[r] = mean(y1[z==1]) - mean(y0[z==1])
  
  causaleffect = ObsCausal.ATT(z, y, x)
  ATT[, r]     = causaleffect[1, ]
  SEboot[, r]  = causaleffect[2, ]
}

apply(ATT, 1, mean) - mean(ATT0)
apply(ATT, 1, sd)
apply(SEboot, 1, mean)


## simulation with incorrect models - homework problem



## an application: the karolinska dataset
karolinska = read.table("karolinska.txt", header = TRUE)
head(karolinska)

z = karolinska$HighVolDiagHosp
y = 1 - (karolinska$YearsSurvivingAfterDiagnosis == 1)
x = as.matrix(karolinska[, c(3, 4, 5)])

ATTs = ObsCausal.ATT(z, y, x, 
                     out.family = binomial, n.boot = 10^3)
ATTs
rbind(ATTs[1, ] - 1.96*ATTs[2, ],
      ATTs[1, ] + 1.96*ATTs[2, ])

 

## an application: the NHANES BMI dataset
library(ATE)
data(nhanes_bmi)
z = nhanes_bmi$School_meal
y = nhanes_bmi$BMI
x = as.matrix(nhanes_bmi[, -c(1, 2)])

ATTs = ObsCausal.ATT(z, y, x, n.boot = 10^3)
ATTs
rbind(ATTs[1, ] - 1.96*ATTs[2, ],
      ATTs[1, ] + 1.96*ATTs[2, ])


## truncated propensity score
ATTs = ObsCausal.ATT(z, y, x, n.boot = 10^3, Utruncpscore = 0.9)
ATTs
rbind(ATTs[1, ] - 1.96*ATTs[2, ],
      ATTs[1, ] + 1.96*ATTs[2, ])

