## IV point estimator
IV_Wald = function(Z, D, Y)
{
       tau_D = mean(D[Z==1]) - mean(D[Z==0])
       tau_Y = mean(Y[Z==1]) - mean(Y[Z==0])
       CACE  = tau_Y/tau_D
       
       return(list(tau_D = tau_D, tau_Y = tau_Y,
                   CACE  = CACE))
}

## IV se via the delta method
IV_Wald_delta = function(Z, D, Y)
{
       est         = IV_Wald(Z, D, Y)
       AdjustedY   = Y - D*est$CACE
       VarAdj      = var(AdjustedY[Z==1])/sum(Z) + 
                          var(AdjustedY[Z==0])/sum(1 - Z)
       return(sqrt(VarAdj)/abs(est$tau_D))
}

##IV se via the bootstrap
IV_Wald_bootstrap = function(Z, D, Y, n.boot = 200)
{
       CACEboot  = replicate(n.boot,
                   {
                   bindex = sample(1:length(Z), replace = TRUE)
                   IV_Wald(Z[bindex], D[bindex], Y[bindex])$CACE
                   })
       
       return(sd(CACEboot))
}


## strong IV simulation
MC       = 500
CACEest  = rep(0, MC)
CACEse   = rep(0, MC)
CACEbse  = rep(0, MC)
n        = 200
for(i in 1:MC)
{
  ## "c" n/2; "a" n/4; "n" n/4
  D0 = c(rep(0, n/2), rep(1, n/4), rep(0, n/4))
  D1 = c(rep(1, n/2), rep(1, n/4), rep(0, n/4))
  Y0 = c(rnorm(n/2, 1), rnorm(n/4, 0), rnorm(n/4, 2))
  Y1 = Y0
  Y1[1:(n/2)] = rnorm(n/2, 3)
  Z  = rbinom(n, 1, 0.5)
  D  = Z*D1 + (1 - Z)*D0
  Y  = Z*Y1 + (1 - Z)*Y0
  
  CACEest[i]   = IV_Wald(Z, D, Y)$CACE
  CACEse[i]    = IV_Wald_delta(Z, D, Y)
  CACEbse[i]   = IV_Wald_bootstrap(Z, D, Y)
  
}

mean(CACEest)
sd(CACEest)
mean(CACEse)
mean(CACEbse)

hist(CACEest, breaks = 15, freq = FALSE,
     xlab = expression(hat(tau)[c]), ylab = "",
     main = "strong IV",      
     border = FALSE, col = "grey")
abline(v = 2)
x = seq(0, 5, 0.01)
y = dnorm(x, mean(CACEest), sd(CACEest))
lines(y ~ x)

## weak IV
for(i in 1:MC)
{
  ## "c" n/5; "a" n*2/5; "n" n*2/5
  D0 = c(rep(0, n/5), rep(1, n*2/5), rep(0, n*2/5))
  D1 = c(rep(1, n/5), rep(1, n*2/5), rep(0, n*2/5))
  Y0 = c(rnorm(n/5, 1), rnorm(n*2/5, 0), rnorm(n*2/5, 2))
  Y1 = Y0
  Y1[1:(n/5)] = rnorm(n/5, 3)
  Z  = rbinom(n, 1, 0.5)
  D  = Z*D1 + (1 - Z)*D0
  Y  = Z*Y1 + (1 - Z)*Y0
  
  CACEest[i]   = IV_Wald(Z, D, Y)$CACE
  CACEse[i]    = IV_Wald_delta(Z, D, Y)
  CACEbse[i]   = IV_Wald_bootstrap(Z, D, Y)
  
}

mean(CACEest)
sd(CACEest)
mean(CACEse)
mean(CACEbse)

hist(CACEest, breaks = 30, freq = FALSE,
     xlab = expression(hat(tau)[c]), ylab = "",
     main = "weak IV",      
     border = FALSE, col = "grey")
abline(v = 2)
x = seq(-40, 40, 0.01)
y = dnorm(x, mean(CACEest), sd(CACEest))
lines(y ~ x)


## weak IV
for(i in 1:MC)
{
  ## "c" n/10; "a" n*2/5; "n" n/2
  D0 = c(rep(0, n/10), rep(1, n*2/5), rep(0, n/2))
  D1 = c(rep(1, n/10), rep(1, n*2/5), rep(0, n/2))
  Y0 = c(rnorm(n/10, 1), rnorm(n*2/5, 0), rnorm(n/2, 2))
  Y1 = Y0
  Y1[1:(n/10)] = rnorm(n/10, 3)
  Z  = rbinom(n, 1, 0.5)
  D  = Z*D1 + (1 - Z)*D0
  Y  = Z*Y1 + (1 - Z)*Y0
  
  CACEest[i]   = IV_Wald(Z, D, Y)$CACE
  CACEse[i]    = IV_Wald_delta(Z, D, Y)
  CACEbse[i]   = IV_Wald_bootstrap(Z, D, Y)
  
}

mean(CACEest)
sd(CACEest)
mean(CACEse)
mean(CACEbse)

hist(CACEest, breaks = 30, freq = FALSE,
     xlab = expression(hat(tau)[c]), ylab = "",
     main = "weak IV",      
     border = FALSE, col = "grey")
abline(v = 2)
x = seq(-400, 600, 1)
y = dnorm(x, mean(CACEest), sd(CACEest))
lines(y ~ x)

## jobs data in "mediation package"
## one sided noncompliance: D(0) = 0
jobsdata = read.csv("jobsdata.csv")
Z = jobsdata$treat
D = jobsdata$comply
Y = jobsdata$job_seek
est = IV_Wald(Z, D, Y)$CACE
dse = IV_Wald_delta(Z, D, Y)
bse = IV_Wald_bootstrap(Z, D, Y)
est
c(est - 1.96*dse, est + 1.96*dse)
c(est - 1.96*bse, est + 1.96*bse)

