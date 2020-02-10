## covariate adjustment in IV analysis
IV_Lin = function(Z, D, Y, X)
{
  X     = scale(as.matrix(X))
  tau_D = lm(D ~ Z + X + Z*X)$coef[2]
  tau_Y = lm(Y ~ Z + X + Z*X)$coef[2]
  names(tau_D) = NULL
  names(tau_Y) = NULL
  CACE  = tau_Y/tau_D
  
  return(list(tau_D = tau_D, tau_Y = tau_Y,
              CACE  = CACE))
}

## IV_adj se via the delta method
IV_Lin_delta = function(Z, D, Y, X)
{
  X      = scale(as.matrix(X))
  est    = IV_Lin(Z, D, Y, X)
  
  betaY1 = lm(Y ~ X, subset = (Z == 1))$coef[-1]
  betaY0 = lm(Y ~ X, subset = (Z == 0))$coef[-1]
  betaD1 = lm(D ~ X, subset = (Z == 1))$coef[-1]
  betaD0 = lm(D ~ X, subset = (Z == 0))$coef[-1]
  
  AdjustedY1   = Y - X%*%betaY1 - 
                     (D - X%*%betaD1)*est$CACE
  AdjustedY0   = Y - X%*%betaY0 - 
                     (D - X%*%betaD0)*est$CACE
  VarAdj       = var(AdjustedY1[Z==1])/sum(Z) + 
                     var(AdjustedY0[Z==0])/sum(1 - Z)
  
  return(sqrt(VarAdj)/abs(est$tau_D))
}

##IV_adj se via the bootstrap
IV_Lin_bootstrap = function(Z, D, Y, X, n.boot = 200)
{
  X         = scale(as.matrix(X))
  CACEboot  = replicate(n.boot,
                        {
                          bindex = sample(1:length(Z), replace = TRUE)
                          IV_Lin(Z[bindex], D[bindex], Y[bindex], X[bindex])$CACE
                        })
  
  return(sqrt(var(CACEboot)))
}


## strong IV simulation
MC       = 200
CACEest  = rep(0, MC)
CACEse   = rep(0, MC)
CACEbse  = rep(0, MC)
n        = 200
for(i in 1:MC)
{
  X  = matrix(rnorm(n*2), n, 2)
  ## "c" n/2; "a" n/4; "n" n/4
  D0 = c(rep(0, n/2), rep(1, n/4), rep(0, n/4))
  D1 = c(rep(1, n/2), rep(1, n/4), rep(0, n/4))
  Y0 = c(rnorm(n/2, 1), rnorm(n/4, 0), rnorm(n/4, 2)) + 
            X%*%c(1, -1)
  Y1 = Y0
  Y1[1:(n/2)] = rnorm(n/2, 3) + X%*%c(-1, 1)
  Z  = rbinom(n, 1, 0.5)
  D  = Z*D1 + (1 - Z)*D0
  Y  = Z*Y1 + (1 - Z)*Y0
  
  CACEest[i]   = IV_Lin(Z, D, Y, X)$CACE
  CACEse[i]    = IV_Lin_delta(Z, D, Y, X)
  CACEbse[i]   = IV_Lin_bootstrap(Z, D, Y, X)
  
}

mean(CACEest)
sd(CACEest)
mean(CACEse)
mean(CACEbse)

hist(CACEest, breaks = 10, freq = FALSE,
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
  X  = matrix(rnorm(n*2), n, 2)
  ## "c" n*2/5; "a" n*2/5; "n" n/5
  D0 = c(rep(0, n*2/5), rep(1, n*2/5), rep(0, n/5))
  D1 = c(rep(1, n*2/5), rep(1, n*2/5), rep(0, n/5))
  Y0 = c(rnorm(n*2/5, 1), rnorm(n*2/5, 0), rnorm(n/5, 2)) + 
           X%*%c(1, -1)
  Y1 = Y0
  Y1[1:(n*2/5)] = rnorm(n*2/5, 3) + X%*%c(-1, 1)
  Z  = rbinom(n, 1, 0.5)
  D  = Z*D1 + (1 - Z)*D0
  Y  = Z*Y1 + (1 - Z)*Y0
  
  CACEest[i]   = IV_Lin(Z, D, Y, X)$CACE
  CACEse[i]    = IV_Lin_delta(Z, D, Y, X)
  CACEbse[i]   = IV_Lin_bootstrap(Z, D, Y, X)
  
}

mean(CACEest)
sd(CACEest)
mean(CACEse)
mean(CACEbse)

hist(CACEest, breaks = 20, freq = FALSE,
     xlab = expression(hat(tau)[c]), ylab = "",
     main = "weak IV",      
     border = FALSE, col = "grey")
abline(v = 2)
x = seq(0, 10, 0.01)
y = dnorm(x, mean(CACEest), sd(CACEest))
lines(y ~ x)


## weak IV
for(i in 1:MC)
{
  X  = matrix(rnorm(n*2), n, 2)
  ## "c" n/10; "a" n*2/5; "n" n/2
  D0 = c(rep(0, n/10), rep(1, n*2/5), rep(0, n/2))
  D1 = c(rep(1, n/10), rep(1, n*2/5), rep(0, n/2))
  Y0 = c(rnorm(n/10, 1), rnorm(n*2/5, 0), rnorm(n/2, 2)) + 
    X%*%c(1, -1)
  Y1 = Y0
  Y1[1:(n/10)] = rnorm(n/10, 3) + X%*%c(-1, 1)
  Z  = rbinom(n, 1, 0.5)
  D  = Z*D1 + (1 - Z)*D0
  Y  = Z*Y1 + (1 - Z)*Y0
  
  CACEest[i]   = IV_Lin(Z, D, Y, X)$CACE
  CACEse[i]    = IV_Lin_delta(Z, D, Y, X)
  CACEbse[i]   = IV_Lin_bootstrap(Z, D, Y, X)
  
}

mean(CACEest)
sd(CACEest)
mean(CACEse)
mean(CACEbse)

hist(CACEest, breaks = 15, freq = FALSE,
     xlab = expression(hat(tau)[c]), ylab = "",
     main = "weak IV",      
     border = FALSE, col = "grey")
abline(v = 2)
x = seq(-6000, 10, 1)
y = dnorm(x, mean(CACEest), sd(CACEest))
lines(y ~ x)


## jobs data
jobsdata = read.csv("jobsdata.csv")
getX     = lm(treat ~ sex + age + marital 
                      + nonwhite + educ + income,
              data = jobsdata)
X = model.matrix(getX)[, -1]
Z = jobsdata$treat
D = jobsdata$comply
Y = jobsdata$job_seek
est = IV_Lin(Z, D, Y, X)$CACE
dse = IV_Lin_delta(Z, D, Y, X)
bse = IV_Lin_bootstrap(Z, D, Y, X)
est
c(est - 1.96*dse, est + 1.96*dse)
c(est - 1.96*bse, est + 1.96*bse)

