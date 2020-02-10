## estimation
Neyman_SRE = function(z, y, x)
{
       xlevels = unique(x)
       K       = length(xlevels)
       PiK     = rep(0, K)
       TauK    = rep(0, K)
       varK    = rep(0, K)
       for(k in 1:K)
       {
             xk         = xlevels[k]
             zk         = z[x == xk]
             yk         = y[x == xk]
             PiK[k]     = length(zk)/length(z)
             TauK[k]    = mean(yk[zk==1]) - mean(yk[zk==0])
             varK[k]    = var(yk[zk==1])/sum(zk) + 
                               var(yk[zk==0])/sum(1 - zk)
       }
       
       return(c(sum(PiK*TauK), sum(PiK^2*varK)))
}

## pennsylvania re-employment bonus experiment
## description of the DATA: 
## Koenker and Xiao 2002 Econometrica 
## "Inference on the Quantile Regression Process" 

penndata = read.table("Penn46_ascii.txt")
head(penndata)

z = penndata$treatment
y = log(penndata$duration)
block = penndata$quarter
est = Neyman_SRE(z, y, block)
est[1]
sqrt(est[2])
