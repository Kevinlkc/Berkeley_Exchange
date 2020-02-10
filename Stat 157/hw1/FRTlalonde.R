library(Matching)
data(lalonde)
head(lalonde)

z = lalonde$treat
y = lalonde$re78
par(mfrow = c(1, 1))
hist(y[z == 0], col = "grey", border = FALSE, 
     freq = FALSE, breaks = 30, 
     xlab = "real earnings in 1978", 
     ylab = "density",
     main = "", yaxt = "n")
hist(y[z == 1], freq = FALSE, 
     breaks = 30, add = TRUE)
abline(v = mean(y[z == 1]))
abline(v = mean(y[z == 0]), lty = 2)

## approximate tests, easy implementation in R
asym.pv = c(t.test(y ~ z, var.equal = TRUE)$p.value, 
t.test(y ~ z, var.equal = FALSE)$p.value,
wilcox.test(y ~ z)$p.value,
ks.test(y[z == 1], y[z == 0])$p.value)

## Fisher randomization test
## sample size 
nn = length(z)
nn
## number of treated units
nn1 = sum(z)
nn1
## number of permutations
choose(nn, nn1)

## Monte Carlo version of the FRT
MC = 10^3
Tauhat   = rep(0, MC)
Student  = rep(0, MC)
Wilcox   = rep(0, MC)
Ks       = rep(0, MC)
for(mc in 1:MC)
{
       zperm       = sample(z)
       Tauhat[mc]  = t.test(y ~ zperm, var.equal = TRUE)$statistic 
       Student[mc] = t.test(y ~ zperm, var.equal = FALSE)$statistic
       Wilcox[mc]  = wilcox.test(y ~ zperm)$statistic 
       Ks[mc]      = ks.test(y[zperm == 1], y[zperm == 0])$statistic
}

par(mfrow = c(2, 2), mar = c(4,1,1,1), mgp = c(1.5, 0.5, 0))
hist(Tauhat, col = "grey", 
     border = FALSE, freq = FALSE,
     xlab = expression(t[eqvar]),
     ylab = "",
     main = "", yaxt = "n")
tauhat = t.test(y ~ z, var.equal = TRUE)$statistic
abline(v = tauhat)
mean(Tauhat <= tauhat)
xx = seq(-5, 5, 0.01)
yy = dnorm(xx)
lines(yy ~ xx)


hist(Student, col = "grey", 
     border = FALSE, freq = FALSE,
     xlab = expression(t[uneqvar]),
     ylab = "",
     main = "", yaxt = "n")
student = t.test(y ~ z, var.equal = FALSE)$statistic
abline(v = student)
mean(Student <= student)
lines(yy ~ xx)

hist(Wilcox, col = "grey", 
     border = FALSE, freq = FALSE,
     xlab = expression(W),
     ylab = "",
     main = "", yaxt = "n")
W = wilcox.test(y ~ z)$statistic
abline(v = W)
mean(Wilcox <= W)
mean.null = nn1*(nn - nn1 + 1)/2 
var.null  = nn1*(nn - nn1)*(nn + 1)/12
se.null   = sqrt(var.null)
xx        = seq(mean.null - 3*se.null, mean.null + 3*se.null, 1)
yy        = dnorm(xx, mean.null, se.null)
lines(yy ~ xx)

hist(Ks, col = "grey", 
     border = FALSE, freq = FALSE,
     xlab = expression(D),
     ylab = "",
     main = "", yaxt = "n")
D = ks.test(y[z == 1], y[z == 0])$statistic
abline(v = D)
mean(Ks >= D)



