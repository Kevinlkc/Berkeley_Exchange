## Imbens and Rubin book: matched pair data
## television program aimed at improving reading skills for children
dataxy = c(12.9, 12.0, 54.6, 60.6,
           15.1, 12.3, 56.5, 55.5,
           16.8, 17.2, 75.2, 84.8,
           15.8, 18.9, 75.6, 101.9,
           13.9, 15.3, 55.3, 70.6,
           14.5, 16.6, 59.3, 78.4,
           17.0, 16.0, 87.0, 84.2,
           15.8, 20.1, 73.7, 108.6)
           
dataxy = matrix(dataxy, 8, 4,  byrow = TRUE)           

diffx = dataxy[, 2] - dataxy[, 1]
diffy = dataxy[, 4] - dataxy[, 3]

dataxy = cbind(dataxy, diffx, diffy)

rownames(dataxy) = 1:8
colnames(dataxy) = c("x.control", "x.treatment", 
                     "y.control", "y.treatment",
                     "diffx", "diffy")
                     
dataxy

## analysis without covariates
n      = dim(dataxy)[1]
tauhat = mean(dataxy[, "diffy"])
vhat   = var(dataxy[, "diffy"])/n
tauhat
sqrt(vhat)

## regression analysis
summary(lm(diffy ~ 1, 
           data = data.frame(dataxy)))$coef

summary(lm(diffy ~ diffx, 
           data = data.frame(dataxy)))$coef







