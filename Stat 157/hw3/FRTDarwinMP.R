## enumerate all possible matched pair assignment 
## using binary representation
MP_enumerate = function(i, n.pairs = 15) 
{
 a = 2^((n.pairs-1):0)
 b = 2*a
 2*sapply(i-1, 
          function(x) 
            as.integer((x %% b)>=a)) - 1
}


## Darwin's data from Fisher's book
ytreat     = c(188, 96, 168, 176, 153, 
               172, 177, 163, 146, 173, 
               186, 168, 177, 184, 96)
ycontrol   = c(139, 163, 160, 160, 147, 
               149, 149, 122, 132, 144, 
               130, 144, 102, 124, 144)
difference = ytreat - ycontrol
n.pairs    = length(difference)
abs.diff   = abs(difference)
t.obs      = mean(difference)
t.ran      = sapply(1:2^15, 
                    function(x){ 
                      sum(MP_enumerate(x, 15)*abs.diff) 
                      })/n.pairs
pvalue     = mean(t.ran>=t.obs)

hist(t.ran, breaks = 50, col = "grey", border = NA,
     xlab = expression(hat(tau)), 
     ylab = "", yaxt = 'n', 
     main = "randomization distribution - Darwin's data")
abline(v = t.obs)
text(30, 400, 
     paste("p-value = ", round(pvalue, 3), sep = ""))



