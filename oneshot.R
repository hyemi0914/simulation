###------------------------------------------------------###
###        R code for one-shot simulation study          ###
###------------------------------------------------------###
rm(list=ls())
set.seed(1)

## load R functions 
source("Gamma-SB_onedim-MCMC.R")

## simulation settings 
m <- 200
Ay <- rep(5, m)
By <- rep(5, m)
om <- 0.05
m_signal <- round(om * m)

mu <- 5    # location of null signals 
sig <- rep(5, m)
sig[1:m_signal] <- rgamma(m_signal, shape=mu*20, rate=2)
hist(sig)
y <- rgamma(m, shape = Ay, rate = By / sig)




## tuning parameters 
mc <- 5000
burn <- 2000
q <- 0.05
const <- NA
S <- 100


## SB prior 
mcmc.SB <- SB_prior(y, Ay, By, a_ta=1, b_ta=1, a_be=mu, b_be=1, a=2, b=0.5, mc=mc, burn=burn)[["La"]]    
est1 <- colMeans(mcmc.SB)


## IRB prior 
mcmc.IRB <- IRB_prior(y, Ay, By, a_ta=1, b_ta=1, a_be=mu, b_be=1, a=2, b=0.5, mc=mc, burn=burn)[["La"]]   
est2 <- colMeans(mcmc.IRB)


## global (gamma) prior
mcmc.G <- global_prior(y, Ay, By, a_ta=1, b_ta=1, a_be=mu, b_be=1, mc=mc, burn=burn)[["La"]]   
est3 <- colMeans(mcmc.G)


## DasGupta's method 
est4 <- DasGupta(y, Ay, By, const)


## VASH method 
library(vashr)
fit <- vash(y, df=10) 
est5 <- fit$sd.post



## Figure of shrinkage effect
pdf("shrink.pdf", height=8, width=16, pointsize=18)
par(mfcol=c(1,2))
# fig1 
plot(est1, y, xlab="point estimate", ylab="observed data", xlim=range(y), pch=4)
points(est2, y, pch=20)
abline(0, 1, lty=2)
points(est3, y, col=2, pch=4)
points(est4, y, col=3, pch=4)
points(est5, y, col=4, pch=4)
legend("bottomright", legend=c("SB", "ISB", "GL", "DG", "VS"), col=c(1,1:4), pch=c(4,1,4,4,4))
# fig2 
plot(est1, y, xlab="point estimate", ylab="observed data", pch=4, xlim=c(0, 10), ylim=c(0, 15))
points(est2, y, pch=20)
abline(0, 1, lty=2)
points(est3, y, col=2, pch=4)
points(est4, y, col=3, pch=4)
points(est5, y, col=4, pch=4)
abline(v=5)
legend("bottomright", legend=c("SB", "ISB", "GL", "DG", "VS"), col=c(1,1:4), pch=c(4,1,4,4,4))
dev.off()



## MAPE 
mean(abs(y-sig)/sig)
mean(abs(est1-sig)/sig)
mean(abs(est2-sig)/sig)
mean(abs(est3-sig)/sig)
mean(abs(est4-sig)/sig)
mean(abs(est5-sig)/sig)



## MAPE (large signals) 
mean((abs(y-sig)/sig)[1:m_signal])
mean((abs(est1-sig)/sig)[1:m_signal])
mean((abs(est2-sig)/sig)[1:m_signal])
mean((abs(est3-sig)/sig)[1:m_signal])
mean((abs(est4-sig)/sig)[1:m_signal])
mean((abs(est5-sig)/sig)[1:m_signal])



# RSMSE
sqrt( mean((y-sig)^2/sig^2) )
sqrt( mean((est1-sig)^2/sig^2) )
sqrt( mean((est2-sig)^2/sig^2) )
sqrt( mean((est3-sig)^2/sig^2) )
sqrt( mean((est4-sig)^2/sig^2) )
sqrt( mean((est5-sig)^2/sig^2) )
