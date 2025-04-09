###------------------------------------------------------###
###        R code for one-shot simulation study          ###
###------------------------------------------------------###
rm(list=ls())
set.seed(1)

## load R functions 
source("Gamma-SB-MCMC.R")

## simulation settings
m <- 200
D <- 2
Ay <- rep(5, m)
By <- rep(5, m)
om <- 0.05
m_signal <- round(om * m)
mu <- 5    # location of null signals
c_values <- seq(0, 100, by = 20)
sim_results <- list()

for (c in c_values) {
  # d = 1
  sig1 <- rep(mu, m)
  slab_idx1 <- sample(1:m, m_signal)
  sig1[slab_idx1] <- rgamma(m_signal, shape = mu * 20, rate = 2)
  
  # d = 2
  sig2 <- rep(mu, m)
  S_shared <- round(m_signal * c / 100)
  S_new <- m_signal - S_shared
  shared_idx <- sample(slab_idx1, S_shared)
  spike_idx2 <- setdiff(1:m, slab_idx1)
  new_idx <- sample(spike_idx2, S_new)
  slab_idx2 <- c(shared_idx, new_idx)
  sig2[slab_idx2] <- rgamma(m_signal, shape = mu * 20, rate = 2)
  
  y1 <- rgamma(m, shape = Ay, rate = By / sig1)
  y2 <- rgamma(m, shape = Ay, rate = By / sig2)
  
  ## tuning parameters 
  mc <- 5000
  burn <- 2000
  q <- 0.05
  const <- NA
  S <- 100
  
  ## SB prior 
  mcmc.SB1 <- SB_prior(y1, Ay, By, a_ta = 1, b_ta = 1, a_be = mu, b_be = 1, a = 2, b = 0.5, mc = mc, burn = burn)[["La"]]
  mcmc.SB2 <- SB_prior(y2, Ay, By, a_ta = 1, b_ta = 1, a_be = mu, b_be = 1, a = 2, b = 0.5, mc = mc, burn = burn)[["La"]]
  est1 <- colMeans(mcmc.SB1)
  est2 <- colMeans(mcmc.SB2)
  
  sim_results[[paste0("c_", c)]] <- list(y1 = y1, y2 = y2, est1 = est1, est2 = est2, slab_idx1 = slab_idx1, slab_idx2 = slab_idx2)
}

# fig1
par(mfrow = c(3, 2))
for (c in c_values) {
  results <- sim_results[[paste0("c_", c)]]
  y1 <- results$y1
  est1 <- results$est1
  y2 <- results$y2
  est2 <- results$est2
  plot(est1, y1, xlab = "point estimate", ylab = "observed data", xlim = range(y1), pch = 4)
  points(est2, y2, col = 2, pch = 4)
  abline(0, 1, lty = 2)
  legend("bottomright", legend = c("SB1", "SB2"), col = c(1, 2), pch = c(4, 4))
  
  plot(est1, y1, xlab = "point estimate", ylab = "observed data", pch = 4, xlim = c(0, 10), ylim = c(0, 15))
  points(est2, y2, col = 2, pch = 4)
  abline(0, 1, lty = 2)
  abline(v = 5)
  legend("bottomright", legend = c("SB1", "SB2"), col = c(1, 2), pch = c(4, 4))
}
dev.off()
