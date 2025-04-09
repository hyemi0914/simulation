library(GIGrvg)
###--------------------------------------------------------------###
###        R code for implementing Scaled beta (SB) prior        ###
###            and inverse rescaled beta (IRB) prior             ###
###--------------------------------------------------------------###

##     Function of Miller's sampling algorithm
## (used in MCMC algorithm in the proposed priors)
miller <- function(x, mu, a0, b0, ep, M){
  n <- length(x)
  R <- sum(log(x))
  S <- sum(x)
  T. <- S / mu - R + n * log(mu) - n
  A <- a0 + n / 2; B <- b0 + T.
  for(j in 1:M){
    a <- A / B
    A <- a0 - n * a + n * a^2 * trigamma(a)
    B <- b0 + (A - a0) / a - n * log(a) + n * digamma(a) + T.
    if(abs(a / (A / B) - 1) < ep){
      return(c(A = A, B = B, j = j))
    }
  }
  return(c(A = A, B = B, j = M + 1))
}



###    Scaled Beta (SB) prior    ###
## model settings 
# y \sim Ga(A_y, B_y/lambda)
# lambda \sim Ga(1+tau*u, beta*tau*u)
# u \sim SB(a, b)
# tau \sim Ga(a_ta, b_ta)
# beta \sim Ga(a_be, b_be)
## IMPUT
# y: sequence of observed values
# A_y, B_y: fixed values in the sampling model 
# (a_ta, b_ta): hyperparameter in prior of tau
# (a_be, b_be): hyperparameter in prior of beta
# a, b: shape parameters of SB prior 
# ep, M: tuning parameters in the Miller's algorithm
# mc: MCMC length 
# burn: burn-in period
## OUTPUT
# posterior samples of parameters 

SB_prior <- function(y1, A_y, B_y, a_ta, b_ta, a_be, b_be, a, b, de=10^(-8), ep=10^(-8), M=10, mc, burn){
  # preparation
  n <- length(y1)
  # initial values
  la <- A_y / (B_y * y1)
  rho. <- rep(1, n)
  nu <- rep(1, n)
  ta <- 1
  be <- 1
  
  # objects to store posterior samples
  La <- matrix(NA, mc, n)
  Rho. <- matrix(NA, mc, n)
  Nu <- matrix(NA, mc, n)
  Ta <- rep(NA, mc)
  Be <- rep(NA, mc)
  
  # MCMC
  for(iota in 1:mc){
    # la
    la <- rgamma(n, shape = A_y + 1 + nu, rate = B_y * y1 + be * nu)
    La[iota, ] <- la
    
    # rho
    rho._proposal <- sapply(2 * (1 + nu / ta) + 2 * de, function(twicerate) rgig(1, lambda = a * D + b, chi = 2 * de, psi = twicerate))
    uniform <- runif(n, min = 0, max = 1)
    logratio <- dgamma(rho._proposal, shape = a * D + b, rate = 1 + nu / ta, log = TRUE) - dgamma(rho., shape = a * D + b, rate = 1 + nu / ta, log = TRUE) - ((a + b - 1) * log(rho._proposal) - ((2 * (1 + nu / ta) + 2 * de) * rho._proposal + (2 * de) / rho._proposal) / 2 - (a + b - 1) * log(rho.) + ((2 * (1 + nu / ta) + 2 * de) * rho. + (2 * de) / rho.) / 2)
    rho. <- ifelse(test = (log(uniform) <= logratio), yes = rho._proposal, no = rho.)
    Rho.[iota, ] <- rho.
    
    # be
    be <- rgamma(1, shape = a_be + sum(nu + 1), rate = b_be + sum(nu * la))
    Be[iota] <- be
    
    # ta
    ta <- rgig(1, lambda = a_ta - n * a, chi = 2 * sum(rho. * nu), psi = 2 * b_ta)
    Ta[iota] <- ta
    
    # nu
    ABJ <- apply(rbind(la, a, rho. / ta), 2, function(s) miller(x = s[1], mu = 1 / be, a0 = s[2], b0 = s[3], ep = ep, M = M))
    A <- ABJ[1, ]
    B <- ABJ[2, ]
    nu_proposal <- rgamma(n, shape = A, rate = B)
    uniform <- runif(n, min = 0, max = 1)
    logratio <- ((a - 1) * log(nu_proposal) - rho. * nu_proposal / ta + nu_proposal * log(be * nu_proposal) - lgamma(nu_proposal) + nu_proposal * log(la) - la * be * nu_proposal - (A - 1) * log(nu_proposal) + B * nu_proposal) - ((a - 1) * log(nu) - rho. * sum(nu / ta) + nu * log(be * nu) - lgamma(nu) + nu * log(la) - la * be * nu - (A - 1) * log(nu) + B * nu)
    nu <- ifelse(test = (log(uniform) <= logratio), yes = nu_proposal, no = nu)
    Nu[iota, ] <- nu
  }
  
  # summary
  Res <- list(La = 1/La[-(1:burn), , drop = FALSE], Be = Be[-(1:burn)], Ta = Ta[-(1:burn)], 
              Rho. = Rho.[-(1:burn), , drop = FALSE], Nu = Nu[-(1:burn), , drop = FALSE])
  return(Res)
}
