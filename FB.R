ForwardBack <- function(obs.vec, pi.mar, tran.m, mean.vec, sd.vec){
  N <- nrow(tran.m)
  Ti <- length(obs.vec)
  density.m <- matrix(rep(as.double(0), Ti*N), nrow = Ti, ncol = N)
  for (i in 1:Ti){
    density.m[i,] <- mapply(
      dnorm, 'x' = obs.vec[i], 'mean' = mean.vec, 'sd' = sd.vec, SIMPLIFY = T)
  }
  # Forward
  alpha_b_h <- as.double(pi.mar)
  alpha.log <- matrix(rep(NA, Ti*N), nrow = Ti)
  alpha.hat <- matrix(rep(0, Ti*N), nrow = Ti)
  Ct.log <- as.double(0) # it is actually the minus log of Ct
  ct_seq <- c(rep(NA, Ti))
  for (i in 1:Ti){
    if (i >1)
      alpha_b_h <- alpha_b_h%*%tran.m
    alpha_b_h <- alpha_b_h*density.m[i,] #now is bar, next becomes hat
    alpha_sum <- sum(alpha_b_h)
    alpha_b_h <- alpha_b_h/alpha_sum
    Ct.log  <-  Ct.log + log(alpha_sum)
    ct_seq[i] <- alpha_sum
    alpha.hat[i,] <- alpha_b_h
    alpha.log[i,] <- log(alpha_b_h) + Ct.log #
  }
  ct_seq <- 1/ct_seq
  ct_seq_log <- log(ct_seq)
  # backward
  beta_b_h <- as.double(rep(1, N))
  beta.log <- matrix(as.double(rep(0, N*Ti)), nrow = Ti)
  beta.hat <- matrix(as.double(rep(1/N, N*Ti)), nrow = Ti)#This is not correct
  Dt.log <- as.double(log(N))
  for (i in (Ti-1):1){
    beta_b_h <- tran.m%*%(density.m[i+1,]*beta_b_h) 
    beta.log[i,] <- log(beta_b_h) + Dt.log
    beta_sum <- sum(beta_b_h)
    beta_b_h <- beta_b_h/beta_sum
    beta.hat[i,] <- beta_b_h
    Dt.log <- Dt.log + log(beta_sum)
  }
  return(list(betalog = beta.log, alphalog = alpha.log, Loglikelihood = Ct.log, alpha_hat = alpha.hat, beta_hat = beta.hat))
}
