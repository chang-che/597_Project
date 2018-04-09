ForwardBack2 <- function(obs.vec, pi.mar, tran.m, mean.vec, sd.vec){
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
    alphabar_sum <- sum(alpha_b_h)
    alpha_b_h <- alpha_b_h/alphabar_sum
    Ct.log  <-  Ct.log + log(alphabar_sum)
    ct_seq[i] <- 1/alphabar_sum 
    alpha.hat[i,] <- alpha_b_h
    alpha.log[i,] <- log(alpha_b_h) + Ct.log #
  }
  
  
  # backward
  beta_b_h <- as.double(rep(1, N))*ct_seq[Ti]
  beta.hat <- matrix(rep(NA, N*Ti), nrow = Ti)
  beta.hat[Ti,] <- beta_b_h
  
  for (i in (Ti-1):1){
    beta_b_h <- tran.m%*%(density.m[i+1,]*beta_b_h) 
    
    beta_b_h <- beta_b_h*ct_seq[i]
    beta.hat[i,] <- beta_b_h
  }
  return(list(alphalog = alpha.log, Loglikelihood = Ct.log, 
              alphahat = alpha.hat, betahat = beta.hat,
              ct = ct_seq, density_m = density.m))
}
