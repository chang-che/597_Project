##We wish to calculate the probability of the observation sequence O,given the model lambda
##we need to consider every possible state sequence. 
##Let the model has N different kind of states.

ForwardBack2 <- function(obs.vec, pi.mar, tran.m, mean.vec, sd.vec){
  N <- nrow(tran.m)
  Ti <- length(obs.vec)
  density.m <- matrix(rep(as.double(0), Ti*N), nrow = Ti, ncol = N)##matrix of bj(Ot)
  for (i in 1:Ti){
     density.m[i,] <- mapply(
       dnorm, 'x' = obs.vec[i], 'mean' = mean.vec, 'sd' = sd.vec, SIMPLIFY = T)
   }
  # Forward
  ##First, we defined the forward variable alpha(j)=P(O1...Ot,qt=Si|lambda)
    ## 1) Initialization: alpha.1(i)=pi.i*bi(O1), 1<=i<=N
  alpha_b_h <- as.double(pi.mar)
  alpha.log <- matrix(rep(NA, Ti*N), nrow = Ti)
  alpha.hat <- matrix(rep(0, Ti*N), nrow = Ti)
  Ct.log <- as.double(0) # it is actually the minus log of Ct
  ct_seq <- c(rep(NA, Ti))
    ## 2) Induction: alpha.t+1=sum of state i from 1 to N of
    ## alpha.t(i)*aij*bj(Ot+1)
  for (i in 1:Ti){
    if (i >1)
      alpha_b_h <- alpha_b_h%*%tran.m #tran.m is the matrix of aij
    alpha_b_h <- alpha_b_h*density.m[i,] #now is bar, next becomes hat
    ## 3) Termination:P(O|lambda)=sum of alpha from state 1 to N,given final time T.
    alphabar_sum <- sum(alpha_b_h)
    alpha_b_h <- alpha_b_h/alphabar_sum
    Ct.log  <-  Ct.log + log(alphabar_sum)
    ct_seq[i] <- 1/alphabar_sum 
    alpha.hat[i,] <- alpha_b_h
    alpha.log[i,] <- log(alpha_b_h) + Ct.log # This output log of alpha is used for debug
  }
  
  
  # backward
  ##For backward variable, we defined beta.t(i)=P(Ot+1...OT|qt=Si,lambda)
    ## 1) Initialization: let beta.T(i)=1
  beta_b_h <- as.double(rep(1, N))*ct_seq[Ti]
  beta.hat <- matrix(rep(NA, N*Ti), nrow = Ti)
  beta.hat[Ti,] <- beta_b_h
  
    ## 2) Induction:beta.t(i)=sum of aij*bj(Ot+1)*beta.t+1(j) from state 1 to N, and j is the state.
    ## account for the observation sequence from time t+1 on
  for (i in (Ti-1):1){
    beta_b_h <- tran.m%*%(density.m[i+1,]*beta_b_h) 
    
    beta_b_h <- beta_b_h*ct_seq[i]
    beta.hat[i,] <- beta_b_h
  }
  return(list(alphalog = alpha.log, Loglikelihood = Ct.log, 
              alphahat = alpha.hat, betahat = beta.hat,
              ct = ct_seq, density_m = density.m))
}
