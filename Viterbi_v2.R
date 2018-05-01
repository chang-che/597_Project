#Problem2:Viterbi algorithm is to find the maximum state sequence probability for the observation
# By induction, delta.t+1(j)=max of delta.t(i)*aij*bj(Ot+1), among state i from 1 to N
#delta.t(i)is the highest probability at time t, which ends in state Si
Vit_hmm <- function(obs.vec, pi.mar, tran.m, mean.vec, sd.vec){
  N <- nrow(tran.m)
  Ti <- length(obs.vec)
  density.m <- matrix(rep(as.double(0), Ti*N), nrow = Ti, ncol = N)
  for (i in 1:Ti){
    density.m[i,] <- mapply(
      dnorm, 'x' = obs.vec[i], 'mean' = mean.vec, 'sd' = sd.vec, 'log' = T, SIMPLIFY = T)
  }
  tran.m.log <- log(tran.m)
  delta_log <- matrix(rep(NA, N*Ti), ncol = N)
  ## 1)Initialization: delta.1(i)=(pi.i)*bi(O1),1<=i<=N
  ##                   Psi.1(i)=0
  delta_log[1,] <- log(pi.mar) + density.m[1,]
  
  #keep track of states
  bt <- rep(NA, Ti)
  delta_ind <- matrix(rep(NA, (Ti-1)*N), ncol = N)
  ## 2) Recursion??delta.t(j)= max of delta.t-1(i)*aij*bj(Ot) among state i from 1 to N
  for (t in 2:Ti){
    delta_max <- sapply(1:N, function(j) tran.m.log[,j]+delta_log[t-1, ]+density.m[t,j])
    delta_ind[t-1, ] <- apply(delta_max, 2, which.max)#find i, where helps delta*aij reaches the max at time t-1. 2 presents for column
    delta_log[t,] <- apply(delta_max,2,max)#maximum value of delta at time t-1
  }
  ## 3) Termination: P=max of delta.T(i), among state i from 1 to N
  ##                 qT=i which let P reaches the max
  bt[Ti] <- which.max(delta_log[Ti, ])
  for (i in (Ti-1):1) {
    bt[i] <- delta_ind[i,bt[i+1]]
  }
  return(bt)
}