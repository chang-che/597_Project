library(invgamma)
# initials = number of random initializations for each possible N = 2,3,4,5
# maxiter = number of maximum iterations for BaumWelch algorithm to stop
# epsilon = threshold to determine whether diff of loglikelihood is small enough for the BM function
hmm_BIC<-function(obs.vec, initials, maxiter = 500, epsilon = 10^(-6)){ 
  n<-length(obs.vec)
  BIC<-matrix(NA,nrow=5,ncol=initials)
  sm <- mean(obs.vec)
  sv <- var(obs.vec)
  for(i in 2:5){
    k<- i*(i-1) + i*2 + i-1  #df of transition matrix + mean&sd + pi
    for(j in 1:initials){
      tran.m <- matrix(runif(i*i), nrow = i)
      tran.m <- t(apply(tran.m, 1, function(x) x/sum(x)))
      ##simulate the marginal distribution at time 1
      a <- runif(i,min = 0,max = 1)
      pi.mar <- a/sum(a) 
      ##simulate the mean and standard deviation
      mean.vec<-rep(NA,i)
      #sd.vec <- rep(NA,i)
      sd.vec <- rep(1,i)
      for(m in 1:i){
          #sd.vec[m]<-sqrt(rinvgamma(1,shape = (n-1)/2, rate = (n-1)*sv/2))
          mean.vec[m]<-rnorm(1, mean = sm, sd = sd.vec[m]/sqrt(n))
        }
      L<-BaumWelch_hmm(obs.vec, pi.mar, tran.m, mean.vec, sd.vec, maxiter, epsilon, silence = T)$LL
      BIC[i,j]<-log(n)*k-2*L
      }
    }
  #BIC1<-min(BIC[-1,], na.rm = T)
  BIC1<-apply(BIC[-1,], 1, function(x) min(x, na.rm = T))
  state<-which.min(BIC1) + 1
  return(list(BIC=BIC1,state=state))
}