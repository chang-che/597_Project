BaumWelch_hmm <- function(obs.vec, pi.mar, tran.m, mean.vec, sd.vec, maxiter, epsilon, silence = F){
  Ti <- length(obs.vec)
  N <- length(pi.mar)
  oldLL <- -Inf
  log.lik.iter <- c(-Inf)
  for (iter in 1:maxiter){
    FB <- ForwardBack2(obs.vec, pi.mar, tran.m, mean.vec, sd.vec)
    LL <- FB$Loglikelihood
    log.lik.iter <- c(log.lik.iter, LL)
    diff <- LL - oldLL
    if (silence == F){
      cat('iter=', iter, '\n')
      cat('LL=', round(LL, 10), '\n')
      cat('diff=', diff, '\n')
    }
    if (is.na(diff)){
      break
    }
    if(abs(diff)<epsilon){
      #stop('Algorithm converges!')
      break
    }
    #can add more criterions when needed
    gamma.log <- FB$alphalog+FB$betalog-LL
    gamma.m <- exp(gamma.log)
    AH <- FB$alphahat
    BH <- FB$betahat
    DM <- FB$density_m
    ct_vec <-FB$ct
    gamma.m <- AH*BH
    gamma.m <- t(sapply(1:Ti, function(x) gamma.m[x,]/ct_vec[x]))
    mean.vec <- apply(gamma.m, 2, function(x) sum(x*obs.vec)/sum(x))
    sd.vec <- sqrt(apply((matrix(obs.vec, nrow = Ti, ncol = N) - matrix(mean.vec,nrow = Ti, ncol = N, byrow = TRUE))^2 * gamma.m, MARGIN = 2,FUN = sum)/apply(gamma.m, MARGIN = 2, FUN = sum))
    pi.mar <- gamma.m[1,]
    tran.m <- (t(AH[1:Ti-1,])%*%(DM[2:Ti,]*BH[2:Ti,]))*tran.m #maybe a parenthesis is redundant.
    #tran.m <- t(scale(t(tran.m),center = F, scale = colSums(t(tran.m))))
    tran.m <- t(apply(tran.m, 1, function(x) x/sum(x)))
    oldLL <- LL
  }
  return(list(tran_m = tran.m, pi_mar = pi.mar, mean = mean.vec, sd = sd.vec, LL = LL, iter_LL = log.lik.iter))
}
