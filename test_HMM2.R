N <- 4
Ti <- length(obs.vec)
maxiter <- 500
epsilon <- 1e-6
pi.mar <- as.double(rep(1/N,N))
oldLL <- -Inf
log.lik.iter <- c(-Inf)
mean.vec <- c(2,3,9,10)
sd.vec <- c(1.5,2,3,4)
tran.m <- matrix(as.double(rep(1/N, N*N)), ncol = N, byrow = T)
FB <- ForwordBack(obs.vec, pi.mar, tran.m, mean.vec, sd.vec)
LL <- FB$Loglikelihood
log.lik.iter <- c(log.lik.iter, LL)
diff <- LL - oldLL
cat('iter=', iter, '\n')
cat('LL=', round(LL, 10), '\n')
cat('diff=', diff, '\n')
if(abs(diff)<epsilon){
  stop('Algorithm converges!')
}
gamma.log <- FB$alphalog+FB$betalog-LL
gamma.m <- exp(gamma.log)
alpha.hat <- FB$alpha_hat
beta.hat <- FB$beta_hat
mean.vec <- apply(gamma.m, 2, function(x) sum(x*obs.vec)/sum(x))
sd.vec <- sqrt(apply((matrix(obs.vec, nrow = Ti, ncol = N) - matrix(mean.vec,nrow = Ti, ncol = N, byrow = TRUE))^2 * gamma.m, MARGIN = 2,FUN = sum)/apply(gamma.m, MARGIN = 2, FUN = sum))
pi.mar <- gamma.m[1,]
tran.m <- (t(alpha.hat[1:Ti-1,])%*%(density.m[2:Ti,]*beta.hat[2:Ti,]))*tran.m #maybe a parenthesis is redundant.
tran.m <- t(scale(t(tran.m),center = F, scale = colSums(t(tran.m))))
oldLL <- LL

head(FB$alphalog)
head(FB$betalog)
FB$Loglikelihood