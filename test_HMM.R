density.m <- matrix(rep(NA, Ti*N), nrow = Ti)
mean.vec <- c(1,4,6,8)
sd.vec <- c(0.5, 1, 0.5, 4)

for (i in Ti){
  density.m[i,] <- mapply(dnorm, 'x' = obs.vec[i], 'mean' = mean.vec, 'sd' = sd.vec, SIMPLIFY = T)
}
density.m
