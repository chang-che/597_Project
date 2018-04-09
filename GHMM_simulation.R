####################simulation of Guassian HMM
##

set.seed(1234)
n_states = 4
n_time = 200
x1 = list()
for (i in 1:n_states) {
  x <- runif(n_states)
  x1[[i]] <- x/sum(x)
}

#tran.m.true <- matrix(unlist(x1), ncol = n_states, byrow = T)
#apply(tran.m, 1, sum)
tran.m.true <- matrix(c(1/2,1/2,0,0,
                        0,1/2,1/2,0,
                        0,0,1/2,1/2,
                        1/2,0,0,1/2), ncol = n_states, byrow = T)
#
listofnormals = list()
listofnormals[[1]] = list(1, 'mean' = 2, 'sd' = 1 )
listofnormals[[2]] = list(1, 'mean' = 4, 'sd' = 3 )
listofnormals[[3]] = list(1, 'mean' = 8, 'sd' = 3 )
listofnormals[[4]] = list(1, 'mean' = 14, 'sd' = 3 )
ini.vec = c(1, 0, 0, 0)
obs.vec = rep(NA, n_time)
obs.states = rep(NA, n_time)

for (i in 1:n_time) {
  x <- sample(1:n_states, 1, prob = ini.vec)
  obs.states[i] <- x
  obs.vec[i] = do.call(rnorm, listofnormals[[x]])
  ini.vec <- ini.vec %*% tran.m.true
}

obs.states
tran.m.true

