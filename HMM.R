###############
#assumptions 1 the number of states is N (it will be tuned later)
# assumptions 2 the obs follows a mixture Guassian model
# number of obs = T

N <- 4
Ti <- 30
Obser <- rnorm(Ti)
# pi0 initial state
x <- runif(N)
pi_0 <- x/sum(x)
# transition matrix initial states
x1 = list()

for (i in 1:N) {
  x <- runif(N)
  x1[[i]] <- x/sum(x)
}
tran.m <- matrix(unlist(x1), ncol = N, byrow = T)
ksi <- c()
# initial of parameter of Guassian
listofnormals = list()
listofnormals[[1]] = list('mean' = 1, 'sd' = 2 )
listofnormals[[2]] = list('mean' = 3, 'sd' = 2 )
listofnormals[[3]] = list('mean' = 5, 'sd' = 3 )
listofnormals[[4]] = list('mean' = 10, 'sd' = 1 )
#forward algorithm
?dnorm
listofobs <- lapply(listofnormals, function(a){a[['x']] = Obser[1]
return(a)})
alpha1 <- pi_0*sapply(listofobs, function(a) do.call(dnorm, a))

list_alpha <- list()
list_alpha[[1]] <- alpha1

for(i in 2:Ti){
  listofobs <- lapply(listofnormals, function(a){a[['x']] = Obser[i]
  return(a)})
  list_alpha[[i]] <- sapply(1:N, function(j) sum(list_alpha[[i-1]]*tran.m[,j]*do.call(dnorm, listofobs[[j]])))
}
# when time is large, the probability will become zero because of truncation of R
list_alpha

#backward algorithm
betaT <- rep(1, 4)

list_beta <- list()
list_beta[[1]] <- betaT

for(i in 2:Ti){
  listofobs <- lapply(listofnormals, function(a){a[['x']] = Obser[i]
  return(a)})
  list_beta[[i]] <- sapply(1:N, function(j) sum(list_beta <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- [[i-1]]*tran.m[,j]*do.call(dnorm, listofobs[[j]])))
}





