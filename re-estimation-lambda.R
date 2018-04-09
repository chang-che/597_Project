###############
#assumptions 1 the number of states is N (it will be tuned later)
# assumptions 2 the obs follows a mixture Guassian model
# number of obs = T

N <- 4
Ti <- length(obs.vec) #use the simulated observations as the one
Obser <- obs.vec
log.lik.iter <- c(-1, 0)

###################################################
#Initialization
###################################################

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

# intial setting of Guassian
listofnormals = list()
listofnormals[[1]] = list('mean' = 1, 'sd' = 2 )
listofnormals[[2]] = list('mean' = 3, 'sd' = 2 )
listofnormals[[3]] = list('mean' = 6, 'sd' = 3 )
listofnormals[[4]] = list('mean' = 10, 'sd' = 3 )

norm.m = matrix(unlist(listofnormals), byrow = T, ncol = 2)
colnames(norm.m) <- c('mean', 'sd')
mu.vec <- norm.m[,'mean']
sd.vec <- norm.m[,'sd']

epsilon <- 1e-6
iter <- 1

mu.list <- list()
mu.list[[1]] <- mu.vec
sd.list <- list()
sd.list[[1]] <- sd.vec

###################################

while(abs(log.lik.iter[iter+1] - log.lik.iter[iter]) > epsilon){
  

###################################################################################
#Both algorithms are based on recursive programming i.e. Time order
#forward algorithm


alpha1 <- pi_0*mapply(dnorm, 'x' = Obser[1], 'mean' = mu.vec, 'sd' = sd.vec)
list_alpha <- list()
list_alpha[[1]] <- alpha1

for(t in 2:Ti){
  
  list_alpha[[t]] <- sapply(1:N, function(j) 
    sum(list_alpha[[t-1]]*tran.m[,j]*mapply(dnorm, 'x' = Obser[t], 'mean' = mu.vec[j], 'sd' = sd.vec[j])))
}
# when time is large, the probability will become zero because of truncation of R
list_alpha.m <- matrix(unlist(list_alpha), byrow = F, ncol = Ti)

###################################################################################
# the log likelihood under current paramters
ct.vec <- apply(list_alpha.m, 2, function(x) 1/sum(x))
log.lik.iter <- c(log.lik.iter, -sum(log(ct.vec)))

#backward algorithm
betaT <- rep(1, N)

list_beta <- list()

list_beta[[Ti]] <- betaT


for(t in rev(1:(Ti-1))){
  list_beta[[t]] <- sapply(1:N, function(i) 
    sum(list_beta[[t+1]]*tran.m[i,]*mapply(dnorm, 'x'= Obser[t+1], 'mean' = mu.vec, 'sd' = sd.vec)))
}
list_beta.m <- matrix(unlist(list_beta), ncol = Ti, byrow = F )

# scaling the matrix of forward and backward
alpha.m = scale(list_alpha.m, center = F, scale = colSums(list_alpha.m))
beta.m = scale(list_beta.m, center = F, scale = colSums(list_beta.m))

# density values matrix of observation
density.l = list()
i = 1
for (x in Obser) {
  density.l[[i]] <- mapply(dnorm, 'x' = x, 'mean' = mu.vec, 'sd' = sd.vec)
  i = i+1
}

density.m <- matrix(unlist(density.l), byrow = F, ncol = Ti)

################################################
# reestimation of transition matrix
tran.new = matrix(rep(NA, N*N), nrow = N)
#tran.l <- list()
for(i in 1:N){
  for(j in 1:N){
    tran.new[i,j] <- sum(alpha.m[i,1:Ti-1]*tran.m[i,j]*density.m[j, 2:Ti]*beta.m[j, 2:Ti]) 
  }
  
}

#x = t(scale(t(tran.new), center = F, scale = colSums(t(tran.new))))
# apply(tran.new, 1, function(x) scale(x, center = F, scale = sum(x)))

y = t(apply(tran.new, 1, function(x) x/sum(x)))
#y
#apply(y, 1, sum)
#y == x

tran.new <- y
tran.m <- tran.new
######################################
#reestimation of parameters of normal distribution

gamma.m <- list_alpha.m*list_beta.m
gamma.m <- scale(gamma.m, center = F, scale = colSums(gamma.m))

sd.vec <- c()
for(i in 1:N){
  sd.vec[i] <- sqrt(sum(gamma.m[i,]*(Obser-mu.vec[i])^2)/sum(gamma.m[i,]))
  sd.vec[i] <- sqrt(sum(gamma.m[i,]*(Obser-mu.vec[i])^2)/sum(gamma.m[i,]))
}
sd.list[[iter]] <- sd.vec# trace of sd


mu.vec <- apply(gamma.m, 1, function(x) sum(x*Obser)/sum(x))

mu.list[[iter]] <- mu.vec # trace of mu
# sd.vec <- apply(gamma.m, 1, function(x) sum(x*(Obser-mu.vec)^2)/sum(x))


pi_0 <- gamma.m[,1]
iter <- iter+1
}

plot(log.lik.iter)

mu.list
sd.list
mu.vec
