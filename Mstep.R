Mstep <- function(){
  tran.m <- (t(alpha.hat[1:Ti-1,])%*%(density.m[2:Ti,]*beta.hat[2:Ti,]))*tran.m
  
}