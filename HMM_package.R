Estep
function (x, Pi, delta, distn, pm, pn = NULL) 
{
  dfunc <- makedensity(distn)
  m <- nrow(Pi)
  n <- length(x)
  y <- forwardback(x, Pi, delta, distn, pm, pn)
  logbeta <- y$logbeta
  logalpha <- y$logalpha
  LL <- y$LL
  u <- exp(logalpha + logbeta - LL)
  v <- array(NA, dim = c(n - 1, m, m))
  for (k in 1:m) {
    logprob <- dfunc(x = x[-1], getj(pm, k), getj(pn, -1), 
                     log = TRUE)
    logPi <- matrix(log(Pi[, k]), byrow = TRUE, nrow = n - 
                      1, ncol = m)
    logPbeta <- matrix(logprob + logbeta[-1, k], byrow = FALSE, 
                       nrow = n - 1, ncol = m)
    v[, , k] <- logPi + logalpha[-n, ] + logPbeta - LL
  }
  v <- exp(v)
  return(list(u = u, v = v, LL = LL))
}
<bytecode: 0x114d91470>
  <environment: namespace:HiddenMarkov>
  
forwardback
function (x, Pi, delta, distn, pm, pn = NULL, fortran = TRUE) 
{
  m <- nrow(Pi)
  n <- length(x)
  dfunc <- makedensity(distn)
  prob <- matrix(as.double(0), nrow = n, ncol = m)
  for (k in 1:m) prob[, k] <- dfunc(x = x, getj(pm, k), pn, 
                                    log = FALSE)
  phi <- as.double(delta)
  logalpha <- matrix(as.double(rep(0, m * n)), nrow = n)
  lscale <- as.double(0)
  if (fortran != TRUE) {
    for (i in 1:n) {
      if (i > 1) 
        phi <- phi %*% Pi
      phi <- phi * prob[i, ]
      sumphi <- sum(phi)
      phi <- phi/sumphi
      lscale <- lscale + log(sumphi)
      logalpha[i, ] <- log(phi) + lscale
    }
    LL <- <-  <- lscale
  }
  else {
    if (!is.double(Pi)) 
      stop("Pi is not double precision")
    if (!is.double(prob)) 
      stop("prob is not double precision")
    memory0 <- rep(as.double(0), m)
    loop1 <- .Fortran("loop1", m, n, phi, prob, Pi, logalpha, 
                      lscale, memory0, PACKAGE = "HiddenMarkov")
    logalpha <- loop1[[6]]
    LL <- loop1[[7]]
  }
  logbeta <- matrix(as.double(rep(0, m * n)), nrow = n)
  phi <- as.double(rep(1/m, m))
  lscale <- as.double(log(m))
  if (fortran != TRUE) {
    for (i in seq(n - 1, 1, -1)) {
      phi <- Pi %*% (prob[i + 1, ] * phi)
      logbeta[i, ] <- log(phi) + lscale
      sumphi <- sum(phi)
      phi <- phi/sumphi
      lscale <- lscale + log(sumphi)
    }
  }
  else {
    memory0 <- rep(as.double(0), m)
    loop2 <- .Fortran("loop2", m, n, phi, prob, Pi, logbeta, 
                      lscale, memory0, PACKAGE = "HiddenMarkov")
    logbeta <- loop2[[6]]
  }
  return(list(logalpha = logalpha, logbeta = logbeta, LL = LL))
}
<bytecode: 0x110110158>
  <environment: namespace:HiddenMarkov>
  
Mstep.norm
function (x, cond, pm, pn) 
{
  nms <- sort(names(pm))
  n <- length(x)
  m <- ncol(cond$u)
  if (all(nms == c("mean", "sd"))) {
    mean <- as.numeric(matrix(x, nrow = 1) %*% cond$u)/apply(cond$u, 
                                                             MARGIN = 2, FUN = sum)
    sd <- sqrt(apply((matrix(x, nrow = n, ncol = m) - matrix(mean, 
                                                             nrow = n, ncol = m, byrow = TRUE))^2 * cond$u, MARGIN = 2, 
                     FUN = sum)/apply(cond$u, MARGIN = 2, FUN = sum))
    return(list(mean = mean, sd = sd))
  }
  if (all(nms == "mean")) {
    mean <- as.numeric(matrix(x, nrow = 1) %*% cond$u)/apply(cond$u, 
                                                             MARGIN = 2, FUN = sum)
    return(list(mean = mean))
  }
  if (all(nms == "sd")) {
    sd <- sqrt(apply((matrix(x, nrow = n, ncol = m) - matrix(pn$mean, 
                                                             nrow = n, ncol = m, byrow = FALSE))^2 * cond$u, MARGIN = 2, 
                     FUN = sum)/apply(cond$u, MARGIN = 2, FUN = sum))
    return(list(sd = sd))
  }
  stop("Invalid specification of parameters")
}
<bytecode: 0x10dd18678>
  <environment: namespace:HiddenMarkov>
  
Baum.Welch
function (x, Pi, delta, distn, pm, pn = NULL, nonstat = TRUE, 
          maxiter = 500, tol = 1e-05, prt = TRUE, posdiff = (distn[1] != 
                                                               "glm")) 
{
  .Deprecated("BaumWelch", package = "HiddenMarkov", msg = "'Baum.Welch' is deprecated.\n          Use 'BaumWelch' instead, see help('BaumWelch').")
  if (distn[1] != "glm") {
    Mstep <- parse(text = paste("Mstep.", distn, "(x, cond, pm, pn)", 
                                sep = ""))
  }
  else {
    Mstep <- parse(text = paste("Mstep.glm", "(x, cond, pm, pn, distn[2], distn[3])", 
                                sep = ""))
  }
  m <- nrow(Pi)
  n <- length(x)
  oldLL <- -Inf
  for (iter in 1:maxiter) {
    cond <- Estep(x, Pi, delta, distn, pm, pn)
    diff <- cond$LL - oldLL
    if (prt) {
      cat("iter =", iter, "\n")
      cat("LL =", formatC(cond$LL, digits = log10(1/tol) + 
                            2, format = "f"), "\n")
      cat("diff =", diff, "\n\n")
    }
    if (diff < 0 & posdiff) 
      stop("Worse log-likelihood")
    if (abs(diff) < tol) 
      break
    Pi <- diag(1/apply(cond$v, MARGIN = 2, FUN = sum)) %*% 
      apply(cond$v, MARGIN = c(2, 3), FUN = sum)
    if (nonstat) 
      delta <- cond$u[1, ]
    else delta <- compdelta(Pi)
    pm <- eval(Mstep)
    oldLL <- cond$LL
  }
  return(list(delta = delta, Pi = Pi, u = cond$u, v = cond$v, 
              pm = pm, LL = cond$LL, iter = iter, diff = diff))
}