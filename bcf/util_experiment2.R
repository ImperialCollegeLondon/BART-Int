library(mvtnorm)

treatment_func <- function(x, mu) {
  alpha <- 1 + x[,4]^2*(2*x[,2]-1)
  alpha <- alpha - min(alpha)
  alpha <- 1*sd(mu)*alpha/sd(alpha)
  return(alpha)
}

prognostic_func <- function(x, linear) {
  lev = c(-1,2,0)
  
  # if (linear) {
  #   linear
  #   result = 1 + x[,1]*(2*x[,2] - 2*(1-x[,2])) + lev[x3]
  # } else {
  #   # nonlinear
  #   result = 10*pnorm(x[,1]/4+x[,5]/10) + lev[x[,3]]
  # }

  # nonlinear
  result = 10*pnorm(x[,1]/4+x[,5]/4) + lev[x[,3]]

  return(result)
}

generate_x <- function(n, p) {
  x1 = sample(c(-2,2),n,replace=TRUE) + 0.25*rnorm(n)
  x2 = rbinom(n,1,0.7)
  x3 = sample(1:3,n,replace=TRUE,prob = c(0.1,0.3,0.6))
  x3 = as.factor(x3)

  x4 = rnorm(n)
  x5 = rbinom(n,4,0.5)
  x = data.frame(x1,x2,x3,x4,x5)
  return(x)
}

prob_density <- function(x, linear=1) {
  n <- dim(x)[1]
  prob <- 0.5 * dnorm(x[, 1], 2, 0.25) + 0.5 * dnorm(x[, 1], -2, 0.25)
  prob[x[, 2] == 0] <- prob[x[, 2] == 0] * 0.3
  prob[x[, 2] == 1] <- prob[x[, 2] == 1] * 0.7
  prob[x[, 3] == 1] <- prob[x[, 3] == 1] * 0.1
  prob[x[, 3] == 2] <- prob[x[, 3] == 2] * 0.3
  prob[x[, 3] == 3] <- prob[x[, 3] == 3] * 0.6
  prob <- prob * dnorm(x[, 4])
  prob <- prob * dbinom(x[, 5], 4, 0.5)
  return(prob)
}