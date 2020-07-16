library(mvtnorm)

prob_density <- function(x, linear=1) {
  n <- dim(x)[1]
  prob <- mvtnorm::dmvnorm(x[, c(1, 2, 3)])
  prob[x[, 4] == 0] <- prob[x[, 4] == 0] * pnorm(1)
  prob[x[, 4] == 1] <- prob[x[, 4] == 1] * (1 - pnorm(1))
  prob <- prob / 3
  mu <- prognostic_func(x[, 1:5], linear)
  propensity <- 0.8 * pnorm(3 * mu / sd(mu) - 0.5 * x[, 1]) + 0.05 + 0.1 * runif(n)
  prob <- prob * dbinom(x[, 6], size=rep(1, n), propensity)
  return (prob)
}