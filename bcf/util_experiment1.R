library(mvtnorm)

g <- function(w) {
  results <- rep(NA, length(w))
  results[w == 1] <- 2
  results[w == 2] <- -1
  results[w == 3] <- -4
  if (any(w != 1 | w != 2 | w != 3))
  return(results)
}

treatment_func <- function(x, homogeneous) {
  # Homogeneous
  if (homogeneous==1) {
    return(rep(3, dim(x)[1]))
  } else {
    # Heterogeneous
    return(1 + 2 * x[, 2] * as.integer(x[, 5]))
  }
}

prognostic_func <- function(x, linear) {
  if (linear == 1) {
    # g1
    results <- 1 + g(x[, 5]) + x[, 1] * x[, 3]
  } else {
    # g2
    results = -6 + g(x[, 5]) + 6 * abs(x[, 3] - 1)
  }
  return(results)
}

generate_x <- function(n, p) {
  x <- matrix(rnorm(n*p), ncol = p)
  x[, 4] <- x[, 4] > 1
  x[, 5] <- sample(c(1, 2, 3), n, replace = TRUE)
  return(data.frame(x[, 1:3], as.factor(x[, 4]), as.factor(x[, 5])))
}

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