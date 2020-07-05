library(bcf)
# Uncomment the following for the modified bcf package
# source("bcf/R/RcppExports.R")
# source("bcf/R/bcf.R")

# Experiment1
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
  if (homogeneous) {
    return(rep(3, dim(x)[1]))
  } else {
    # Heterogeneous
    return(1 + 2 * x[, 2] * x[, 5])
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


args <- commandArgs(TRUE)
n <- as.double(args[1])
sigma <- as.double(args[2]) # observation noise of y
homogeneous <- as.double(args[3]) # 1 if true
linear <- as.character(args[4]) # 1 if linear

# n <- 250
# homogeneous <- 0
# linear <- 1
# sigma <- 0.1

set.seed(1)
p <- 5

# Create covariates x
x <- matrix(rnorm(n*p), ncol = p)
x[, 4] <- x[, 4] > 1
x[, 5] <- sample(c(1, 2, 3), n, replace = TRUE)
x_input <- dbarts::makeModelMatrixFromDataFrame(data.frame(x[, 1:3], as.factor(x[, 4]), as.factor(x[, 5])))

# Compute tau, mu and pi
tau <- treatment_func(x, homogeneous)
mu <- prognostic_func(x, linear)
propensity <- 0.8 * pnorm(3 * mu / sd(mu) - 0.5 * x[, 1]) + 0.05 + 0.1 * runif(n)

# Model: 
# Y_i = f(x_i, Z_i) + epsilon_i
# epsilon_i ~ N(0, sigma^2)
# f(x_i, z_i) = mu(x_i) + tau(x_i) * z_i
z <- rbinom(n, rep(1, n), propensity)
y <- mu + tau * z + sigma * rnorm(n)

# If pi is unknown, we would need to estimate it here
pihat <- propensity

bcf_fit <- bcf(y, z, x_input, x_input, pihat, nburn=2000, nsim=2000, nthin = 5)

# Get posterior samples of treatment effects
tau_post <- bcf_fit$tau
tauhat <- colMeans(tau_post)
plot(tau, tauhat); abline(0,1)

# the posterior of the averaged treatment effects
print(paste0("Mean of |CATE - CATE_hat|: ", mean(abs(tau-tauhat))))

if (homogeneous){
  # Sample ATE
  y_treated <- y
  y_controled <- y
  y_treated[z == 0] <- (y + tauhat)[z == 0]
  y_controled[z == 1] <- (y - tauhat)[z == 1]
  sate <- mean(y_treated - y_controled)
  print(paste0("True SATE: ", 3, " . Estimated: ", sate, ". Abs. diff: ", abs(sate - 3)))
}


# Compute coverage

