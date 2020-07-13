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

generate_x <- function(n, p) {
  x <- matrix(rnorm(n*p), ncol = p)
  x[, 4] <- x[, 4] > 1
  x[, 5] <- sample(c(1, 2, 3), n, replace = TRUE)
  return(data.frame(x[, 1:3], as.factor(x[, 4]), as.factor(x[, 5])))
}


args <- commandArgs(TRUE)
n <- as.double(args[1])
sigma <- as.double(args[2]) # observation noise of y
homogeneous <- as.double(args[3]) # 1 if true
linear <- as.character(args[4]) # 1 if linear

ntrain <- 250
ncandidate <- 50
n <- ntrain + ncandidate
sigma <- 0.1
homogeneous <- 1
linear <- 1

set.seed(1)
p <- 5

# Create covariates x
x <- generate_x(n, p)
x_input <- dbarts::makeModelMatrixFromDataFrame(x)

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

bcf_fit <- bcf(y[1:ntrain], z[1:ntrain], x_input[1:ntrain, ], x_input[1:ntrain, ], pihat[1:ntrain], nburn=2000, nsim=2000, nthin= 5)

# Get posterior samples of treatment effects
tau_post <- bcf_fit$tau
tauhat <- colMeans(tau_post)
plot(tau[1:ntrain], tauhat); abline(0,1)

# the posterior of the averaged treatment effects
print(paste0("Mean of |CATE - CATE_hat|: ", mean(abs(tau[1:ntrain] - tauhat))))

if (homogeneous){
  # Sample ATE
  y_treated <- bcf_fit$yhat
  y_controled <- bcf_fit$yhat
  y_treated[z == 0] <- (y_treated + tauhat)[z == 0]
  y_controled[z == 1] <- (y_controled - tauhat)[z == 1]
  sate <- mean(y_treated - y_controled)
  print(paste0("True SATE: ", 3, " . Estimated: ", sate, ". Abs. diff: ", abs(sate - 3)))
}



# BART
source("bcf/BART.R")
x_input_bart <- dbarts::makeModelMatrixFromDataFrame(data.frame(x[, 1:3], as.factor(x[, 4]), as.factor(x[, 5]), z))
BARTResults <- computeBART(x_input_bart, y, x_input_bart, y, 1)


# Gaussian processes
makeGPBQModelMatrix <- function(df, treatment) {
  return(dbarts::makeModelMatrixFromDataFrame(data.frame(df[, 1:3], as.factor(df[, 4]), as.factor(df[, 5]), treatment)))
}
# Uncomment the following to exclude treatment z in the data input 
# makeGPBQModelMatrix <- function(df) {
#   return(dbarts::makeModelMatrixFromDataFrame(df))
# }

# Needed for GPBQ to estimate \Pi[k_{*, X}] and \Pi\Pi[K] 
generate_x_GPBQ <- function(n) {
  x <- generate_x(n, 5)
  z <- rbinom(n, rep(1, n), propensity)
  return(makeGPBQModelMatrix(x, z))
}
# Uncomment the following to exclude treatment z in the data input 
# generate_x_GPBQ <- function(n) {
#   return(makeGPBQModelMatrix(generate_x(n, 5)))
# }

source("bcf/GPBQ.R")
# Uncomment the following to include treatment z in the data input 
trainX <- makeGPBQModelMatrix(x[1:ntrain,], z[1:ntrain])
trainY <- y[1:ntrain]
candidateX <- makeGPBQModelMatrix(x[-(1:ntrain),], z[-(1:ntrain)])
candidateY <- y[-(1:ntrain)]
# Uncomment the following to exclude treatment z in the data input 
# trainX <- makeGPBQModelMatrix(x[1:ntrain,])
# trainY <- y[1:ntrain]
# candidateX <- makeGPBQModelMatrix(x[-(1:ntrain),])
# candidateY <- y[-(1:ntrain)]

library(reticulate)
source("src/optimise_gp.R")
lengthscale <- optimise_gp_r(trainX, trainY, kernel = "rbf", epochs = 500)
print("...Finished training for the lengthscale")

GPBQResults <- computeGPBQWeighted(trainX, trainY, candidateX, candidateY, "rbf", lengthscale, 1)
GPBQResults$meanValueGP

