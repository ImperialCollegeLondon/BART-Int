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
ntrain <- as.double(args[1])
ncandidate <- as.double(args[2])
n_seqential <- as.double(args[3])
sigma <- as.double(args[4]) # observation noise of y
homogeneous <- as.double(args[5]) # 1 if true
linear <- as.character(args[6]) # 1 if linear
num_cv <- as.character(args[7]) # 1 if linear

# ntrain <- 500
# ncandidate <- 2000
# sigma <- 0.1
# homogeneous <- 1
# linear <- 1
# n_seqential <- 20
n <- ntrain + ncandidate
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

if (homogeneous == 1){
  # Sample ATE
  y_treated <- bcf_fit$yhat
  y_controled <- bcf_fit$yhat
  y_treated[z == 0] <- (y_treated + tauhat)[z == 0]
  y_controled[z == 1] <- (y_controled - tauhat)[z == 1]
  ate <- mean(y_treated - y_controled)
  print(paste0("True ATE: ", 3, " . Estimated: ", ate, ". Abs. diff: ", abs(ate - 3)))
}

# Gaussian processes
makeGPBQModelMatrix <- function(df, treatment) {
  return(data.frame(df[, 1:3], as.factor(df[, 4]), as.factor(df[, 5]), treatment))
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
lengthscale <- optimise_gp_r(dbarts::makeModelMatrixFromDataFrame(trainX), trainY, kernel = "matern32", epochs = 500)
print("...Finished training for the lengthscale")

source("bcf/GPBQ.R")
GPBQResults <- computeGPBQEmpirical(
  X=trainX, 
  Y=trainY, 
  candidateX=candidateX, 
  candidateY=candidateY, 
  kernel="matern32", 
  lengthscale=lengthscale, 
  epochs=n_seqential
)
 
# BART
source("bcf/BART.R")
BARTResults <- computeBARTWeighted(trainX, trainY, candidateX, candidateY, n_seqential, num_cv=num_cv, linear=linear)
print(paste0("True ATE: ", 3, " . Estimated: ", BARTResults$meanValueBART, ". Abs. diff: ", abs(BARTResults$meanValueBART - 3)))

# Bayesian Quadrature methods: with BART, Monte Carlo Integration and Gaussian Process respectively
print("Final ATE:")
print(c("Actual integral:", 3))
print(c("BCF integral:", ate))
print(c("BART integral:", BARTResults$meanValueBART[(n_seqential+1)]))
print(c("GP integral:", GPBQResults$meanValueGP[(n_seqential+1)]))

results <- data.frame(
  "epochs" = c(1:(n_seqential+1)),
  "BARTMean" = BARTResults$meanValueBART, "BARTsd" = BARTResults$standardDeviationBART,
  "GPMean" = GPBQResults$meanValueGP, "GPsd" = sqrt(GPBQResults$varianceGP),
  "actual" = rep(3, n_seqential+1)
)

write.csv(results, paste("bcf/results/experiment_", ntrain, "_", num_cv,  ".csv"))