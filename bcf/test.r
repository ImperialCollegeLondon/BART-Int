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

n <- 250
sigma <- 0.1
homogeneous <- 1
linear <- 1

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
  y_treated <- bcf_fit$yhat
  y_controled <- bcf_fit$yhat
  y_treated[z == 0] <- (y_treated + tauhat)[z == 0]
  y_controled[z == 1] <- (y_controled - tauhat)[z == 1]
  sate <- mean(y_treated - y_controled)
  print(paste0("True SATE: ", 3, " . Estimated: ", sate, ". Abs. diff: ", abs(sate - 3)))
}

# Z = 1, 
# Y = Y(1)
# Y(0)
# E[Y(1) - Y(0)] = ATE
# Y(0) = Y(1) - ATE

# Y(1) = Y(0) + ATE


# Gaussian processes
computeGPBQ <- function(X, Y, dim, epochs, kernel="rbf", FUN, lengthscale=1, sequential=TRUE, measure) 
  #'Gaussian Process with Bayesian Quadrature
  #' 
  #'@description This function calculates the approxiamtion of integration using
  #'Gaussian Process, Bayesian Quadrature and Sequential Design
  #' 
  #'@param dim Integer; Dimension of input X 
  #'@param epochs Integer; Number of new data points
  #'@param FUN Function; The function to be integrated
  #'@param lengthscale Integer; The parameter in standard kernel
  #'
  #'@return List; A list containing meanValue (apprimation) and variance of GP method

{
  #define genz function
  genz <- FUN
  meanValueGP <- c()
  varianceGP <- c()
   
  N <- dim(X)[1]

  K <- matrix(0,nrow=N,ncol=N)
  jitter = 1e-7

  if (kernel == "rbf") {
    kernel <- rbfdot(.5/lengthscale^2)
  } else if (kernel == "matern32") {
     kernel <- maternKernelWrapper(lengthscale)
  }  
  
  # K = kernel(X)
  K = kernelMatrix(X, X, kernel = kernel)
  # compute the variance
  if (measure == "uniform"){
    int.points.1 <- X
    int.points.2 <- X
  } else {
    warning(paste("Only uniform (empirical) measure is allowed but got", measure))
  }
  cov <- kernel(int.points.1, int.points.2)
  var.firstterm <- mean(cov[upper.tri(cov)])
  cov <- kernel(int.points.1, X)
  z <- colMeans(cov) 
  covInverse <- chol2inv(chol(K + diag(jitter, nrow(K))))
  meanValueGP[1] <- t(z) %*% covInverse %*% Y
  tmp <- t(z)%*% covInverse %*% z 
  varianceGP[1] <- var.firstterm - tmp
  # cat(var.firstterm, tmp,"\n")

  # train
  if (epochs == 1){
    return (list("meanValueGP" = meanValueGP, "varianceGP" = varianceGP, "X" = X, "Y" = Y, "K" = K))
  }
}



whichKernel <- "rbf"
x_input_gp <- cbind(x_input, pihat, z)
setwd("../")
library(reticulate)
source("src/optimise_gp.R")
lengthscale <- optimise_gp_r(x_input_gp, matrix(y, nrow = n), kernel = whichKernel, epochs = 500)
print("...Finished training for the lengthscale")
source("src/GPBQ.R")
setwd("bcf/")
# Get GP posterior
predictionGPBQ <- computeGPBQ(
  x_input_gp,
  matrix(y, nrow = n),
  dim = dim(x_input_gp)[2],
  epochs = 1,
  kernel = whichKernel,
  FUN = mean,   # FUN is not used if epoch = 1
  lengthscale,
  sequential = FALSE,
  measure = "uniform"
)







