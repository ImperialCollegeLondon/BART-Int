#Bayesian Quadrature with Gaussian Process
############
library(mvtnorm)
library(MASS)
library(kernlab)
library(rdist)

maternKernelWrapper <- function(lengthscale = 1, sigma = 1) {
  maternKernel <- function(x, y) 
    #'Matern 3/2 Kernel Function in GP
    #' 
    #'@description This function calculates the covariance k(x,y)
    #' 
    #'@param x
    #'@param y
    #'@param lengthscale Real; parameter in standard kernel
    #'@param sigma Real; parameter in standard kernel
    #'
    #'@return K; Covariance matrix between xprime and X
  {
    d <- sqrt(sum((x - y)^2))
    k <- sigma^2 * (1 + sqrt(3) * d / lengthscale) * exp(- sqrt(3) * d / lengthscale)
    return (k)
  }
  return (maternKernel)
}

rescale <- function(x) {x * attr(x, 'scaled:scale') + attr(x, 'scaled:center')}
maternKernelWrapper_2 <- function(lengthscale=1, sigma=1) {
  matern <- function (X, Y=NA) {
    if (is.na(Y)) {
      d <- pdist(X)
    } else {
      d <- cdist(X, Y)
    }
    k <- sigma^2 * (1 + sqrt(3) * d / lengthscale) * exp(- sqrt(3) * d / lengthscale)
    return (k)
  }
}


computeGPBQ_matern <- function(X, Y, dim, epochs, kernel="rbf", FUN, lengthscale=1, sequential=TRUE, measure) 
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

  kernel <- maternKernelWrapper_2(lengthscale)
  
  K = kernel(X)
  # compute the variance
  if (measure == "uniform"){
    int.points.1 <- replicate(dim, runif(10000))
    int.points.2 <- replicate(dim, runif(10000))
  } else if (measure == "gaussian") {
    int.points.1 <- replicate(dim, rtnorm(10000, mean = 0.5, lower=0, upper=1))
    int.points.2 <- replicate(dim, rtnorm(10000, mean = 0.5, lower=0, upper=1))
  } else if (measure == "exponential") {
    int.points.1 <- replicate(dim, rexp(10000))
    int.points.2 <- replicate(dim, rexp(10000))
  }
  cov <- kernel(int.points.1, int.points.2)
  var.firstterm <- mean(cov[upper.tri(cov)])
  if (measure == "uniform"){
    int.points.1 <- replicate(dim, runif(1000000))
  } else if (measure == "gaussian") {
    int.points.1 <- replicate(dim, rtnorm(1000000, mean = 0.5, lower=0, upper=1))
  } else if (measure == "exponential") {
    int.points.1 <- replicate(dim, rexp(10000))
  }
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
  for (p in 2:epochs) {
   
    print(paste("GPBQ: Epoch =", p))
    candidateSetNum <- 100
    if (measure == "uniform") {
      candidateSet <- replicate(dim, runif(candidateSetNum))
      weights <- 1
    } else if (measure == "gaussian") {
      candidateSet <- replicate(dim, rtnorm(candidateSetNum, mean = 0.5, lower = 0, upper = 1))
      weights <- dtnorm(candidateSet, mean=0.5, lower = 0, upper = 1)
    } else if (measure == "exponential") {
      candidateSet <- replicate(dim, rexp(candidateSetNum))
      weights <- dexp(candidateSet)      
    }
    K_prime <- diag(N+p-1)
    K_prime[1:(N+p-2), 1:(N+p-2)] <- K
    
    if (sequential){
      K_star_star <- kernel(candidateSet)
      K_star <- kernel(candidateSet, X)
      candidate_Var <- diag(K_star_star - K_star %*% chol2inv(chol(K + diag(jitter, N+p-2))) %*% t(K_star))
      candidate_Var <- candidate_Var * weights
      index <- which(candidate_Var == max(candidate_Var))[1]
    }
    else {    
      index <- sample(1:candidateSetNum, 1)
    }
    kernel_new_entry <- kernel(matrix(candidateSet[index,], nrow=1), X)   
    
    K_prime[N+p-1,1:(N+p-2)] <- kernel_new_entry
    K_prime[1:(N+p-2),N+p-1] <- kernel_new_entry
    X <- rbind(X,candidateSet[index,])
    additionalResponse <- as.matrix( t(candidateSet[index,]), ncol = length(candidateSet[index,]) )

    Y <- c(Y, genz(additionalResponse))
    K <- K_prime
    
    # add in extra term obtained by integration
    cov <- kernel(int.points.1, X)
    z <- colMeans(cov)
    covInverse <- chol2inv(chol(K + diag(jitter, N+p-1)))
    meanValueGP[p] <- t(z) %*% covInverse %*% Y
    varianceGP[p] <- var.firstterm - t(z) %*% covInverse %*% z
  }

  return (list("meanValueGP" = meanValueGP, "varianceGP" = varianceGP, "X" = X, "Y" = Y, "K" = K))
}


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
  }
  else if (kernel == "matern32") {
     kernel <- maternKernelWrapper(lengthscale)
  }  
  
  K = kernelMatrix(kernel, X)
  # compute the variance
  if (measure == "uniform"){
    int.points.1 <- replicate(dim, runif(10000))
    int.points.2 <- replicate(dim, runif(10000))
  } else if (measure == "gaussian") {
    int.points.1 <- replicate(dim, rtnorm(10000, mean = 0.5, lower=0, upper=1))
    int.points.2 <- replicate(dim, rtnorm(10000, mean = 0.5, lower=0, upper=1))
  }
  cov <- kernelMatrix(kernel, int.points.1, int.points.2)
  var.firstterm <- mean(cov[upper.tri(cov)])
  if (measure == "uniform"){
    int.points.1 <- replicate(dim, runif(1000000))
  } else if (measure == "gaussian") {
    int.points.1 <- replicate(dim, rtnorm(1000000, mean = 0.5, lower=0, upper=1))
  }
  cov <- kernelMatrix(kernel, int.points.1, X)
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
  for (p in 2:epochs) {
   
    print(paste("GPBQ: Epoch =", p))
    candidateSetNum <- 100
    if (measure == "uniform") {
      candidateSet <- replicate(dim, runif(candidateSetNum))
      weights <- 1
    } else if (measure == "gaussian") {
      candidateSet <- replicate(dim, rtnorm(candidateSetNum, mean = 0.5, lower = 0, upper = 1))
      weights <- dtnorm(candidateSet, mean=0.5, lower = 0, upper = 1)
    }
	  candidateSet <- replicate(dim, runif(candidateSetNum))
    K_prime <- diag(N+p-1)
    K_prime[1:(N+p-2), 1:(N+p-2)] <- K
    
    if (sequential){
      K_star_star <- kernelMatrix(kernel, candidateSet)
      K_star <- kernelMatrix(kernel, candidateSet, X)
      candidate_Var <- diag(K_star_star - K_star %*% chol2inv(chol(K + diag(jitter, N+p-2))) %*% t(K_star))
      candidate_Var <- candidate_Var * weights
      index <- which(candidate_Var == max(candidate_Var))[1]
    }
    else {    
      index <- sample(1:candidateSetNum, 1)
    }
    
    kernel_new_entry <- kernelMatrix(kernel, matrix(candidateSet[index,], nrow=1), X)
    
    K_prime[N+p-1,1:(N+p-2)] <- kernel_new_entry
    K_prime[1:(N+p-2),N+p-1] <- kernel_new_entry
    X <- rbind(X,candidateSet[index,])
    additionalResponse <- as.matrix( t(candidateSet[index,]), ncol = length(candidateSet[index,]) )

    Y <- c(Y, genz(additionalResponse))
    K <- K_prime
    
    # add in extra term obtained by integration
    cov <- kernelMatrix(kernel, int.points.1, X)
    z <- colMeans(cov)
    covInverse <- chol2inv(chol(K + diag(jitter, N+p-1)))
    meanValueGP[p] <- t(z) %*% covInverse %*% Y
    varianceGP[p] <- var.firstterm - t(z) %*% covInverse %*% z
  }

  return (list("meanValueGP" = meanValueGP, "varianceGP" = varianceGP, "X" = X, "Y" = Y, "K" = K))
}
