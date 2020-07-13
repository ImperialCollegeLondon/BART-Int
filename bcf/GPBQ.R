library(mvtnorm)
library(MASS)
library(kernlab)
library(rdist)
# source("bcf/BART.R")   # importing density function for x


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
  matern <- function (trainX, trainY=NA) {
    if (is.na(trainY)) {
      d <- pdist(trainX)
    } else {
      d <- cdist(trainX, trainY)
    }
    k <- sigma^2 * (1 + sqrt(3) * d / lengthscale) * exp(- sqrt(3) * d / lengthscale)
    return (k)
  }
}


computeGPBQWeighted <- function(trainX, trainY, candidateX, candidateY, kernel="rbf", lengthscale, num_iterations, save_k=FALSE, save_posterior_dir="bcf/synthetic", num_cv="default") 
{
  N <- dim(trainX)[1]
  xdim <- dim(trainX)[2]
  
  print(c("Adding number of new training data:", num_iterations))
  # outputs
  meanValue <- rep(0, num_iterations)
  standardDeviation <- rep(0, num_iterations)
  trainData <- cbind(trainX, trainY)
  fullData <- rbind(trainX, candidateX)

  meanValueGP <- c()
  varianceGP <- c()
   
  K <- matrix(0, nrow=N, ncol=N)
  jitter = 1e-7

  if (kernel == "rbf") {
    kernel <- rbfdot(.5/lengthscale^2)
  }
  else if (kernel == "matern32") {
     kernel <- maternKernelWrapper_2(lengthscale)
  }  
  
  int.points.1 <- generate_x_GPBQ(5000)
  int.points.2 <- generate_x_GPBQ(5000)
  K <- kernelMatrix(kernel, trainX)
  # compute the variance weighted by probability density
  cov <- kernelMatrix(kernel, int.points.1, int.points.2)
  var.firstterm <- mean(cov[upper.tri(cov)])
  cov <- kernelMatrix(kernel, int.points.1, trainX)
  z <- colMeans(cov) 
  covInverse <- chol2inv(chol(K + diag(jitter, nrow(K))))
  meanValueGP[1] <- t(z) %*% covInverse %*% trainY
  tmp <- t(z)%*% covInverse %*% z 
  varianceGP[1] <- var.firstterm - tmp
  # cat(var.firstterm, tmp,"\n")

  # train
  if (num_iterations == 1){
    return (list("meanValueGP" = meanValueGP, "varianceGP" = varianceGP, "X" = trainX, "Y" = trainY, "K" = K))
  }
  for (i in 2:num_iterations) {
   
    print(paste("GPBQ: Epoch =", i))
    K_prime <- diag(N+i-1)
    K_prime[1:(N+i-2), 1:(N+i-2)] <- K
    
    K_star_star <- kernelMatrix(kernel, candidateX)
    K_star <- kernelMatrix(kernel, candidateX, trainX)
    candidate_var <- diag(K_star_star - K_star %*% chol2inv(chol(K + diag(jitter, N+i-2))) %*% t(K_star)) * prob_density(candidateX)
    index <- sample(which(candidate_var == max(candidate_var)), 1)
    
    # Update kernel matrix
    kernel_new_entry <- kernelMatrix(kernel, matrix(candidateX[index,], nrow=1), trainX)
    K_prime[N+i-1, 1:(N+i-2)] <- kernel_new_entry
    K_prime[1:(N+i-2), N+i-1] <- kernel_new_entry
    K <- K_prime
  
    # Update train and candidate X
    trainX <- rbind(trainX, candidateX[index,])
    candidateX <- candidateX[-index,]
    # Update train and candidate Y (for output)
    trainY <- c(trainY, candidateY[index])
    candidateY <- candidateY[-index]
    
    # add in extra term obtained by integration
    cov <- kernelMatrix(kernel, int.points.1, trainX)
    z <- colMeans(cov)
    covInverse <- chol2inv(chol(K + diag(jitter, N+i-1)))
    meanValueGP[i] <- t(z) %*% covInverse %*% trainY
    varianceGP[i] <- var.firstterm - t(z) %*% covInverse %*% z
  }
  
  trainData <- cbind(trainX, trainY)

  if (save_k == TRUE) {
    save(list("K" = K), file = paste0(save_posterior_dir, "/GPBQ_k_synthetic_%s_%s" %--% c(i, num_cv), ".RData"))
  }
  return (list("meanValueGP" = meanValueGP, "varianceGP" = varianceGP, "trainData" = trainData, "K" = K))
}
