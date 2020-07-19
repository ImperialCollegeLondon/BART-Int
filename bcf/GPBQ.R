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
  matern <- function (trainX, trainY=NA) {
    if (any(is.na(trainY))) {
      d <- pdist(trainX)
    } else {
      d <- cdist(trainX, trainY)
    }
    k <- sigma^2 * (1 + sqrt(3) * d / lengthscale) * exp(- sqrt(3) * d / lengthscale)
    return (k)
  }
}

computeGPBQEmpirical <- function(X, Y, candidateX, candidateY, epochs, kernel="rbf", lengthscale, sequential=TRUE, linear=1) 
{
  meanValueGP <- c()
  varianceGP <- c()
  condition <- c()
  
  N <- dim(X)[1]
  
  K <- matrix(0,nrow=N,ncol=N)
  jitter = 1e-6
  
  if (kernel == "rbf") {
    kernel <- rbfdot(.5/lengthscale^2)
  }
  else if (kernel == "matern32") {
    kernel <- maternKernelWrapper_2(lengthscale)
  }  
  
  # compute the variance
  notOneHotX <- X
  X <- dbarts::makeModelMatrixFromDataFrame(X)
  X <- as.matrix(X)
  notOneHotCandidateX <- candidateX
  candidateX <- dbarts::makeModelMatrixFromDataFrame(candidateX)
  candidateX <- as.matrix(candidateX)
  
  allSet <- as.matrix(rbind(X, candidateX))
  xDim <- dim(X)[2] - 1
  allSetTreated <- allSet
  allSetTreated[, (xDim + 1)] <- 1
  allSetControl <- allSet
  allSetControl[, (xDim + 1)] <- 0
  
  K_allTreated <- kernel(allSetTreated)
  var.firsttermTreated <- sum(K_allTreated)/(nrow(allSetTreated)^2)
  
  K = kernel(X)
  cov <- kernel(allSetTreated, X)
  z <- colMeans(cov) 
  covInverse <- chol2inv(chol(K + diag(jitter, nrow(K))))
  predTreated <- t(z) %*% covInverse %*% Y
  
  K_allControl <- kernel(allSetControl)
  var.firsttermControl <- sum(K_allControl)/(nrow(allSetControl)^2)
  
  covControl <- kernel(allSetControl, X)
  zControl <- colMeans(covControl) 
  predControl <- t(zControl) %*% covInverse %*% Y
  treatment_effects <- predTreated - predControl
  meanValueGP[1] <- treatment_effects
  varianceGP[1] <- var.firsttermTreated + var.firsttermControl - t(z) %*% covInverse %*% z - t(zControl) %*% covInverse %*% zControl
  condition[1] <- rcond(K)
  # train
  if (epochs == 0){
    return (list("meanValueGP" = meanValueGP, "varianceGP" = varianceGP, "X" = X, "Y" = Y, "K" = K, "cond" = condition))
  }
  
  K_star_star <- kernel(candidateX)
  K_star <- kernel(candidateX, X)
  for (p in 1:epochs) {
    print(paste("GPBQ: Epoch =", p))
    K_prime <- diag(N+p)
    K_prime[1:(N+p-1), 1:(N+p-1)] <- K
    candidate_Var <- diag(K_star_star - K_star %*% solve(K + diag(jitter, nrow(K)), t(K_star))) * prob_density(notOneHotX, linear)
    
    index <- which(candidate_Var == max(candidate_Var))[1]
    kernel_new_entry <- kernel(matrix(candidateX[index,], nrow=1), X)
    K_prime[N+p,1:(N+p-1)] <- kernel_new_entry
    K_prime[1:(N+p-1),N+p] <- kernel_new_entry
    notOneHotX <- rbind(notOneHotX, notOneHotCandidateX[index,])
    X <- dbarts::makeModelMatrixFromDataFrame(notOneHotX)
    
    Y <- c(Y, candidateY[index])
    K <- K_prime
    
    # update kernel matrices
    K_star_star <- K_star_star[-index, -index]
    K_star_new <- matrix(nrow=(nrow(K_star)-1), ncol=(ncol(K_star)+1))
    K_star_new[1:(nrow(K_star)-1), 1:ncol(K_star)] <- K_star[-index,]
    K_star_new[1:(nrow(K_star)-1),(ncol(K_star)+1)] <- kernel(candidateX[-index,], matrix(candidateX[index,], nrow=1))
    K_star <- K_star_new
    
    # update candidateX
    notOneHotCandidateX <- notOneHotCandidateX[-index,]
    candidateX <- dbarts::makeModelMatrixFromDataFrame(notOneHotCandidateX)
    candidateY <- candidateY[-index]
    
    # add in extra term obtained by integration
    # cov <- kernel(allSet, X)
    cov_new <- matrix(nrow=nrow(cov), ncol=nrow(X))
    cov_new[1:nrow(cov),1:(nrow(X)-1)] <- cov 
    cov_new[1:nrow(cov), nrow(X)] <- kernel(allSetTreated, matrix(X[nrow(X),], ncol=ncol(X)))
    cov <- cov_new
    # z <- colMeans(K_star)
    z <- colMeans(cov)
    covInverse <- chol2inv(chol(K + diag(jitter, N+p)))
    predTreated <- t(z) %*% covInverse %*% Y
    
    covControl <- kernel(allSetControl, X)
    zControl <- colMeans(covControl) 
    predControl <- t(zControl) %*% covInverse %*% Y
    treatment_effects <- predTreated - predControl
    
    meanValueGP[p+1] <- treatment_effects
    varianceGP[p+1] <- var.firsttermTreated + var.firsttermControl - t(z) %*% covInverse %*% z - t(zControl) %*% covInverse %*% zControl
    condition[p+1] <- rcond(K)
  }
  return (list("meanValueGP" = meanValueGP, "varianceGP" = varianceGP, "X" = X, "Y" = Y, "K" = K, "cond" = condition))
}