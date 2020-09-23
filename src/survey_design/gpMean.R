#Bayesian Quadrature with Gaussian Process
############
library(mvtnorm)
library(MASS)
library(kernlab)
library(MCMCglmm)
library(rdist)
maternKernelWrapper <- function(lengthscale = 1, sigma = 1) {
  maternKernel <- function(x, y) 
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
    if (any(is.na(Y))) {
      d <- pdist(X)
    } else {
      d <- cdist(X, Y)
    }
    k <- sigma^2 * (1 + sqrt(3) * d / lengthscale) * exp(- sqrt(3) * d / lengthscale)
    return (k)
  }
}

computeGPBQEmpirical <- function(nonOneHotX, X, Y, nonOneHotCandidateSet, candidateSet, candidateY, epochs, kernel="rbf", lengthscale, sequential=TRUE) 
{
  meanValueGP <- c()
  varianceGP <- c()
  condition <- c()

  emp_density_func <- build_density(rbind(nonOneHotX, nonOneHotCandidateSet))
  weights <- emp_density_func(nonOneHotCandidateSet)
   
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
  allSet <- as.matrix(rbind(X, candidateSet))
  X <- as.matrix(X)
  candidateSet <- as.matrix(candidateSet)
  
  K_all <- kernel(allSet)
  
  var.firstterm <- sum(K_all)/(nrow(allSet)^2)
  
  K = kernel(X)
  cov <- kernel(allSet, X)
  z <- colMeans(cov) 
  covInverse <- chol2inv(chol(K + diag(jitter, nrow(K))))
  meanValueGP[1] <- t(z) %*% covInverse %*% Y
  varianceGP[1] <- var.firstterm - t(z) %*% covInverse %*% z
  condition[1] <- rcond(K)
  # train
  if (epochs == 1){
    return (list("meanValueGP" = meanValueGP, "varianceGP" = varianceGP, "X" = X, "Y" = Y, "K" = K, "cond" = condition))
  }
  
  K_star_star <- kernel(candidateSet)
  K_star <- kernel(candidateSet, X)
  for (p in 2:epochs) {
    print(paste("GPBQ: Epoch =", p))
    K_prime <- diag(N+p-1)
    K_prime[1:(N+p-2), 1:(N+p-2)] <- K
    candidate_Var <- diag(K_star_star - K_star %*% solve(K + diag(jitter, nrow(K)), t(K_star)))
    candidate_Var <- candidate_Var * weights
    
    index <- which(candidate_Var == max(candidate_Var))[1]
    kernel_new_entry <- kernel(matrix(candidateSet[index,], nrow=1), X)
    K_prime[N+p-1,1:(N+p-2)] <- kernel_new_entry
    K_prime[1:(N+p-2),N+p-1] <- kernel_new_entry
    X <- rbind(X,candidateSet[index,])

    Y <- c(Y, candidateY[index])
    K <- K_prime
    
    # update kernel matrices
    K_star_star <- K_star_star[-index, -index]
    K_star_new <- matrix(nrow=(nrow(K_star)-1), ncol=(ncol(K_star)+1))
    K_star_new[1:(nrow(K_star)-1), 1:ncol(K_star)] <- K_star[-index,]
    K_star_new[1:(nrow(K_star)-1),(ncol(K_star)+1)] <- kernel(candidateSet[-index,], matrix(candidateSet[index,], nrow=1))
    K_star <- K_star_new
    
    # update candidateSet
    candidateSet <- candidateSet[-index,]
    weights <- weights[-index]
    candidateY <- candidateY[-index]
    
    # add in extra term obtained by integration
    # cov <- kernel(allSet, X)
    cov_new <- matrix(nrow=nrow(cov), ncol=nrow(X))
    cov_new[1:nrow(cov),1:(nrow(X)-1)] <- cov 
    cov_new[1:nrow(cov), nrow(X)] <- kernel(allSet, matrix(X[nrow(X),], ncol=ncol(X)))
    cov <- cov_new
    # z <- colMeans(K_star)
    z <- colMeans(cov)
    covInverse <- chol2inv(chol(K + diag(jitter, N+p-1)))
    meanValueGP[p] <- t(z) %*% covInverse %*% Y
    varianceGP[p] <- var.firstterm - t(z) %*% covInverse %*% z
    condition[p] <- rcond(K)
  }
  return (list("meanValueGP" = meanValueGP, "varianceGP" = varianceGP, "X" = X, "Y" = Y, "K" = K, "cond" = condition))
}
