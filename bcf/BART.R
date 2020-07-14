library(lhs)
library(dbarts)
library(data.tree)
library(BART)
library(matrixStats)
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

computeBARTWeighted <- function(trainX, trainY, candidateX, candidateY, num_iterations, save_posterior = FALSE, save_posterior_dir = "bcf/results", num_cv = "default", linear=1) {
  dim <- ncol(trainX)
  print(c("Adding number of new training data:", num_iterations))
  # outputs
  meanValue <- rep(0, num_iterations)
  standardDeviation <- rep(0, num_iterations)
  trainData <- cbind(trainX, trainY)
  fullData <- rbind(trainX, candidateX)
  oppositeFullData <- fullData
  oppositeFullData[, dim] <- (fullData[, dim] + 1) %% 2
  
  colnames(trainData)[dim + 1] <- "response"

  # generate extra training data using the scheme (see pdf)
  if (num_iterations == 0) {
    
    # set seed to enable reproduction of the results
    print(c("BART: Epoch=", 1))
    # first build BART model
    sink("/dev/null")
    model <- bart(trainData[, 1:dim], trainData[, dim + 1], keeptrees = TRUE, keepevery = 3L, nskip = 1000, ndpost = 5000, ntree = 50, k = 2, usequant = FALSE)
    # model <- wbart(bartModelMatrix(trainData[,1:dim]), trainData[,(dim+1)], bartModelMatrix(candidateX), keepevery=3L, nskip=500L, ndpost=5000L, ntree=50)
    sink()
    
    pred <- predict(model, fullData)
    if (save_posterior == TRUE) {
      posterior_samples <- list("posterior_samples" = rowMeans(pred))
      save(posterior_samples, file = paste(save_posterior_dir, "/posterior_BART_synthetic_%s_%s" %--% c(i, num_cv), ".RData", sep = ""))
    }

    pred_opposite <- predict(model, oppositeFullData)
    treatment_effects <- pred - pred_opposite
    treatment_effects[, fullData[, dim]==0] <- -treatment_effects[, fullData[, dim]==0] 
    meanValue <- mean(treatment_effects)
    # standardDeviation[i] <- sd(trainData[, dim+1]) 
    standardDeviation <- sum((rowMeans(treatment_effects) - meanValue) ^ 2) / (nrow(treatment_effects) - 1)
    # remove newly added value from candidate set
    return(list("meanValueBART" = meanValue, "standardDeviationBART" = standardDeviation, "trainData" = trainData))
  }
  
  # generate extra training data using the scheme (see pdf)
  for (i in 1:num_iterations) {

    # set seed to enable reproduction of the results
    print(c("BART: Epoch=", i))
    # first build BART model
    sink("/dev/null")
    model <- bart(trainData[, 1:dim], trainData[, dim + 1], keeptrees = TRUE, keepevery = 3L, nskip = 1000, ndpost = 5000, ntree = 50, k = 2, usequant = FALSE)
    # model <- wbart(bartModelMatrix(trainData[,1:dim]), trainData[,(dim+1)], bartModelMatrix(candidateX), keepevery=3L, nskip=500L, ndpost=5000L, ntree=50)
    sink()

    # predict the values
    fValues <- predict(model, candidateX)

    # select the best candidate, find its response
    var <- colVars(fValues) * prob_density(candidateX, linear)
    index <- sample(which(var == max(var)), 1)
    response <- candidateY[index]

    # add new data to train set
    trainData <- rbind(trainData, cbind(candidateX[index,], response))
    # remove newly added value from candidate set
    candidateX <- candidateX[-index,]
    candidateY <- candidateY[-index]
    
    # Integral with respect to \Pi_n
    pred <- predict(model, fullData)
    if (save_posterior == TRUE) {
      posterior_samples <- list("posterior_samples" = rowMeans(pred))
      save(posterior_samples, file = paste(save_posterior_dir, "/posterior_BART_synthetic_%s_%s" %--% c(i, num_cv), ".RData", sep = ""))
    }
    pred_opposite <- predict(model, oppositeFullData)
    treatment_effects <- pred - pred_opposite
    treatment_effects[, fullData[, dim]==0] <- -treatment_effects[, fullData[, dim]==0] 
    meanValue[i] <- mean(treatment_effects)
    # standardDeviation[i] <- sd(trainData[, dim+1]) 
    standardDeviation[i] <- sum((rowMeans(treatment_effects) - meanValue[i]) ^ 2) / (nrow(treatment_effects) - 1)
  }
  
  # set seed to enable reproduction of the results
  print(c("BART: Epoch=", i+1))
  # first build BART model
  sink("/dev/null")
  model <- bart(trainData[, 1:dim], trainData[, dim + 1], keeptrees = TRUE, keepevery = 3L, nskip = 1000, ndpost = 5000, ntree = 50, k = 2, usequant = FALSE)
  # model <- wbart(bartModelMatrix(trainData[,1:dim]), trainData[,(dim+1)], bartModelMatrix(candidateX), keepevery=3L, nskip=500L, ndpost=5000L, ntree=50)
  sink()
  # Integral with respect to \Pi_n
  pred <- predict(model, fullData)
  if (save_posterior == TRUE) {
    posterior_samples <- list("posterior_samples" = rowMeans(pred))
    save(posterior_samples, file = paste(save_posterior_dir, "/posterior_BART_synthetic_%s_%s" %--% c(i, num_cv), ".RData", sep = ""))
  }
  pred_opposite <- predict(model, oppositeFullData)
  treatment_effects <- pred - pred_opposite
  treatment_effects[,fullData[, dim]==0] <- -treatment_effects[,fullData[, dim]==0] 
  meanValue[i+1] <- mean(treatment_effects)
  # standardDeviation[i] <- sd(trainData[, dim+1]) 
  standardDeviation[i+1] <- sum((rowMeans(treatment_effects) - meanValue[i]) ^ 2) / (nrow(treatment_effects) - 1)

  return(list("meanValueBART" = meanValue, "standardDeviationBART" = standardDeviation, "trainData" = trainData))
}