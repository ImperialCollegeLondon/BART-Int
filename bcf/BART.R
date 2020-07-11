library(lhs)
library(dbarts)
library(data.tree)
library(BART)
library(matrixStats)
library(mvtnorm)

prob_density <- function(x) {
  prob <- mvtnorm::dmvnorm(x[, c(1, 2, 3)])
  prob <- prob / 3
  prob[x[, 4] == 0] <- prob[x[, 4] == 0] * pnorm(1)
  prob[x[, 4] == 1] <- prob[x[, 4] == 1] * (1 - pnorm(1))
  return (prob)
}

computeBARTWeighted <- function(trainX, trainY, candidateX, candidateY, num_iterations, save_posterior = FALSE, save_posterior_dir = "bcf/synthetic", num_cv = "default") {
  dim <- ncol(trainX)
  print(c("Adding number of new training data:", num_iterations))
  # outputs
  meanValue <- rep(0, num_iterations)
  standardDeviation <- rep(0, num_iterations)
  trainData <- cbind(trainX, trainY)
  fullData <- rbind(trainX, candidateX)

  colnames(trainData)[dim + 1] <- "response"

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
    fValues <- predict(model, candidateX) * prob_density(candidateX)

    # select the best candidate, find its response
    var <- colVars(fValues)
    index <- sample(which(var == max(var)), 1)
    response <- candidateY[index]

    # add new data to train set
    trainData <- rbind(trainData, c(candidateX[index,], response))
    # Integral with respect to \Pi_n
    pred <- predict(model, fullData)
    if (save_posterior == TRUE) {
      posterior_samples <- list("posterior_samples" = rowMeans(pred))
      save(posterior_samples, file = paste(save_posterior_dir, "/posterior_BART_synthetic_%s_%s" %--% c(i, num_cv), ".RData", sep = ""))
    }

    meanValue[i] <- mean(pred)
    # standardDeviation[i] <- sd(trainData[, dim+1]) 
    standardDeviation[i] <- sum((rowMeans(pred) - meanValue[i]) ^ 2) / (nrow(pred) - 1)

    # remove newly added value from candidate set
    candidateX <- candidateX[-index,]
    candidateY <- candidateY[-index]

  }

  return(list("meanValueBART" = meanValue, "standardDeviationBART" = standardDeviation, "trainData" = trainData))
}