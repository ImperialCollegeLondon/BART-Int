library(lhs)
library(dbarts)
library(data.tree)
library(BART)
library(matrixStats)
# library(docstring)

# define string formatting
`%--%` <- function(x, y) 
# from stack exchange:
# https://stackoverflow.com/questions/46085274/is-there-a-string-formatting-operator-in-r-similar-to-pythons
{
  do.call(sprintf, c(list(x), y))
}

computeBART <- function(nonOneHotTrainX, trainX, trainY, nonOneHotCandidateX, candidateX, candidateY, num_iterations, save_posterior=FALSE, save_posterior_dir="results/survey_design", num_cv="default") 
{
  dim <- ncol(trainX)
  print( c("Adding number of new training data:", num_iterations) )
  # outputs
  meanValue <- rep(0, num_iterations)
  standardDeviation <- rep(0, num_iterations)
  eduMeanValue <- matrix(0, 2*num_iterations, nrow=2, ncol=num_iterations)
  eduStandardDeviation <- matrix(0, 2*num_iterations, nrow=2, ncol=num_iterations)
  trainData <- cbind(trainX, trainY)
  fullData <- rbind(trainX, candidateX)
  emp_density_func <- build_density(rbind(nonOneHotTrainX, nonOneHotCandidateX))
  weights <- emp_density_func(nonOneHotCandidateX)
  
  colnames(trainData)[dim+1] <- "response"

  # generate extra training data using the scheme (see pdf)
  for (i in 1:num_iterations) {
    
    # set seed to enable reproduction of the results
    print(c("BART: Epoch=", i))
    # first build BART model
    sink("/dev/null")
    # model <- bart(trainData[, 1:dim], trainData[, dim+1], keeptrees=TRUE, keepevery=3L, 
    #               nskip=500, ndpost=2000, ntree=50, k=2, usequant=FALSE)
    model <- bart(trainData[,1:dim], trainData[, dim+1], keeptrees=TRUE, keepevery=3L, 
                  nskip=1000, ndpost=5000, ntree=50, k=2, usequant=FALSE)             
    # model <- wbart(bartModelMatrix(trainData[,1:dim]), trainData[,(dim+1)], bartModelMatrix(candidateX), keepevery=3L, nskip=500L, ndpost=5000L, ntree=50)
    sink()

    # predict the values
    fValues <- predict(model, candidateX)

    # select the best candidate, find its response
    var <- colVars(fValues) * weights
    index <- sample(which(var==max(var)), 1)
    response <- candidateY[index]

    # # data segmentation by education level
    # mat <- matrix(c(1,16,17,24), nrow = 2, ncol = 2)

    # for (cat in 1:2){

    #   eduCandidateX <- candidateX[(candidateX$Education >= mat[1, cat] & candidateX$Education <= mat[2, cat]), ]
      
    #   # make prediction
    #   fValues <- predict(model, eduCandidateX)

    #   # posterior mean and variance
    #   integrals <- rowMeans(fValues)
    #   eduMeanValue[cat, i] <- mean(integrals)
    #   eduStandardDeviation[cat, i] <- sd(integrals)

    # }
    
    # add new data to train set
    trainData <- rbind(trainData, cbind(candidateX[index, ], response))
    # Integral with respect to \Pi_n
    pred <- predict(model, fullData)
    if (save_posterior == TRUE) {
      posterior_samples <- list("posterior_samples" = rowMeans(pred))
      save(posterior_samples, file = paste(save_posterior_dir, "/posterior_BART_survey_%s_%s" %--% c(i, num_cv), ".RData", sep=""))
    }

    meanValue[i] <- mean(pred)
    # standardDeviation[i] <- sd(trainData[, dim+1]) 
    standardDeviation[i] <- sum((rowMeans(pred) - meanValue[i])^2) / (nrow(pred) - 1)

    # remove newly added value from candidate set
    candidateX <- candidateX[-index,]
    weights <- weights[-index]
    candidateY <- candidateY[-index]

  }

  return(list("meanValueBART"=meanValue, "standardDeviationBART"=standardDeviation, 
              "eduMeanValueBART"=eduMeanValue, "eduStandardDeviationBART"=eduStandardDeviation, "trainData"=trainData))
}

computeMI <- function(trainX, trainY, candidateX, candidateY, num_iterations, seed=NA)
{

  if (!is.na(seed)){
    set.seed(seed)
    candidateY <- sample(candidateY)
  }
  MImean <- rep(NA, num_iterations)
  combinedY <- c(trainY, candidateY[1:num_iterations])
  combinedN <- length(combinedY)
  
  MImean <- cumsum(combinedY) / (1:combinedN)
  cumvar <- (cumsum(combinedY^2) - MImean^2 * (1:combinedN)) / (0:(combinedN - 1))
  MIstandardDeviation <- sqrt(cumvar) / sqrt(1:combinedN)
  
  MImean <- MImean[-(1:length(trainY))]
  MIstandardDeviation <- MIstandardDeviation[-c(1:length(trainY))]

  return(list("meanValueMI"=MImean, "standardDeviationMI"=MIstandardDeviation))

}