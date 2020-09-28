###
# Functions for BART
###
library(lhs)
library(dbarts)
library(data.tree)
library(matrixStats)
library(mvtnorm)
library(msm)
terminalProbability <- function(currentNode)

#'Terminal Node Probability
#' 
#'@description This function calculates the probability assigned to the 
#' input terminal node.
#' 
#'@param currentNode Node; Terminal node of the tree
#'
#'@return The probability assigned to the input terminal node.

{
  prob <- currentNode$probability

  while (!isRoot(currentNode$parent)) {
    currentNode <- currentNode$parent
    prob <- prob * currentNode$probability
  }

  return(prob)
}

fillProbabilityForNode <- function(oneTree, cutPoints, cut, measure)

#'Fill Non-Terminal Node Probability
#' 
#'@description This function calculates the prior probability of all the non-terminal tree nodes 
#'and fill that node with the prior probability.
#' 
#'@param currentNode Node; Non-terminal nodes in the tree
#'
#'@return The non-terminal tree nodes filled with prior probability.

{

  if (!is.null(oneTree$leftChild)) {

    decisionRule <- cutPoints[[oneTree$splitVar]][oneTree$splitIndex]

    if (measure == "uniform") {
      oneTree$leftChild$probability <- (decisionRule - cut[1, oneTree$splitVar]) / (cut[2, oneTree$splitVar] - cut[1, oneTree$splitVar])
      oneTree$rightChild$probability <- (cut[2, oneTree$splitVar] - decisionRule) / (cut[2, oneTree$splitVar] - cut[1, oneTree$splitVar])
    } else if (measure == "gaussian") {
      normalizingConst <- pmvnorm(cut[1, oneTree$splitVar], cut[2, oneTree$splitVar], mean = 0.5, sigma = 1)
      oneTree$leftChild$probability <- pmvnorm(cut[1, oneTree$splitVar], decisionRule, sigma = 1) / normalizingConst
      oneTree$rightChild$probability <- 1 - oneTree$leftChild$probability
    } else if (measure == "exponential") {
      normalizingConst <- pexp(cut[1, oneTree$splitVar]) - pexp(cut[2, oneTree$splitVar])
      oneTree$leftChild$probability <- pexp(decisionRule) - pexp(cut[1, oneTree$splitVar]) / normalizingConst
      oneTree$rightChild$probability <- 1 - oneTree$leftChild$probability      
    }

    cut[, oneTree$splitVar] = c(0, decisionRule)

    fillProbabilityForNode(oneTree$leftChild, cutPoints, cut, measure)

    cut[, oneTree$splitVar] = c(decisionRule, 1)

    fillProbabilityForNode(oneTree$rightChild, cutPoints, cut, measure)

  } else if (is.null(oneTree$probability)) {

    oneTree$probability <- 1

  }

  return(oneTree)
}

terminalProbabilityStore <- function(Tree)
#'Fill Terminal Node Probability
#' 
#'@description This function store probabilities into all terminal nodes of the tree.
#'@param currentNode Node; Root node of a tree
#'
#'@return Root node with terminal-node probabilities stored.
{
  terminalNodes = Traverse(Tree, filterFun = isLeaf)

  for (i in 1:length(terminalNodes)) {
    probability2 <- terminalProbability(terminalNodes[[i]])
    terminalNodes[[i]]$terminal_probability <- probability2
  }

  return(Tree)
}

getTree <- function(sampler, chainNum, sampleNum, treeNum)

#'Build tree
#' 
#'@description The function builds tree and fills attributes to each node of the tree
#' 
#'@param sampler List; Bart model
#'@param chainNum Integer; The index of the chain 
#'@param sampleNum Integer; The index of the posterior sample draw
#'@param treenNum Integer; The index of the tree sample extracted
#'
#'@return A tree sample in BART fitting at the position specifized
#'by chainNum, sampleNUm and treeNum, with attributeds filled into
#'each node of the tree

{
  cutPoints <- dbarts:::createCutPoints(sampler)
  if (sampler$control@keepTrees) {
    treeString <- sampler$state[[chainNum]]@savedTrees[treeNum, sampleNum]
    treeFits <- sampler$state[[chainNum]]@savedTreeFits[, treeNum, sampleNum]
  }
  else {
    treeString <- sampler$state[[chainNum]]@trees[treeNum]
    treeFits <- sampler$state[[chainNum]]@treeFits[, treeNum]
  }

  tree <- dbarts:::buildTree(strsplit(gsub("\\.", "\\. ", treeString),
                                      " ", fixed = TRUE)[[1]])
  tree$remainder <- NULL

  tree$indices <- seq_len(length(sampler$data@y))
  tree <- dbarts:::fillObservationsForNode(tree, sampler, cutPoints)

  tree <- dbarts:::fillPlotInfoForNode(tree, sampler, treeFits)
  maxDepth <- dbarts:::getMaxDepth(tree)

  tree <- dbarts:::fillPlotCoordinatesForNode(tree, maxDepth, 1L, 1L)
  numEndNodes <- tree$index - 1L
  tree$index <- NULL

  return(tree)
}

singleTreeSum <- function(treeNum, model, drawNum, dim, measure)

#'Sum over Terminal Nodes
#' 
#'@description This function sums over the product of a tree terminal node's mu value
#'and probability. Tree's position is specified by treeNum and drawNum
#' 
#'@param treeNum Integer; The index of tree sample extracted
#'@param moodel List; The BART fitting model
#'@drawNum Integer; The index of the posterior draw
#'@param dim Integer; Dimension of inputs X
#'
#'@return Real; Sum over all nodes of a single tree

{
  cutPoints <- dbarts:::createCutPoints(model$fit)
  cut <- array(c(0, 1), c(2, dim))
  
  if (measure %in% c("exponential")) {
    cut <- array(c(0, Inf), c(2, dim))
  }

  treeList <- getTree(model$fit, 1, drawNum, treeNum)

  selectedTree <- FromListSimple(treeList)

  #Modify tree by the functions written above
  selectedTree <- fillProbabilityForNode(selectedTree, cutPoints, cut, measure)
  selectedTree <- terminalProbabilityStore(selectedTree)


  terminalNodeList <- Traverse(selectedTree, filterFun = isLeaf)

  #Calculate approximation of integreal in the single tree 
  integral <- 0
  for (node in terminalNodeList) {
    #We use the mean of prediction value Y's in the terminal node as u
    #rescale this step
    integral <- integral + node$terminal_probability * node$mu
  }
  return(integral)
}

posteriorSum <- function(drawNum, model, dim, measure)

#'Sum over Trees
#' 
#'@description This function sums the output of singleTreeSum
#'over all the tree of the posterior draw specified by drawNum. 
#' 
#'@drawNum Integer; The index of the posterior draw
#'@param moodel List; The BART fitting model
#'@param dim Integer; Dimension of inputs X
#'
#'@return Real; Sum over a posterior draw

{
  nTree <- ncol(model$fit$state[[1]]@treeFits)
  treeNum <- seq(1, nTree, length.out = nTree)

  #Extra variables
  var <- list(model, drawNum, dim, measure)

  #Calculate integration over all trees in the draw by mapply
  integral <- sum(unlist(mapply(singleTreeSum, treeNum, MoreArgs = var)))

  return(integral)
}


sampleIntegrals <- function(model, dim, measure)

#'Sum over Posterior Draws
#' 
#'@description This function sums the values of posteriorSum over all posterior draws 
#' 
#'@param moodel List; The BART fitting model
#'@param dim Integer; Dimension of inputs X
#'
#'@return Real; Sum over a posterior draw

{
  nDraw <- dim(model$fit$state[[1]]@treeFits)[2]
  drawNum <- seq(1, nDraw, length.out = nDraw)

  #Extra Variables
  var <- list(model, dim, measure)
  integrals <- mapply(posteriorSum, drawNum, MoreArgs = var)
  return(integrals)
}

BARTBQSequential <- function(dim, trainX, trainY, numNewTraining, FUN, sequential, measure, save_posterior = FALSE, save_posterior_dir = "results/genz", save_posterior_filename = "default")

#'BART-BQ with Sequential Design
#' 
#'@description This function approxiamtes the integration of target function using
#'BART and Sequential Design.
#' 
#'@param dim Integer; Dimension of inputs X
#'@param trainX Matrix; covariates of training data
#'@param trainY Numeric Arraay; response of training data
#'@param numNewTraining Integer; number of new training points to be added
#'@param FUN Function; function that we are integrating over
#'@param sequential boolean; whether or not to use sequential design
#'
#'@return List; list of mean integral value, standard deviation of integral value and new traiing set

{

  print(c("Adding number of new training data:", numNewTraining))
  # outputs
  meanValue <- rep(0, numNewTraining)
  standardDeviation <- rep(0, numNewTraining)
  trainData <- data.frame(trainX, trainY)
  print(paste("BART: Epoch = ", 1))
  # find the min and max range of y
  ymin <- min(trainData[, (dim + 1)]);
  ymax <- max(trainData[, (dim + 1)])
  # first build BART and scale mean and standard deviation
  sink("/dev/null")
  # model <- bart(trainData[,1:dim], trainData[,dim+1], keeptrees=TRUE, keepevery=5L, nskip=500, ndpost=1500, ntree=50, k = 2)
  # sigest = 0.0001 ensures that there is almost no noise
  model <- bart(trainData[, 1:dim], trainData[, dim + 1], keeptrees = TRUE, keepevery = 5L, nskip = 500, ndpost = 1500, ntree = 200, k = 2, sigest=0.0001)
  sink()
  # obtain posterior samples
  integrals <- sampleIntegrals(model, dim, measure)
  integrals <- (integrals + 0.5) * (ymax - ymin) + ymin
  meanValue[1] <- mean(integrals)
  standardDeviation[1] <- sqrt(sum((integrals - meanValue[1]) ^ 2) / (length(integrals) - 1))
  if (save_posterior == TRUE) {
    posterior_samples <- list("posterior_samples" = integrals)
    save(posterior_samples, file = paste(save_posterior_dir, "/posterior_BART_Dim%s_%s_1" %--% c(dim, save_posterior_filename), ".RData", sep = ""))
  }
  if (numNewTraining == 1) {
    return(list("meanValueBART" = meanValue, "standardDeviationBART" = standardDeviation,
                "trainData" = trainData, "model" = model))
  }

  # generate extra training data using the scheme (see pdf)
  for (i in 2:numNewTraining) {

    print(paste("BART: Epoch =", i))
    # find the min and max range of y
    ymin <- min(trainData[, (dim + 1)]);
    ymax <- max(trainData[, (dim + 1)])
    # first build BART and scale mean and standard deviation
    sink("/dev/null")
    # sigest = 0.0001 ensures that there is almost no noise
    model <- bart(trainData[, 1:dim], trainData[, dim + 1], keeptrees = TRUE, keepevery = 5L, nskip = 500, ndpost = 1500, ntree = 200, k = 2, sigest=0.0001)
    # model <- bart(trainData[,1:dim], trainData[,dim+1], keeptrees=TRUE, keepevery=5L, nskip=500, ndpost=1500, ntree=50, k = 2)
    # model <- bart(trainData[,1:dim], trainData[,dim+1], keeptrees=TRUE, keepevery=20L, nskip=1000, ndpost=2000)
    # model <- bart(trainData[,1:dim], trainData[,dim+1], keeptrees=TRUE, keepevery=20L, nskip=1000, ndpost=10000, ntree=50, k = 2)

    sink()
    # obtain posterior samples
    integrals <- sampleIntegrals(model, dim, measure)
    integrals <- (integrals + 0.5) * (ymax - ymin) + ymin
    meanValue[i] <- mean(integrals)
    standardDeviation[i] <- sqrt(sum((integrals - meanValue[i]) ^ 2) / (length(integrals) - 1))
    if (save_posterior == TRUE) {
      posterior_samples <- list("posterior_samples" = integrals)
      save(posterior_samples, file = paste(save_posterior_dir, "/posterior_BART_Dim%s_%s_%s" %--% c(dim, save_posterior_filename, i), ".RData", sep = ""))
    }

    # sequential design section, where we build the new training data
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

    # predict the values
    fValues <- predict(model, candidateSet)
    if (sequential) {
      var <- colVars(fValues) * weights
      index <- sample(which(var == max(var)), 1)
    }
    else {
      index <- sample(1:candidateSetNum, 1)
    }
    value <- FUN(t(candidateSet[index,]))
    trainData <- rbind(trainData, c(candidateSet[index,], value))

  }

  return(list("meanValueBART" = meanValue, "standardDeviationBART" = standardDeviation,
               "trainData" = trainData, "model" = model))
}

BART_posterior <- function(dim, trainX, trainY, numNewTraining, FUN, sequential, measure = "uniform")

#'BART-BQ with Sequential Design
#' 
#'@description This function approxiamtes the integration of target function using
#'BART and Sequential Design.
#' 
#'@param dim Integer; Dimension of inputs X
#'@param trainX Matrix; covariates of training data
#'@param trainY Numeric Arraay; response of training data
#'@param numNewTraining Integer; number of new training points to be added
#'@param FUN Function; function that we are integrating over
#'@param sequential boolean; whether or not to use sequential design
#'
#'@return List; list of mean integral value, standard deviation of integral value and new traiing set

{

  print(c("Adding number of new training data:", numNewTraining))
  # outputs
  meanValue <- rep(0, numNewTraining)
  standardDeviation <- rep(0, numNewTraining)
  trainData <- data.frame(trainX, trainY)

  # generate extra training data using the scheme (see pdf)
  for (i in 1:numNewTraining) {
    print(paste("BART: Epoch =", i))
    # find the min and max range of y
    ymin <- min(trainData[, (dim + 1)]);
    ymax <- max(trainData[, (dim + 1)])
    # first build BART and scale mean and standard deviation
    sink("/dev/null")
    # sigest = 0.0001 ensures that there is almost no noise
    model <- bart(trainData[, 1:dim], trainData[, dim + 1], keeptrees = TRUE, keepevery = 5L, nskip = 1000, ndpost = 5000, ntree = 50, k = 2, sigest = 0.0001)
    sink()
    # obtain posterior samples
    integrals <- sampleIntegrals(model, dim, measure)
    integrals <- (integrals + 0.5) * (ymax - ymin) + ymin

    # sequential design section, where we build the new training data
    candidateSetNum <- 100
    if (measure == "uniform") {
      candidateSet <- replicate(dim, runif(candidateSetNum))
    } else if (measure == "gaussian") {
      candidateSet <- replicate(dim, rtnorm(candidateSetNum, mean = 0.5, lower = 0, upper = 1))
    } else if (measure == "exponential") {
      candidateSet <- replicate(dim, rexp(candidateSetNum))
    }

    # predict the values
    fValues <- predict(model, candidateSet)

    probability = 1 #uniform probability
    #expectedValue <- colMeans(fValues*probability)

    if (sequential) {
      var <- colVars(fValues)
      index <- sample(which(var == max(var)), 1)
    }
    else {
      index <- sample(1:candidateSetNum, 1)
    }
    value <- FUN(t(candidateSet[index,]))
    trainData <- rbind(trainData, c(candidateSet[index,], value))
  }
  return(list("model" = model, "trainData" = trainData))
}

mainBARTBQ <- function(dim, num_iterations, FUN, trainX, trainY, sequential = TRUE, measure = "uniform", save_posterior = FALSE, save_posterior_dir = "results/genz", save_posterior_filename = "default")

#'BART-BQ with Sequential Design
#' 
#'@description This function approxiamtes the integration of target function using
#'BART and Sequential Design.
#' 
#'@param dim Integer; Dimension of inputs X
#'@param num_iterations Integer; Number of new data points
#'@param FUN Function; Integral Function
#'@param trainX Matrix; covariates of training data
#'@param trainY Numeric Arraay; response of training data
#'@param sequential Boolean: whether to use sequential design or not. If FALSE, we simply just increase
#'the number of design points
#'
#'@return List; list of mean integral value, standard deviation of integral value and new traiing set

{
  # prepare training data and parameters
  genz <- FUN #select genz function
  numNewTraining <- num_iterations
  prediction <- BARTBQSequential(dim, trainX, trainY, numNewTraining, FUN = genz, sequential, measure, save_posterior = save_posterior, save_posterior_dir = save_posterior_dir, save_posterior_filename = save_posterior_filename)

  return(prediction)
}
