# !/usr/bin/env R
# Load required packages
library(MASS)
library(cubature)
library(lhs)
library(data.tree)
library(dbarts)
library(matrixStats)
library(mvtnorm)
library(doParallel)
library(kernlab)
library(msm)
library(MCMCglmm)

# define string formatting
`%--%` <- function(x, y)
  # from stack exchange:
  # https://stackoverflow.com/questions/46085274/is-there-a-string-formatting-operator-in-r-similar-to-pythons
{
  do.call(sprintf, c(list(x), y))
}

# global parameters: dimension
args <- commandArgs(TRUE)
dim <- as.double(args[1])
num_iterations <- 1
whichRare <- as.double(args[2])
whichKernel <- as.character(args[3])
# turn on/off sequential design
# 1 denotes TRUE to sequential
# 0 denotes FALSE to sequential
cat("\nBegin testing:\n")
sequential <- FALSE
cat("Sequantial design set to", sequential, "\n")
# prior measure over the inputs
# uniform by default
measure="exponential"
cat("Integral Measure:", measure, "\n")


# save posterior samples
if (as.double(args[5]) == 1 & !is.na(args[5])) {
  save_posterior <- TRUE
} else {
  save_posterior <- FALSE
}
num_cv <- as.double(args[6])
print(c(dim, num_iterations, whichRare))
source("src/rareFunctions.R") # rare function to test

if (whichRare < 1 | whichRare > 3) { stop("undefined rare function. Change 3rd argument to 1") }
if (whichRare == 1) { 
  rareFunction <- function(xx) { return(indicator_greater(xx, threshold = 3)) }
  rareFunctionName <- deparse(substitute(indicator_greater))
}
if (whichRare == 2) { 
  rareFunction <- function(xx) { return(indicator_square_greater(xx, threshold = 5)) }
  rareFunctionName <- deparse(substitute(indicator_square_greater))
}

if (whichRare == 3) {
  rareFunction <- function(xx) { return(portfolio_loss(xx, gamma=2)) }
  rareFunctionName <- deparse(substitute(portfolio_loss))
}

print("Testing with: %s" %--% rareFunctionName)

# prepare training dataset
if (measure == "uniform") {
  trainX <- replicate(dim, runif(500 * dim))
  trainY <- rareFunction(trainX)
} else if (measure == "gaussian") {
  trainX <- replicate(dim, rtnorm(500 * dim, mean = 0.5, lower = 0, upper = 1))
  trainY <- rareFunction(trainX)
} else if (measure == "exponential") {
  trainX <- replicate(dim, rexp(500 * dim))
  trainY <- rareFunction(trainX)
}

for (num_cv in num_cv:num_cv) {
  # set new seed
  set.seed(num_cv)
  cat("NUM_CV", num_cv, "\n")
  # BART-Int method
  # set number of new query points using sequential design
  source("src/BARTInt.R")
  t0 <- proc.time()
  predictionBART <- mainBARTInt(
    dim, 
    num_iterations,
    FUN = rareFunction,
    trainX,
    trainY,
    sequential,
    measure,
    save_posterior = save_posterior,
    save_posterior_dir = "results/rares",
    save_posterior_filename = paste(rareFunctionName, num_cv, sep = "_")
  )
  t1 <- proc.time()
  bartTime <- (t1 - t0)[[1]]
  
  # Bayesian Quadrature with Monte Carlo integration method
  print("Begin Monte Carlo Integration")
  source("src/monteCarloIntegration.R")
  
  t0 <- proc.time()
  predictionMonteCarlo <- monteCarloIntegrationUniform(
    FUN = rareFunction, 
    trainX, 
    trainY, 
    numSamples = num_iterations, 
    dim, 
    measure
  )
  # t1 <- proc.time()
  MITime <- (t1 - t0)[[1]]
  
  # Bayesian Quadrature with Gaussian Process
  print("Begin Gaussian Process Integration")
  library(reticulate)
  source("src/optimise_gp.R")
  lengthscale <- optimise_gp_r(trainX, trainY, kernel = whichKernel, epochs = 500)
  print("...Finished training for the lengthscale")
  
  source("src/GPBQ.R")
  t0 <- proc.time()
  # need to add in function to optimise the hyperparameters
  predictionGPBQ <- computeGPBQ_matern(
    trainX,
    trainY,
    dim,
    epochs = num_iterations,
    kernel = whichKernel,
    FUN = rareFunction,
    lengthscale,
    sequential,
    measure
  )
  t1 <- proc.time()
  GPTime <- (t1 - t0)[[1]]
  
  # Bayesian Quadrature methods: with BART, Monte Carlo Integration and Gaussian Process respectively
  print("Final Results:")
  print(c("BART integral:", predictionBART$meanValueBART[num_iterations]))
  print(c("MI integral:", predictionMonteCarlo$meanValueMonteCarlo[num_iterations]))
  print(c("GP integral:", predictionGPBQ$meanValueGP[num_iterations]))
  
  print("Writing full results to results/rares/%s" %--% c(whichRare))
  results <- data.frame(
    "epochs" = c(1:num_iterations),
    "BARTMean" = predictionBART$meanValueBART, "BARTsd" = predictionBART$standardDeviationBART,
    "MIMean" = predictionMonteCarlo$meanValueMonteCarlo[1:num_iterations], "MIsd" = predictionMonteCarlo$standardDeviationMonteCarlo[1:num_iterations],
    "GPMean" = predictionGPBQ$meanValueGP, "GPsd" = sqrt(predictionGPBQ$varianceGP),
    "runtimeBART" = rep(bartTime, num_iterations),
    "runtimeMI" = rep(MITime, num_iterations),
    "runtimeGP" = rep(GPTime, num_iterations)
  )
  results_models <- list("BART" = predictionBART, "GP" = predictionGPBQ, "MC" = predictionMonteCarlo)
  if (!sequential) {
    csvName <- "results/rares/%s/%sDim%sNoSequential%s_%s.csv" %--% c(
      whichRare,
      rareFunctionName,
      dim,
      tools::toTitleCase(measure),
      num_cv
    )
    save(results_models, file = "results/rares/%s/%sDim%sNoSequential%s_%s.RData" %--% c(
      whichRare,
      rareFunctionName,
      dim,
      tools::toTitleCase(measure),
      num_cv
    ))
  } else {
    csvName <- "results/rares/%s/%sDim%s%s_%s.csv" %--% c(
      whichRare,
      rareFunctionName,
      dim,
      tools::toTitleCase(measure),
      num_cv
    )
    save(results_models, file = "results/rares/%s/%sDim%s%s_%s.RData" %--% c(
      whichRare,
      rareFunctionName,
      dim,
      tools::toTitleCase(measure),
      num_cv
    ))
  }
  
  write.csv(results, file = csvName, row.names = FALSE)
}
