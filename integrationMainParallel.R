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
set.seed(0)

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
num_iterations <- as.double(args[2])
num_iterations <- num_iterations*dim
whichGenz <- as.double(args[3])
whichKernel <- as.character(args[4])
# turn on/off sequential design
# 1 denotes TRUE to sequential
# 0 denotes FALSE to sequential
cat("\nBegin testing:\n")
if (as.double(args[5]) == 1 | is.na(as.double(args[5]))) {
  sequential <- TRUE
} else {
  sequential <- FALSE
}
cat("Sequantial design set to", sequential, "\n")
# prior measure over the inputs
# uniform by default
if (as.character(args[6]) != "gaussian" | is.na(args[6])) {
  measure <- "uniform"
} else{
  measure <- as.character(args[6])
}
cat("Prior measure:", measure, "\n")

# extra parameter for step function
# 1 by default
if (whichGenz == 7 & is.na(args[7])) {
  jumps <- 1
  cat("Number of jumps for step function:", jumps, "\n")
} else if (whichGenz == 7){
  jumps <- as.double(args[7])
  cat("Number of jumps for step function:", jumps, "\n")
}

# extra parameter for additive Gaussian function
if (whichGenz == 9){ add_gauss_a <- NA}


print(c(dim, num_iterations, whichGenz))
source("src/genz/genz.R") # genz function to test

if (whichGenz < 1 | whichGenz > 9) { stop("undefined genz function. Change 3rd argument to 1-9") }

if (whichGenz == 1) { genz <- cont; genzFunctionName <-  deparse(substitute(cont)) }
if (whichGenz == 2) { genz <- copeak; genzFunctionName <-  deparse(substitute(copeak)) }
if (whichGenz == 3) { genz <- disc; genzFunctionName <-  deparse(substitute(disc)) }
if (whichGenz == 4) { genz <- gaussian; genzFunctionName <-  deparse(substitute(gaussian)) }
if (whichGenz == 5) { genz <- oscil; genzFunctionName <-  deparse(substitute(oscil)) }
if (whichGenz == 6) { genz <- prpeak; genzFunctionName <-  deparse(substitute(prpeak)) }
if (whichGenz == 7) { genz <- function(xx){return(step(xx, jumps=jumps))}; genzFunctionName <-  deparse(substitute(step)) }
if (whichGenz == 8) { genz <- mix; genzFunctionName <-  deparse(substitute(mix)) }
if (whichGenz == 9) { genz <- function(xx){return(additive_gaussian(xx, a=add_gauss_a))}; genzFunctionName <-  deparse(substitute(additive_gaussian)) }

print("Testing with: %s" %--% genzFunctionName)

num_cv_start <- as.double(args[8])
num_cv_end <- as.double(args[9])

# save posterior samples
if (as.double(args[10]) == 1 | !is.na(as.double(args[8]))) {
  save_posterior <- TRUE
} else {
  save_posterior <- FALSE
}

# prepare training dataset
if (measure == "uniform") {
  trainX <- replicate(dim, runif(20*dim))
  trainY <- genz(trainX)
} else if (measure == "gaussian") {
  trainX <- replicate(dim, rtnorm(20*dim, mean=0.5, lower=0, upper=1))
  trainY <- genz(trainX)
}

for (num_cv in num_cv_start:num_cv_end) {
  # set new seed
  set.seed(num_cv)
  cat("NUM_CV", num_cv, "\n")
  # Bayesian Quadrature method
  # set number of new query points using sequential design
  source("src/BARTBQ.R")
  t0 <- proc.time()
  predictionBART <- mainBARTBQ(dim, num_iterations, FUN = genz, trainX, trainY, sequential, measure, save_posterior=save_posterior, save_posterior_filename=paste(genzFunctionName, num_cv, sep="_"))
  t1 <- proc.time()
  bartTime <- (t1 - t0)[[1]]
  
  # Bayesian Quadrature with Monte Carlo integration method
  print("Begin Monte Carlo Integration")
  source("src/monteCarloIntegration.R")
  
  t0 <- proc.time()
  predictionMonteCarlo <- monteCarloIntegrationUniform(FUN = genz, trainX, trainY, numSamples=num_iterations, dim, measure)
  t1 <- proc.time()
  
  MITime <- (t1 - t0)[[1]]
  
  # Bayesian Quadrature with Gaussian Process
  print("Begin Gaussian Process Integration")
  library(reticulate)
  source("src/optimise_gp.R")
  lengthscale <- optimise_gp_r(trainX, trainY, kernel = whichKernel, epochs=500)
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
    FUN = genz, 
    lengthscale,
    sequential, 
    measure
  )  
  t1 <- proc.time()
  GPTime <- (t1 - t0)[[1]]
  
  # Read in analytical integrals
  source("src/genz/analyticalIntegrals.R")
  dimensionsList <- c(1,2,3,5,10,20)
  whichDimension <- which(dim == dimensionsList)
  if (whichGenz <= 6){
    analyticalIntegrals <- read.csv("results/genz/integrals.csv", header = FALSE)
    real <- analyticalIntegrals[whichGenz, whichDimension]
  } else if (whichGenz == 7) {
    real <- stepIntegral(dim, jumps)
  } else if (whichGenz == 8) {
    if (dim ==1){ real <- 0.008327796}
    if (dim ==2){ real <- 0.008327796 * 2}
    if (dim ==3){ real <- 0.008327796 * 3}
  } else if (whichGenz == 9) {
    real <- additiveGaussianIntegral(dim, a = add_gauss_a)
  }
  
  # Bayesian Quadrature methods: with BART, Monte Carlo Integration and Gaussian Process respectively
  print("Final Results:")
  print(c("Actual integral:", real))
  print(c("BART integral:", predictionBART$meanValueBART[num_iterations]))
  print(c("MI integral:", predictionMonteCarlo$meanValueMonteCarlo[num_iterations]))
  print(c("GP integral:", predictionGPBQ$meanValueGP[num_iterations]))
  
  print("Writing full results to results/genz/%s" %--% c(whichGenz))
  results <- data.frame(
    "epochs" = c(1:num_iterations),
    "BARTMean" = predictionBART$meanValueBART, "BARTsd" = predictionBART$standardDeviationBART,
    "MIMean" = predictionMonteCarlo$meanValueMonteCarlo, "MIsd" = predictionMonteCarlo$standardDeviationMonteCarlo,
    "GPMean" = predictionGPBQ$meanValueGP, "GPsd" = sqrt(predictionGPBQ$varianceGP),
    "actual" = rep(real, num_iterations),
    "runtimeBART" = rep(bartTime, num_iterations),
    "runtimeMI" = rep(MITime, num_iterations),
    "runtimeGP" = rep(GPTime, num_iterations)
  )
  results_models <- list("BART"=predictionBART, "GP"=predictionGPBQ, "MC"=predictionMonteCarlo)
  if (!sequential){
    csvName <- "results/genz/%s/%sDim%sNoSequential%s_%s.csv" %--% c(
      whichGenz, 
      genzFunctionName,
      dim,
      tools::toTitleCase(measure),
      num_cv
    )
    figName <- "Figures/%s/%sDim%sNoSequential%s_%s.pdf" %--% c(
      whichGenz,
      genzFunctionName,
      dim,
      tools::toTitleCase(measure),
      num_cv
    )
    save(results_models, file = "results/genz/%s/%sDim%sNoSequential%s_%s.RData" %--% c(
      whichGenz,
      genzFunctionName,
      dim,
      tools::toTitleCase(measure),
      num_cv
    ))
  } else {
    csvName <- "results/genz/%s/%sDim%s%s_%s.csv" %--% c(
      whichGenz, 
      genzFunctionName,
      dim,
      tools::toTitleCase(measure),
      num_cv
    )
    figName <- "Figures/%s/%sDim%s%s_%s.pdf" %--% c(
      whichGenz,
      genzFunctionName,
      dim,
      tools::toTitleCase(measure),
      num_cv
    )
    save(results_models, file = "results/genz/%s/%sDim%s%s_%s.RData" %--% c(
      whichGenz,
      genzFunctionName,
      dim,
      tools::toTitleCase(measure),
      num_cv
    ))
  }
  print(csvName)
  write.csv(results, file = csvName, row.names=FALSE)
}
