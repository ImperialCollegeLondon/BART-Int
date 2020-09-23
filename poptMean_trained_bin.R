library(lhs)
library(dbarts)
library(data.tree)
library(matrixStats)
library(EnvStats)
library(ks)
library(caret)
set.seed(0)
source("src/survey_design/bartMean.R")
source("src/survey_design/gpMean.R")
source("src/optimise_gp.R")

build_density <- function(data) {
  # data: non-one-hot encoded dataframe
  # usually this is rbind(nonOneHotTrainX, nonOneHotCandidateX)
  data[, c(1,2,4,5,6,7)] <- sapply(data[, c(1,2,4,5,6,7)], as.character); 
  data$category <- apply(data[, c(1,2,4,5,6,7)], 1, paste, collapse = "-"); 
  data$category <- sapply(sapply(data$category, as.factor), as.numeric)
  
  density_func <- function(x) {
    # x is a non-one-hot encoded dataframe
    p_education <- predict(ks::kde(data$Education, h=1), x=x$Education)
    
    x[, c(1,2,4,5,6,7)] <- sapply(x[, c(1,2,4,5,6,7)], as.character); 
    x$category <- apply(x[, c(1,2,4,5,6,7)], 1, paste, collapse = "-"); 
    x$category <- sapply(sapply(x$category, as.factor), as.numeric)
    p_discrete <- demp(x$category, data$category)    
    
    return (p_education * p_discrete)
  }
  return (density_func)  
}

# paths to save results and plots
resultPath <- "results/survey_design/"
plotPath <- "Figures/survey_design/"

args <- as.double(commandArgs(TRUE))
num_new_surveys <- args[1]
num_cv_start <- args[2]
num_cv_end <- args[3]
num_data <- args[4] # size of candidate set
num_design <- args[5] # size of initial training set
jitter <- 1e-5
cutoff <- 10


# read in data
trainData <- read.csv("data/train2.csv")
candidateData <- read.csv("data/candidate2.csv")

# convert num to factor, log income
convert <- function(data) {
  log_Total_person_income <- log(data[, ncol(data)])
  data <- sapply(data[, 2:(ncol(data) - 1)], as.factor)
  data <- data.frame(data)
  data <- cbind(data, log_Total_person_income)
}

trainData <- convert(trainData)
candidateData <- convert(candidateData)

trainData <- trainData[!is.infinite(trainData$log_Total_person_income),]
candidateData <- candidateData[!is.infinite(candidateData$log_Total_person_income),]
# encode as continuous
trainData$Education <- as.double(as.character(trainData$Education)) + jitter
candidateData$Education <- as.double(as.character(candidateData$Education)) + jitter
# remove missing values
trainData_full <- trainData[complete.cases(trainData),]
candidateData <- candidateData[complete.cases(candidateData),]
# keep the first num_data observations as candidate set
candidateData <- candidateData[1:num_data,]
# change response to binary
binary_response <- ifelse(trainData_full$log_Total_person_income >= cutoff, 1 - jitter, jitter)
trainData_full$log_Total_person_income <- binary_response
binary_response <- ifelse(candidateData$log_Total_person_income >= cutoff, 1 - jitter, jitter)
candidateData$log_Total_person_income <- binary_response

for (num_cv in num_cv_start:num_cv_end) {
  # set new seed
  set.seed(num_cv)
  print(num_cv)
  trainData <- trainData_full[sample(c(1:dim(trainData_full)[1]), num_design),]

  # extract covariates and response
  cols <- ncol(trainData)
  trainX <- trainData[, - cols]
  trainY <- trainData[, cols]
  candidateX <- candidateData[, - cols]
  candidateY <- candidateData[, cols]
  trainX <- data.frame(trainX)
  candidateX <- data.frame(candidateX)
  # one-hot encoding
  trainX.num <- trainX
  candidateX.num <- candidateX
  trainX.num[, c(1,2,4,5,6,7)] <- sapply(trainX[, c(1,2,4,5,6,7)], as.numeric)
  candidateX.num[, c(1,2,4,5,6,7)] <- sapply(candidateX[, c(1,2,4,5,6,7)], as.numeric)
  dummyFullData <- dummyVars("~.", data = rbind(trainX, candidateX))
  trainX <- data.frame(predict(dummyFullData, newdata = trainX))
  candidateX <- data.frame(predict(dummyFullData, newdata = candidateX))
  # compute population average income estimates by BARTBQ
  t0 <- proc.time()
  BARTresults <- computeBART(trainX.num, trainX, trainY, candidateX.num, candidateX, candidateY, num_iterations = num_new_surveys, save_posterior = TRUE, num_cv = num_cv)
  t1 <- proc.time()
  bartTime <- (t1 - t0)[[1]]
  # population average income estimation by Monte Carlo
  MIresults <- computeMI(trainX.num, trainY, candidateX.num, candidateY, num_iterations = num_new_surveys, seed = num_cv)

  # GPBQ
  lengthscale <- optimise_gp_r(as.matrix(trainX), trainY, kernel = "matern32", epochs = 500)
  t0 <- proc.time()
  GPresults <- computeGPBQEmpirical(trainX.num, as.matrix(trainX), trainY, candidateX.num, as.matrix(candidateX), candidateY, epochs = num_new_surveys, kernel = "matern32", lengthscale = lengthscale)
  t1 <- proc.time()
  GPTime <- (t1 - t0)[[1]]

  # store results
  results <- data.frame(
    "epochs" = c(1:num_new_surveys),
    "BARTMean" = BARTresults$meanValueBART, "BARTsd" = BARTresults$standardDeviationBART,
    "MIMean" = MIresults$meanValueMI, "MIsd" = MIresults$standardDeviationMI,
    "GPMean" = GPresults$meanValueGP, "GPsd" = sqrt(GPresults$varianceGP),
  # "PoptMean" = ground_truths$mi_ground_truths[1], "BpoptMean" = ground_truths$bart_ground_truths[1],
    "runtimeBART" = rep(bartTime, num_new_surveys),
    "runtimeGP" = rep(GPTime, num_new_surveys)
  )
  write.csv(results, file = paste0(resultPath, "results", num_cv, ".csv"), row.names = FALSE)
  results_models <- list("BART" = BARTresults, "MI" = MIresults, "GP" = GPresults)
  save(results_models, file = paste0(resultPath, "results", num_cv, ".RData"))

  print(c("BART-Int: ", results$BARTMean[num_new_surveys]))
  print(c("GP-BQ: ", results$GPMean[num_new_surveys]))
  print(c("MI: ", results$MIMean[num_new_surveys]))
}
