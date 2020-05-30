library(lhs)
library(dbarts)
library(data.tree)
library(matrixStats)
library(caret)
set.seed(0)
source("src/survey_design/bartMean.R")
source("src/survey_design/gpMean.R")
source("src/optimise_gp.R")

# paths to save results and plots
resultPath <- "results/survey_design/"
plotPath <- "Figures/survey_design/"

args <- as.double(commandArgs(TRUE))
num_cv_start <- args[1]
num_cv_end <- args[2]
num_data <- args[3] # set to 2000 for this lengthscale
num_design = args[4]
which_response <- args[5] # set to 1 for binary log-income response
jitter <- 1e-6

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
trainData_full <- trainData[complete.cases(trainData),]
candidateData <- candidateData[complete.cases(candidateData),]

# change response to binary
if (which_response == 1) {
  binary_response <- ifelse(trainData_full$log_Total_person_income >= 10, 1 - jitter, jitter)
  trainData_full$log_Total_person_income <- binary_response
  binary_response <- ifelse(candidateData$log_Total_person_income >= 10, 1 - jitter, jitter)
  candidateData$log_Total_person_income <- binary_response
} else {
}

# compute the real population mean log income
bart_ground_truths <- c()
mi_ground_truths <- c()

for (num_cv in num_cv_start:num_cv_end) {

  # set new seed
  set.seed(num_cv)
  print(num_cv)
  trainData <- trainData_full[sample(c(1:dim(trainData)[1]), num_design),]
  candidateData <- candidateData[1:num_data,]
  poptMean <- mean(c(trainData$log_Total_person_income, candidateData$log_Total_person_income))

  dim <- ncol(trainData)
  allData <- rbind(trainData, candidateData)
  model <- bart(allData[, 1:(dim - 1)], allData[, dim], keeptrees = TRUE, keepevery = 3L, nskip = 1000, ndpost = 10000, ntree = 50, k = 3, usequant = FALSE)
  BARTpoptMean <- mean(model$yhat.train.mean)
  rm(model)
  gc()

  mi_ground_truths[num_cv] <- poptMean
  bart_ground_truths[num_cv] <- BARTpoptMean
}
ground_truths <- data.frame(mi_ground_truths)
ground_truths <- cbind(ground_truths, bart_ground_truths)
write.csv(ground_truths, paste("results/survey_design/popt_", num_design, "_", num_data, "_", num_cv, ".csv", sep = ""))



