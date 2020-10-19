# Store all estimates and best methods in separate csv's.
# To run:
#     In Unix terminal: 
#     Rscript {ROOT}/src/genz/drawGraphs.R
# To see results:
#     Results stored in {ROOT}/report/genz/resultsAll.csv and bestmethods.csv

methods <- c("BART", "MI", "GP")
bestMethod <- matrix(NA, ncol = 6, nrow = 7)
rownames(bestMethod) <- c("cont", "copeak", "disc", "gaussian", "oscil", "prpeak", "step")
colnames(bestMethod) <- c("1", "2", "3", "5", "10", "20")

# Store all results in a single csv file
resultsAll <- matrix(NA, ncol = 7, nrow = 28)
rownames(resultsAll) <- matrix(c(rep("cont",4), rep("copeak", 4), rep("disc", 4), rep("gaussian", 4), 
                                 rep("oscil", 4), rep("prpeak", 4), rep("step", 4)), ncol = 1)
colnames(resultsAll) <- c("Methods", "1", "2", "3", "5", "10", "20")   
resultsAll[, 1] <- c("True", "BART", "MI", "GP")  

mapeValues <- matrix(NA, ncol = 7, nrow = 21)
rownames(mapeValues) <- matrix(c(rep("cont",3), rep("copeak", 3), rep("disc", 3), rep("gaussian", 3), 
                                 rep("oscil", 3), rep("prpeak", 3),  rep("step", 3)), ncol = 1)
colnames(mapeValues) <- c("Methods", "1", "2", "3", "5", "10", "20") 
mapeValues[, 1] <- c("BART", "MI", "GP") 

# RMSE function
RMSE <- function(x, y){
  return( sqrt(mean((x - y)^2)) )
}


for (i in c(7)){
  for (j in c(1)){
    meanabsMape <- 0
    resultsAllEntry <- 0
    
    # global parameters: dimension
    dimensionsList <- c(1,2,3,5,10,20)
    dim <- dimensionsList[j]
    whichGenz <- i
    
    # Skip if dim = 1 or 2 for discontinuous
    # if (whichGenz == 3 & dim == 1 ) { 
    #   next 
    # } else if (whichGenz == 3 & dim == 2) {
    #   next
    # } else if (whichGenz == 3 & dim == 3) {
    #   next
    # }
    # Find Genz function
    if (whichGenz == 1) { genzFunctionName <-  "cont" }
    if (whichGenz == 2) { genzFunctionName <-  "copeak" }
    if (whichGenz == 3) { genzFunctionName <-  "disc" }
    if (whichGenz == 4) { genzFunctionName <-  "gaussian" }
    if (whichGenz == 5) { genzFunctionName <-  "oscil" }
    if (whichGenz == 6) { genzFunctionName <-  "prpeak" }
    if (whichGenz == 7) { genzFunctionName <-  "step" }
    
    for (num_cv in 1:20) {
      # Set path for estimated integral values
      fileName <- paste(toString(genzFunctionName), 'Dim', toString(dim), "", "Gaussian", "_", toString(num_cv),  '.csv', sep='')
      # fileName <- paste(toString(genzFunctionName), 'Dim', toString(dim), "Uniform", "_", toString(num_cv),  '.csv', sep='')
      filePath <- paste('results/genz', toString(whichGenz), fileName, sep='/')
      
      # Retrieve estimated integral values
      integrals <- read.csv(filePath, header=TRUE, sep=",", stringsAsFactors = FALSE)
      predictionBART <- data.frame("meanValueBART" = integrals[, 2], "standardDeviationBART" = integrals[, 3])
      predictionMonteCarlo <- data.frame("meanValueMonteCarlo" = integrals[, 4], "standardDeviationMonteCarlo" = integrals[, 5])
      predictionGPBQ <-  data.frame("meanValueGP" = integrals[, 6], "standardDeviationGP" = integrals[, 7])
      
      # Retrieve analytical integral values
      whichDimension <- which(dim == dimensionsList)
      if (whichGenz == 7) {
        real <- 0.5
      } else {
        analyticalIntegrals <- read.csv("results/genz/integrals.csv", header = FALSE)
        real <- analyticalIntegrals[whichGenz, whichDimension]
      }
      
      
      # Compute abs error
      epoch <- nrow(integrals)
      absMape <- abs((c(predictionBART$meanValueBART[epoch], 
                             predictionMonteCarlo$meanValueMonteCarlo[epoch], 
                             predictionGPBQ$meanValueGP[epoch]) - real) / real)
      meanabsMape <- meanabsMape + absMape 
      # Store all estimates in a single csv
      resultsAllEntry <- c(
            real, 
            predictionBART$meanValueBART[epoch], 
            predictionMonteCarlo$meanValueMonteCarlo[epoch], 
            predictionGPBQ$meanValueGP[epoch]
          ) + resultsAllEntry
      
      # Store all RMSE
      # rmse <- c(RMSE(predictionBART$meanValueBART[epoch], real),
      #           RMSE(predictionMonteCarlo$meanValueMonteCarlo[epoch], real),
      #           RMSE(predictionGPBQ$meanValueGP[epoch], real)
      # )
      # mapeValues[( (3*i - 2) : (3*i) ), (j + 1)] <- signif(rmse, 3)
    }
    meanabsMape <- meanabsMape / 20
    bestMethod[i, j] <- methods[which(meanabsMape == min(meanabsMape))[1]]
    mapeValues[( (3*i - 2) : (3*i) ), (j + 1)] <- signif(meanabsMape, 3)
    resultsAll[( (4*i - 3) : (4*i) ), (j + 1)] <- signif(resultsAllEntry / 20, 3)
  }
  cat(genzFunctionName, "done", '\n')
}

print(bestMethod)
print(mapeValues)
write.csv(bestMethod, file = "results/genz/bestMethods.csv")
write.csv(resultsAll, file = "results/genz/allEstimates.csv")
write.csv(mapeValues, file = "results/genz/mapeValues.csv")
