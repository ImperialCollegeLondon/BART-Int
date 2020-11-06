monteCarloIntegrationUniform <- function(FUN, trainX, trainY, numSamples, dim, measure)
  #'Crude Monte Carlo Approximation
  #' 
  #'@description The function approximates the integral of interest using curde monte carlo
  #'There is no sequential sampling/adaptive Bayesian quadrature as there are no posterior samples
  #' 
  #'@param FUN Function; The function to be integrated 
  #'@param numSamples Integer; The number of samples used in calculate mean
  #'@param dim Integer; The dimension of the input X
  #'
  #'@return List; A list containing meanValue (appximation) and the variance of crude monte Carlo
{
	meanValueMonteCarlo <- rep(0, numSamples)
  	standardDeviationMonteCarlo <- rep(0, numSamples)
  	trainData <- data.frame(trainX, trainY)
	
	meanValueMonteCarlo[1] <- mean(trainData[, dim+1])
	standardDeviationMonteCarlo[1] <- sqrt(var(trainData[, dim+1]))
	if (numSamples == 1) {
		return(list("meanValueMonteCarlo" = meanValueMonteCarlo, 
					"standardDeviationMonteCarlo" = standardDeviationMonteCarlo, 
					"trainData"=trainData))
	}

	for (i in 2:numSamples) {
		
		if (measure == "uniform") {
			candidateX <- matrix(replicate(dim, runif(1, 0, 1)), ncol=dim)
		} else if (measure == "gaussian") {
			candidateX <- matrix(replicate(dim, rtnorm(1, mean=0.5, lower=0, upper=1)), ncol=dim)
		} else if (measure == "exponential") {
		  candidateX <- matrix(replicate(dim, rexp(1)), ncol=dim)
		}
		candidateY <- FUN(candidateX)
		trainData <- rbind(trainData, c(candidateX, candidateY))
		meanValueMonteCarlo[i] <- mean(trainData[, dim+1])
		standardDeviationMonteCarlo[i] <- sqrt(var(trainData[,dim+1]))
	}
	return(list("meanValueMonteCarlo" = meanValueMonteCarlo, 
				"standardDeviationMonteCarlo" = standardDeviationMonteCarlo, 
				"trainData"=trainData))
}
