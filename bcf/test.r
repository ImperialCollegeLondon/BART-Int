library(bcf)
# Uncomment the following for the modified bcf package
# source("bcf/R/RcppExports.R")
# source("bcf/R/bcf.R")

args <- commandArgs(TRUE)
ntrain <- as.double(args[1])
ncandidate <- as.double(args[2])
n_seqential <- as.double(args[3])
sigma <- as.double(args[4]) # observation noise of y
homogeneous <- as.double(args[5]) # 1 if true
linear <- as.character(args[6]) # 1 if linear
num_cv <- as.character(args[7]) # 1 if linear
which_experiment <- as.double(args[8]) # 1 or 2

# ntrain <- 250
# ncandidate <- 2000
# sigma <- 0.1
# homogeneous <- 1
# linear <- 1
# n_seqential <- 1
# which_experiment <- 2
n <- ntrain + ncandidate
set.seed(1)
p <- 5

# Experiment
if (which_experiment == 1) {
  source("bcf/util_experiment1.R")

  # Create covariates x
  x <- generate_x(n, p)
  x_input <- dbarts::makeModelMatrixFromDataFrame(x)

  # Compute tau, mu and pi
  mu <- prognostic_func(x, linear)
  tau <- treatment_func(x, homogeneous)
  propensity <- 0.8 * pnorm(3 * mu / sd(mu) - 0.5 * x[, 1]) + 0.05 + 0.1 * runif(n)
  z <- rbinom(n, rep(1, n), propensity)

} else if (which_experiment) {
  source("bcf/util_experiment2.R")

  # Create covariates x
  x <- generate_x(n, p)
  x_input <- dbarts::makeModelMatrixFromDataFrame(x)

  # Compute tau, mu and pi
  # RIC
  mu <- prognostic_func(x, linear)
  tau <- treatment_func(x, mu)
  propensity <- 0.6 * pnorm(mu - mean(mu), 0, 0.75 * sd(mu)) + 0.2
  z <- rbinom(n, rep(1, n), propensity)

}


# Model: 
# Y_i = f(x_i, Z_i) + epsilon_i
# epsilon_i ~ N(0, sigma^2)
# f(x_i, z_i) = mu(x_i) + tau(x_i) * z_i
Ey <- mu + tau * z
noiseVar <- sigma * sd(Ey)
y <- Ey + noiseVar*rnorm(n)

# If pi is unknown, we would need to estimate it here
# pihat <- propensity
library(nnet)
x.mod <- dbarts::makeModelMatrixFromDataFrame(data.frame(x))
fitz <- nnet(z~., data = x.mod, size = 4, rang = 0.1, maxit = 1000, abstol = 1.0e-8, decay = 5e-2, trace=FALSE)
pihat <- fitz$fitted.values

# BCF
# bcf_fit <- bcf(y[1:ntrain], z[1:ntrain], x_input[1:ntrain, ], x_input[1:ntrain, ], pihat[1:ntrain], nburn=2000, nsim=2000, nthin= 5)
bcf_fit <- bcf(y[1:ntrain], z[1:ntrain], x_input[1:ntrain, ], x_input[1:ntrain, ], pihat[1:ntrain], 10000, 3000)

# Get posterior samples of treatment effects
tau_post <- bcf_fit$tau
tauhat <- colMeans(tau_post)
# plot(tau[1:ntrain], tauhat); abline(0,1)

# Posterior of the averaged treatment effects
print(paste0("Mean of |CATE - CATE_hat|: ", mean(abs(tau[1:ntrain] - tauhat))))

if (homogeneous == 1 & which_experiment == 1){
  # Sample ATE
  y_treated <- bcf_fit$yhat
  y_controled <- bcf_fit$yhat
  y_treated[z == 0] <- (y_treated + tauhat)[z == 0]
  y_controled[z == 1] <- (y_controled - tauhat)[z == 1]
  sate <- mean(y_treated - y_controled)
  true_sate <- mean(tau[1:ntrain])
  print(paste0("True ATE: ", true_sate, " . Estimated: ", sate, ". Abs. diff: ", abs(sate - 3)))
} else {
  y_treated <- bcf_fit$yhat
  y_controled <- bcf_fit$yhat
  y_treated[z == 0] <- (y_treated + tauhat)[z == 0]
  y_controled[z == 1] <- (y_controled - tauhat)[z == 1]
  sate <- mean(y_treated - y_controled)
  true_sate <- mean(tau[1:ntrain])
  print(paste0("True SATE: ", true_sate, ". Estimated: ", sate, ". Abs. diff: ", abs(sate - true_sate))) 
}

# Gaussian processes
# makeGPBQModelMatrix <- function(df, treatment) {
#   return(data.frame(df, treatment))
# }
makeGPBQModelMatrix <- function(df, treatment) {
  return(data.frame(df, treatment))
}
# Uncomment the following to exclude treatment z in the data input 
# makeGPBQModelMatrix <- function(df) {
#   return(dbarts::makeModelMatrixFromDataFrame(df))
# }


# Uncomment the following to include treatment z in the data input 
trainX <- makeGPBQModelMatrix(x[1:ntrain,], z[1:ntrain])
trainY <- y[1:ntrain]
candidateX <- makeGPBQModelMatrix(x[-(1:ntrain),], z[-(1:ntrain)])
candidateY <- y[-(1:ntrain)]
# Uncomment the following to exclude treatment z in the data input 
# trainX <- makeGPBQModelMatrix(x[1:ntrain,])
# trainY <- y[1:ntrain]
# candidateX <- makeGPBQModelMatrix(x[-(1:ntrain),])
# candidateY <- y[-(1:ntrain)]

library(reticulate)
source("src/optimise_gp.R")
lengthscale <- optimise_gp_r(dbarts::makeModelMatrixFromDataFrame(trainX), trainY, kernel = "matern32", epochs = 500)
print("...Finished training for the lengthscale")

source("bcf/GPBQ.R")
GPBQResults <- computeGPBQEmpirical(
  X=trainX, 
  Y=trainY, 
  candidateX=candidateX, 
  candidateY=candidateY, 
  kernel="matern32", 
  lengthscale=lengthscale, 
  epochs=n_seqential
)
 
# BART
source("bcf/BART.R")
BARTResults <- computeBARTWeighted(trainX, trainY, candidateX, candidateY, n_seqential, num_cv=num_cv, linear=linear)
print(paste0("True ATE: ", true_sate, " . Estimated: ", BARTResults$meanValueBART, ". Abs. diff: ", abs(BARTResults$meanValueBART - sate)))

# Bayesian Quadrature methods: with BART and Gaussian Process
print("Final ATE:")
print(c("Actual integral:", true_sate))
print(c("BCF integral:", sate))
print(c("BART integral:", BARTResults$meanValueBART[(n_seqential+1)]))
print(c("GP integral:", GPBQResults$meanValueGP[(n_seqential+1)]))

results <- data.frame(
  "epochs" = c(1:(n_seqential+1)),
  "BARTMean" = BARTResults$meanValueBART, "BARTsd" = BARTResults$standardDeviationBART,
  "GPMean" = GPBQResults$meanValueGP, "GPsd" = sqrt(GPBQResults$varianceGP),
  "actual" = rep(3, n_seqential+1)
)

write.csv(results, paste0("bcf/results/homo", homogeneous, "_", "sigma", sigma, "_", ntrain, "_", num_cv, ".csv"))
