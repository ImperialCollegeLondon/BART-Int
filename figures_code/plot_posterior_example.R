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
dim <- 1
num_iterations <- 1
whichGenz <- 7
whichKernel <- "matern32"
jumps = 1
# turn on/off sequential design
# 1 denotes TRUE to sequential
# 0 denotes FALSE to sequential
cat("\nBegin testing:\n")
sequential <- TRUE
cat("Sequantial design set to", sequential, "\n")
# prior measure over the inputs
# uniform by default
measure <- "uniform"
cat("Prior measure:", measure, "\n")
# extra parameter for step function
print(c(dim, num_iterations, whichGenz))
source("src/genz/genz.R") # genz function to test

if (whichGenz < 1 | whichGenz > 8) { stop("undefined genz function. Change 3rd argument to 1-8") }
if (whichGenz == 3 & dim == 1) { stop("incorrect dimension. Discrete Genz function only defined for dimension >= 2") }

if (whichGenz == 1) { genz <- cont; genzFunctionName <- deparse(substitute(cont)) }
if (whichGenz == 2) { genz <- copeak; genzFunctionName <- deparse(substitute(copeak)) }
if (whichGenz == 3) { genz <- disc; genzFunctionName <- deparse(substitute(disc)) }
if (whichGenz == 4) { genz <- gaussian; genzFunctionName <- deparse(substitute(gaussian)) }
if (whichGenz == 5) { genz <- oscil; genzFunctionName <- deparse(substitute(oscil)) }
if (whichGenz == 6) { genz <- prpeak; genzFunctionName <- deparse(substitute(prpeak)) }
if (whichGenz == 7) { genz <- function(xx) { return(step(xx, jumps = jumps)) }; genzFunctionName <- deparse(substitute(step)) }
if (whichGenz == 8) { genz <- mix; genzFunctionName <- deparse(substitute(mix)) }

print("Testing with: %s" %--% genzFunctionName)

# prepare training dataset
if (measure == "uniform") {
  trainX <- replicate(dim, runif(20))
} else if (measure == "gaussian") {
  trainX <- replicate(dim, rtnorm(20, lower = 0, upper = 1))
}
trainY <- genz(trainX)
source("src/BARTBQ.R")
t0 <- proc.time()
posterior_model <- BART_posterior(dim, trainX, trainY, num_iterations, FUN = genz, sequential, measure)
t1 <- proc.time()
bartTime <- (t1 - t0)[[1]]

x_plot <- replicate(dim, runif(500))
x_plot <- x_plot[order(x_plot),]
y_pred <- predict(posterior_model$model, x_plot)
y_pred_mean <- colMeans(y_pred)
y_pred_sd <- sqrt(colVars(y_pred))

# obtain posterior samples
integrals <- sampleIntegrals(posterior_model$model, dim, measure)
ymin <- min(posterior_model$trainData[, (dim + 1)]);
ymax <- max(posterior_model$trainData[, (dim + 1)])
integrals <- (integrals + 0.5) * (ymax - ymin) + ymin

plot(x_plot, y_pred_mean)
hist(integrals)
# 
# if (!sequential){
#   figName <- "Figures/%s/drawBART%s%sDimNoSequential.pdf" %--% c(whichGenz, genzFunctionName, dim)
#   csvName <- "Figures/%s/drawBART%s%sDimNoSequential.csv" %--% c(whichGenz, genzFunctionName, dim)
#   groundTruthName <- "Figures/%s/trainDrawBart%s%sDimNoSequential.csv" %--% c(whichGenz, genzFunctionName, dim)
# } else {
#   figName <- "Figures/%s/drawBART%s%sDim.pdf" %--% c(whichGenz, genzFunctionName, dim)
#   csvName <- "Figures/%s/drawBART%s%sDim.csv" %--% c(whichGenz, genzFunctionName, dim)
#   groundTruthName <- "Figures/%s/trainDrawBart%s%sDim.csv" %--% c(whichGenz, genzFunctionName, dim)
# }

results <- data.frame(
  "x_plot" = x_plot,
  "y_pred" = y_pred_mean
)
groundTruth <- data.frame(
  "trainX" = trainX,
  "trainY" = trainY
)
# write.csv(results, file = csvName, row.names=FALSE)
# write.csv(groundTruth, file = groundTruthName, row.names=FALSE)

# Bayesian Quadrature with Gaussian Process
print("Begin Gaussian Process Integration")
library(reticulate)
source("src/optimise_gp.R")
lengthscale <- optimise_gp_r(trainX, trainY, kernel = "matern32", epochs = 500)

source("src/GPBQ.R")
t0 <- proc.time()
# need to add in function to optimise the hyperparameters
predictionGPBQ <- computeGPBQ_matern(trainX, trainY, dim, epochs = num_iterations, kernel = "matern32", FUN = genz, lengthscale, sequential, measure)
t1 <- proc.time()
GPTime <- (t1 - t0)[[1]]

results_models <- list("BART" = posterior_model, "GP" = predictionGPBQ, "trainX" = trainX, "trainY" = trainY, "lengthscale" = lengthscale)
# save(results_models, file = "plot_posterior_example.RData")

K <- predictionGPBQ$K
X <- predictionGPBQ$X
# Y <- Y[order(X)]
Y <- predictionGPBQ$Y
maternKernel <- maternKernelWrapper_2(lengthscale = lengthscale)

k_xstar_x <- kernelMatrix(maternKernel, matrix(x_plot, ncol = 1), X)
k_xstar_xstar <- kernelMatrix(maternKernel,
                              matrix(x_plot, ncol = 1),
                              matrix(x_plot, ncol = 1))
jitter = 1e-6
K_inv <- solve(K + diag(jitter, nrow(K)))
gp_post_mean <- k_xstar_x %*% K_inv %*% Y
gp_post_cov <- k_xstar_xstar - k_xstar_x %*% K_inv %*% t(k_xstar_x)
gp_post_sd <- sqrt(diag(gp_post_cov))

#plot of integrals
GPdensity <- dnorm(
  seq(0, 1, 0.01),
  mean = predictionGPBQ$meanValueGP[1],
  sd = sqrt(predictionGPBQ$varianceGP[1])
)
KDE_BART <- density(integrals)

pdf("figures_code/posterior_step.pdf", width=9, height=3)
par(mfrow = c(1,3), pty = "s")
plot(
  seq(0, 1, 0.01), 
  GPdensity, 
  ty="l", 
  col = "dodgerblue", 
  xlim = c(0.3,0.7), 
  ylim = c(0, 60),
  xlab = "x",
  ylab = "Posterior density",
  cex.lab = 1.7,
  cex.axis = 1.7,
  lwd=3
)
points(KDE_BART, ty="l", col = "orangered", lwd=3)
abline(v=0.5, lwd=3)
legend("topleft", legend=c("BART", "GP", expression(Pi*"[f]")),
       col=c("orangered", "dodgerblue", "black"), lty = c(1,1,1), cex=2.4, bty="n")

plot(x_plot, y_pred_mean, col = "orangered", ty="l", lwd=3,ylim=c(-0.2, 1.5),
     ylab = expression(f(x)), xlab = "x",      cex.lab = 1.7,
     cex.axis = 1.7)
# points(trainX, trainY, col = "black")
points(c(0,0.5,0.5, 1), c(0,0,1,1), ty="l", lwd=1, col="black")

for (i in seq(1, 500, 10)) {
  points(rep(x_plot[i], 1000), y_pred[,i], col="orangered",bg='orangered',
         cex=0.2, alpha=0.01, pch=16)
}
# legend("topright", legend=c("BART-Int", "GP-BQ"),
#        col=c("orangered", "dodgerblue"), cex=1, lty = c(1,1,1,1))


a <-density(integrals)$y 
# plot(trainX, trainY, ylim=c(-0.2, 1.3))
plot(x_plot, 
     gp_post_mean, 
     col = "dodgerblue", 
     cex=0.2, 
     ty="l", 
     ylim=c(-0.2, 1.5),
     xlab = "x",
     ylab = expression(f(x)),
     cex.lab = 1.7,
     cex.axis = 1.7,
     lwd=3
)


points(c(0,0.5,0.5, 1), c(0,0,1,1), ty="l", lwd=1, col="black")
# points(trainX[order(trainX),], trainY[order(trainX), ], col = "black", bg='black', pch=21, lwd=3)
polygon(c(x_plot, rev(x_plot)), 
        c(
          gp_post_mean + 1.96*gp_post_sd, 
          rev(gp_post_mean - 1.96*gp_post_sd)
        ), 
        col = adjustcolor("dodgerblue", alpha.f = 0.10), 
        border = "dodgerblue", lty = c("dashed", "solid"))
# points(trainX, trainY, col = "black")

dev.off()
