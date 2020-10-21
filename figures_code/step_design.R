# load("results/genz/7/stepDim1Uniform_5.RData")
load("results/genz/7/stepDim1Gaussian_5.RData")

nInit <- 20
nSeq <- 20
colInit <- rgb(0, 0, 0, 0.7)
orangered <- rgb(1, 0.271, 0, 0.5)
doggerblue <- rgb(0.118, 0.565, 1, 0.5)
crossWidth <- 3

trainX <- results_models$BART$trainData$trainX[1:nInit]
trainY <- results_models$BART$trainData$trainY[1:nInit]
newX <- results_models$BART$trainData$trainX[(nInit + 1):(nInit + nSeq)]
newY <- results_models$BART$trainData$trainY[(nInit + 1):(nInit + nSeq)]
gp_newX <- results_models$GP$X[(nInit + 1):(nInit + nSeq)]
gp_newY <- results_models$GP$Y[(nInit + 1):(nInit + nSeq)]

dim <- 1
source("src/genz/genz.R")
plotX <- as.matrix(seq(0, 1, 1/1000))

pdf(paste("Figures/genz", "/step_design.pdf", sep = ""), width = 9, height = 5)
par(mfrow=c(1,3), pty="s")
plot(plotX, step(plotX, jumps=1), ty="l", ylim = c(-0.5, 1.5), xlab="x", ylab = "y",
     cex.axis = 2, cex.lab = 2, cex.main=2)
points(trainX[order(trainX)], trainY[order(trainX)], col=colInit, bg=colInit, pch=4, lwd=crossWidth, cex=1.5)
legend("topleft", legend="Design Points",
       col=colInit, cex=1.6, pch=4, bty="n", pt.lwd=crossWidth)

plot(plotX, step(plotX, jumps=1), ty="l", ylim=c(-0.5, 1.5), xlab="x", ylab = "y",
     cex.axis = 2, cex.lab = 2, cex.main=2)
points(trainX[order(trainX)], trainY[order(trainX)], col=colInit, bg=colInit, pch=4, lwd=crossWidth, cex=1.5)
points(newX[order(newX)], newY[order(newX)], col=orangered, bg=orangered, pch=21, cex=2)
legend("topleft", legend=c("Design Points", "BART-Int"),
       col=c(colInit, orangered), cex=1.6, pch=c(4, 19), bty="n", pt.lwd=c(crossWidth, 0))

plot(plotX, step(plotX, jumps=1), ty="l", ylim=c(-0.5, 1.5), xlab="x", ylab = "y",
     cex.axis = 2, cex.lab = 2, cex.main = 2)
points(trainX[order(trainX)], trainY[order(trainX)], col=colInit, bg=colInit, pch=4, lwd=crossWidth, cex=1.5)
points(gp_newX[order(gp_newX)], gp_newY[order(gp_newX)], col=doggerblue, bg=doggerblue, pch=21, cex=2)
legend("topleft", legend=c("Design Points", "GP-BQ"),
       col=c(colInit, doggerblue), cex=1.6, pch=c(4, 19), bty="n", pt.lwd=c(crossWidth, 0))
dev.off()
