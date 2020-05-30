load("results/genz/7/stepDim1Uniform_5.RData")

trainX <- results_models$BART$trainData$trainX[1:50]
trainY <- results_models$BART$trainData$trainY[1:50]
newX <- results_models$BART$trainData$trainX[51:70]
newY <- results_models$BART$trainData$trainY[51:70]
gp_newX <- results_models$GP$X[51:70]
gp_newY <- results_models$GP$Y[51:70]

dim <- 1
source("src/genz/genz.R")
plotX <- as.matrix(seq(0, 1, 1/1000))

pdf(paste("Figures/genz", "/step_design.pdf", sep = ""), width = 9, height = 5)
par(mfrow=c(1,2), pty="s")
plot(plotX, step(plotX, jumps=1), ty="l", ylim = c(-0.5, 1.5), xlab="x", ylab = "y",
     cex.axis = 2, cex.lab = 2, main="BART-Int", cex.main=2)
points(trainX[order(trainX)], trainY[order(trainX)], col = "black", bg='black', pch=21, lwd=3)
points(newX[order(newX)], newY[order(newX)], col = "orangered", bg='orangered', pch=21, lwd=3)
legend("topleft", legend=c("Design Points", "BART-Int"),
       col=c("black", "orangered"), cex=1.6, pch = c(19, 19), bty="n")

plot(plotX, step(plotX, jumps=1), ty="l", ylim = c(-0.5, 1.5), xlab="x", ylab = "y",
     cex.axis = 2, cex.lab = 2, main="GP-BQ", cex.main = 2)
points(trainX[order(trainX)], trainY[order(trainX)], col = "black", bg='black', pch=21, lwd=3)
points(gp_newX[order(gp_newX)], gp_newY[order(gp_newX)], col = "dodgerblue", bg='dodgerblue', pch=21, lwd=3)
legend("topleft", legend=c("Design Points", "GP-BQ"),
       col=c("black", "dodgerblue"), cex=1.6, pch = c(19, 19), bty="n")
dev.off()
