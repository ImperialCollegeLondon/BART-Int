# Plot results for each Genz integral after 500 epochs.
# To run:
#     In Unix terminal: 
#     Rscript {ROOT}/src/genz/drawGraphs.R
# To manually adjust y limits of the plots:
#     Uncomment ylim = ylims[1:2] and ylim = ylims[3:3] in the two plot functions, 
#     and input the four limits when running the file
# To see results:
#     Plots stored in {ROOT}/report/Figures for plots
# define string formatting
`%--%` <- function(x, y) 
  # from stack exchange:
  # https://stackoverflow.com/questions/46085274/is-there-a-string-formatting-operator-in-r-similar-to-pythons
{
  do.call(sprintf, c(list(x), y))
}
ylims <- as.double(commandArgs(TRUE))
mape_error <- read.csv("results/genz/mapeValues.csv")
num_cv_plot=18
for (i in c(7)){
  for (j in c(1)){
    for (sequential in c("")){
      # global parameters: dimension
      dimensionsList <- c(1,10,20)
      dim <- dimensionsList[j]
      whichGenz <- i
      
      # Skil if dim = 1 for discontinuous
      if (whichGenz == 3 & dim == 1) { break } 
      
      # Find Genz function
      if (whichGenz == 1) { genzFunctionName <-  "cont" }
      if (whichGenz == 2) { genzFunctionName <-  "copeak" }
      if (whichGenz == 3) { genzFunctionName <-  "disc" }
      if (whichGenz == 4) { genzFunctionName <-  "gaussian" }
      if (whichGenz == 5) { genzFunctionName <-  "oscil" }
      if (whichGenz == 6) { genzFunctionName <-  "prpeak" }
      if (whichGenz == 7) { genzFunctionName <-  "step" }
      if (whichGenz == 8) { genzFunctionName <-  "additive_gaussian" }
      
      for (num_cv in c(num_cv_plot)) {
        # Set path for estimated integral values
        fileName <- paste(toString(genzFunctionName), 'Dim', toString(dim), "Uniform_", toString(num_cv),'.csv', sep='')
        filePath <- paste('results/genz', toString(whichGenz), fileName, sep='/')
        
        # Retrieve estimated integral values
        integrals <- read.csv(filePath, header=TRUE, sep=",", stringsAsFactors = FALSE)
        predictionBART <- data.frame("meanValueBART" = integrals[, 2], "standardDeviationBART" = integrals[, 3])
        predictionMonteCarlo <- data.frame("meanValueMonteCarlo" = integrals[, 4], "standardDeviationMonteCarlo" = integrals[, 5])
        predictionGPBQ <-  data.frame("meanValueGP" = integrals[, 6], "standardDeviationGP" = integrals[, 7])
        
        # Retrieve analytical integral values
        whichDimension <- which(dim == dimensionsList)
        analyticalIntegrals <- read.csv("results/genz/integrals.csv", header = FALSE)
        if (whichGenz==7) {
          real = 0.5
        } else {
          real <- analyticalIntegrals[whichGenz, whichDimension]
        }
      }
    }
    plot_points <- seq(0, 50, 5)
    pdf(paste("Figures/genz", "/combined_", toString(whichGenz), ".pdf", sep = ""), width = 5, height = 5)
    par(pty = "s")
    plot(integrals$epochs+50, 
         integrals$MIMean, 
         ty="l", 
         ylab = "Integral",
         xlab = expression("n'"+n[seq]),
         col = "chartreuse4",
         ylim = c(0.45, 0.59),
         cex.lab = 1.7,
         cex.axis = 1.7
    )
    points(integrals$epochs[plot_points]+50, integrals$MIMean[plot_points], col = "chartreuse4", bg='chartreuse4', pch=21, lwd=3)
    points(integrals$epochs+50, integrals$GPMean, ty="l", col = "dodgerblue", lwd=3)
    points(integrals$epochs[plot_points]+50, integrals$GPMean[plot_points], col = "dodgerblue", bg='dodgerblue', pch=21, lwd=3)
    polygon(c(integrals$epochs+50, rev(integrals$epochs+50)), 
            c(
              integrals$GPMean + 2*integrals$GPsd, 
              rev(integrals$GPMean - 2*integrals$GPsd)
            ), 
            col = adjustcolor("dodgerblue", alpha.f = 0.10), 
            border = "dodgerblue", lty = c("dashed", "solid"))
    points(integrals$epochs+50, integrals$BARTMean, ty="l", col = "orangered", lwd=3)
    points(integrals$epochs[plot_points]+50, integrals$BARTMean[plot_points], col = "orangered",bg='orangered', pch=21, lwd=3)
    # polygon(c(integrals$epochs+50, rev(integrals$epochs+50)),
    #         c(
    #           integrals$BARTMean + 2*integrals$BARTsd,
    #           rev(integrals$BARTMean - 2*integrals$BARTsd)
    #         ),
    #         col = adjustcolor("orangered", alpha.f = 0.10),
    #         border = "orangered", lty = c("dashed", "solid"))
    for (n_seq in 1:20) {
      bart_posterior <- load("results/genz/posterior_BART_Dim1_step_%s_%s.RData" %--% c(num_cv_plot,n_seq))
      if (n_seq == 1) {
        num_posterior_samples <- length(posterior_samples$posterior_samples)
        points(rep(n_seq+50, num_posterior_samples), 
             posterior_samples$posterior_samples,
             col = "orangered",bg='orangered',
             cex=0.1)
      } else {
        points(rep(n_seq+50, num_posterior_samples), posterior_samples$posterior_samples, 
               col = "orangered", bg='orangered', cex=0.1)
      }
    }
    legend("topleft", legend=c("BART", "GP", "MI", expression(Pi*"[f]")),
           col=c("orangered", "dodgerblue", "chartreuse4", "black"), cex=1.4, lty = c(1,1,1,1), bty="n", ncol=2)
    abline(h=integrals$actual)
    dev.off()
  }
}