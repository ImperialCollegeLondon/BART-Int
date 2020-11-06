for (num_cv in 3:3) {
  pdf(paste("bcf/plot1", ".pdf", sep = ""), width = 5, height = 5)
  results <- read.csv(paste("bcf/results/results/exp1_sigma0.1_ntrain50_nseq200_homo1_sequential1_", num_cv, ".csv", sep=""))
  for (k in 0:20) {
    n_seq <- k*10
    if (k==0) {
      bart_posterior <- load(paste("bcf/results/results/posterior_BART_synthetic_sequential1_", 1, "_", num_cv, ".RData", sep=""))
      plot(rep(50, length(posterior_samples$posterior_samples)), posterior_samples$posterior_samples,
           col="orangered", ylim=c(2.2,4), ylab = "ATE", xlim=c(50, 250), 
           xlab = expression(italic(n)[ini] + italic(n)[seq]), cex.lab = 1.5, cex.axis = 1.5, cex=0.1)
      points(50, results$BARTMean[1], col = "orangered", bg='orangered', pch=21, lwd=3)
      
    } else {
      bart_posterior <- load(paste("bcf/results/results/posterior_BART_synthetic_sequential1_", n_seq, "_", num_cv, ".RData", sep=""))
      points(n_seq+50, results$BARTMean[n_seq], col = "orangered", bg='orangered', pch=21, lwd=3)
      points(rep(n_seq+50, length(posterior_samples$posterior_samples)), 
             posterior_samples$posterior_samples, 
             col = "orangered", bg='orangered', cex=0.1)
    }
  }
  abline(h=3)
   
  # random
  results <- read.csv(paste("bcf/results/results/exp1_sigma0.1_ntrain50_nseq200_homo1_sequential2_", num_cv, ".csv", sep=""))
  for (k in 0:19) {
    n_seq <- k*10 + 5
    if (k==0) {
      bart_posterior <- load(paste("bcf/results/results/posterior_BART_synthetic_sequential2_", 2, "_", num_cv, ".RData", sep=""))
      points(rep(50, length(posterior_samples$posterior_samples)), posterior_samples$posterior_samples,
           col="chartreuse4",  cex=0.1)
      points(50, results$BARTMean[2], col = "chartreuse4", bg='chartreuse4', pch=21, lwd=3)
      
    } else {
      bart_posterior <- load(paste("bcf/results/results/posterior_BART_synthetic_sequential2_", n_seq, "_", num_cv, ".RData", sep=""))
      points(n_seq+50, results$BARTMean[n_seq], col = "chartreuse4", bg='chartreuse4', pch=21, lwd=3)
      
      points(rep(n_seq+50, length(posterior_samples$posterior_samples)), 
             posterior_samples$posterior_samples, 
             col = "chartreuse4", bg='chartreuse4', cex=0.1)
    }
  }
  legend("topright", legend=c("BART Active", "BART Random"),
         col=c("orangered", "chartreuse4"), cex=1.3, lty = c(1,1), bty="n")
  dev.off()
  
  pdf(paste("bcf/plot2", ".pdf", sep = ""), width = 5, height = 5)
  results <- read.csv(paste("bcf/results/results/exp1_sigma0.1_ntrain50_nseq200_homo1_sequential1_", num_cv, ".csv", sep=""))
  plot(results$epochs+50, results$GPMean, ty="l", col="dodgerblue", 
       ylim=c(2.2,3.5), ylab = "ATE", xlab = expression(italic(n)[ini] + italic(n)[seq]), cex.lab = 1.5, cex.axis = 1.5)
  polygon(c(results$epochs+50, rev(results$epochs+50)), 
          c(
            results$GPMean + 2*results$GPsd, 
            rev(results$GPMean - 2*results$GPsd)
          ), 
          col = adjustcolor("dodgerblue", alpha.f = 0.10), 
          border = "dodgerblue", lty = c("dashed", "solid"))
  abline(h=3)
  results <- read.csv(paste("bcf/results/results/exp1_sigma0.1_ntrain50_nseq200_homo1_sequential2_", num_cv, ".csv", sep=""))
  points(results$epochs+50, results$GPMean, ty="l", col="brown", ylim=c(2,4))
  polygon(c(results$epochs+50, rev(results$epochs+50)), 
          c(
            results$GPMean + 2*results$GPsd, 
            rev(results$GPMean - 2*results$GPsd)
          ), 
          col = adjustcolor("brown", alpha.f = 0.10), 
          border = "brown", lty = c("dashed", "solid"))  
  legend("topright", legend=c("GP Active", "GP Random"),
         col=c("dodgerblue", "brown"), cex=1.3, lty = c(1,1), bty="n")
  dev.off()
}