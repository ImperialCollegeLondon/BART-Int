# define string formatting
`%--%` <- function(x, y) 
  # from stack exchange:
  # https://stackoverflow.com/questions/46085274/is-there-a-string-formatting-operator-in-r-similar-to-pythons
{
  do.call(sprintf, c(list(x), y))
}
bart <- c()
mi <- c()
gp <- c()
# generate the plots for all 20 runs
for (num_cv in 1:20) {

  # read in results outputted
  df <- read.csv(paste("results/survey_design/results", num_cv, ".csv", sep=""))

  # estimated by running BART on all the labels
  real <- 0.5804115

  # compute MAPE 
  bart <- c(bart, abs((df$BARTMean[200] - real)/ real))
  gp <- c(gp, abs((df$GPMean[200] - real) / real))
  mi <- c(mi, abs((df$MIMean[200] - real) / real))
  
  # create some simple jpg plots
  jpeg(paste("Figures/survey_design/", num_cv, ".jpg"), width=500, height=500)
  plot(df$epochs, df$BARTMean,ty="l", ylim = c(0.4, 0.8), col="orangered")
  points(df$epochs, df$MIMean, ty="l", col="chartreuse4")
  points(df$epochs, df$GPMean, ty="l", col="dodgerblue")
  abline(h=real)
  dev.off()
}

# final mape results
bart_mape = signif(mean(bart), 7)
bart_sd = signif(sd(bart) / sqrt(20), 7)
gp_mape = signif(mean(abs(gp)), 7)
gp_sd = signif(sd(gp) / sqrt(20), 7)
mi_mape = signif(mean(abs(mi)), 7)
mi_sd = signif(sd(mi)/sqrt(20), 7)
print(c("BART MAPE and Std Error:", bart_mape, bart_sd))
print(c("GP MAPE and Std Error:", gp_mape, gp_sd))
print(c("MI MAPE and Std Error:", mi_mape, mi_sd))

# we used 20 design points
num_design=20
## make better looking plots
for (num_cv in c(1:20)) {
  results <- read.csv(paste("results/survey_design/results", num_cv, ".csv", sep=""))
  real <- 0.5804115
  
  # 1. Open jpeg file
  pdf(paste0("Figures/survey_design", "resultsBin", num_cv, ".pdf"), width = 9, height = 5)
  # jpeg(paste0(plotPath, "results", num_cv, ".jpg"), width = 1000, height = 1000)
  par(mfrow = c(1,2), pty = "s")
  ymax <- max(c(abs(results$BARTMean - real), abs(results$GPMean - real)))
  plot(results$epochs+num_design,
       abs(results$BARTMean - real),
       ty="l",
       xlab = expression(n[seq]),
       ylab = "Absolute Error",
       col = "orangered",
       ylim = c(0.001, 1),
       log="y",
       cex.lab=1.6,
       cex.axis=2,
  )
  
  points(results$epochs+num_design, abs(results$MIMean - real), ty="l", col = "chartreuse4")
  points(results$epochs+num_design, abs(results$GPMean - real), ty="l", col = "dodgerblue")
  
  ymin <- min(c(results$BARTMean - 2*results$BARTsd, results$GPMean - 2*sqrt(results$GPsd)))
  ymax <- max(c(results$BARTMean + 2*results$BARTsd, results$GPMean + 2*sqrt(results$GPsd)))    
  legend("topleft", legend=c("Monte Carlo", "BART-Int", "GPBQ"),
         col=c("chartreuse4", "orangered", "dodgerblue"), cex=1.6, lty = c(1,1,1), bty="n")
  plot(results$epochs+num_design, 
       results$BARTMean, 
       ty="l", 
       ylab = "Proportion Log(income) > 10",
       xlab = expression("n'"+n[seq]),
       col = "orangered", 
       lwd=3,
       ylim = c(0, 1),
       cex.lab=1.5,
       cex.axis=2,
  )
  
  # now read in the posterior samples of BART-Int at each iteration 
  for (n_seq in seq(1, 200, 10)) {
    bart_posterior <- load("results/survey_design/posterior_BART_survey_%s_%s.RData" %--% c(n_seq, num_cv))
    if (n_seq == 1) {
      num_posterior_samples <- length(posterior_samples$posterior_samples)
      points(rep(n_seq+num_design, num_posterior_samples), 
             posterior_samples$posterior_samples,
             col = "orangered",bg='orangered',
             cex=0.2, pch=19, alpha=0.01)
    } else {
      points(rep(n_seq+num_design, num_posterior_samples), posterior_samples$posterior_samples, 
             col = "orangered", bg='orangered', cex=0.2, pch=19, alpha=0.01)
    }
  }
  points(results$epochs+num_design, results$MIMean, ty="l", col = "chartreuse4", lwd=3)
  points(results$epochs+num_design, results$GPMean, ty="l", col = "dodgerblue", lwd=3)
  polygon(c(results$epochs+num_design, rev(results$epochs+num_design)), 
          c(
            results$GPMean + 1.96*sqrt(results$GPsd), 
            rev(results$GPMean - 1.96*sqrt(results$GPsd))
          ), 
          col = adjustcolor("dodgerblue", alpha.f = 0.10), 
          border = "dodgerblue", lty = c("dashed", "solid"))
  abline(h=real)
  dev.off()
}
