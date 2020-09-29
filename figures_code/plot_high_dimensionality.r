# define string formatting
`%--%` <- function(x, y) 
  # from stack exchange:
  # https://stackoverflow.com/questions/46085274/is-there-a-string-formatting-operator-in-r-similar-to-pythons
{
  do.call(sprintf, c(list(x), y))
}


results <- data.frame(matrix(0, ncol = 6, nrow = 4))
colnames(results) <- c("BART-Int", "BART-Int-se", "GPBQ", "GPBQ-se", "MI", "MI-se")
dims = c(1, 10, 20, 100)
l = 1
for (dim in dims) {
  bart = c()
  gp = c()
  mi = c()
  for (num_cv in 1:20) {
    df <- read.csv("results/genz/9/additive_gaussianDim%sNoSequentialUniform_%s.csv" %--% c(dim, num_cv))
    bart = c(bart, df$BARTMean)
    gp = c(gp, df$GPMean)
    mi = c(mi, df$MIMean)
  }
  results[l, 1] = signif(mean(abs(bart - df$actual) / df$actual), 3)
  results[l, 2] = signif(sd(abs(bart - df$actual) / df$actual) / sqrt(20), 3)
  results[l, 3] = signif(mean(abs(gp - df$actual) / df$actual), 3)
  results[l, 4] = signif(sd(abs(gp - df$actual) / df$actual) / sqrt(20), 3)
  results[l, 5] = signif(mean(abs(mi - df$actual) / df$actual), 3)
  results[l, 6] = signif(sd(abs(mi - df$actual) / df$actual) / sqrt(20),3)
  
  l = l+1
}
write.csv(results, "results/genz/high_dimensionality.csv")

