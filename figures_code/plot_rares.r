# define string formatting
`%--%` <- function(x, y) 
  # from stack exchange:
  # https://stackoverflow.com/questions/46085274/is-there-a-string-formatting-operator-in-r-similar-to-pythons
{
  do.call(sprintf, c(list(x), y))
}


results <- data.frame(matrix(0, ncol = 6, nrow = 3))
colnames(results) <- c("BART-Int", "BART-Int-se", "GPBQ", "GPBQ-se", "MI", "MI-se")
dims = c(30, 40, 50)
l = 1
for (dim in dims) {
  bart = c()
  gp = c()
  mi = c()
  for (num_cv in 1:20) {
    df <- read.csv("results/rares/3/portfolio_lossDim%sExponential_%s.csv" %--% c(dim, num_cv))
    bart = c(bart, df$BARTMean)
    gp = c(gp, df$GPMean)
    mi = c(mi, df$MIMean)
  }
  if (dim == 1) {actual<-1e-05}
  if (dim == 3) {actual <- 1e-05}
  if (dim == 5) {actual <- 0.0171328}
  if (dim == 10) {actual<-0.0681249}
  if (dim == 20) {actual<-0.068942}
  if (dim == 30) {actual<-0.068855}
  results[l, 1] = signif(mean(abs(bart - actual) / actual), 3)
  results[l, 2] = signif(sd(abs(bart - actual) / actual) / sqrt(20), 3)
  results[l, 3] = signif(mean(abs(gp - actual) / actual), 3)
  results[l, 4] = signif(sd(abs(gp - actual) / actual) / sqrt(20), 3)
  results[l, 5] = signif(mean(abs(mi - actual) / actual), 3)
  results[l, 6] = signif(sd(abs(mi - actual) / actual) / sqrt(20),3)
  
  l = l+1
}
print(results)
write.csv(results, "results/rares/rares.csv")

