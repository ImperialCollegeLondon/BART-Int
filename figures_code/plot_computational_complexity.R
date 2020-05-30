# define string formatting
`%--%` <- function(x, y) 
  # from stack exchange:
  # https://stackoverflow.com/questions/46085274/is-there-a-string-formatting-operator-in-r-similar-to-pythons
{
  do.call(sprintf, c(list(x), y))
}

nums <- c(100,500,1000,2000,3000,7000,10000)
bart_times <- c()
mi_times <- c()
gp_times <- c()
bart_times_10 <- c()
bart_times_20 <- c()

i=1
for (num in nums) {
  filename <- "results/genz/1/computational_complexity_stepDim1NoSequentialUniform_num%s.csv" %--% num
  df <- read.csv(filename)
  bart_times[i] <- df$runtimeBART
  mi_times[i] <- df$runtimeMI
  gp_times[i] <- df$runtimeGP
  i <- i+1
}

#d=10
i=1
for (num in nums) {
  filename <- "results/genz/1/computational_complexity_stepDim10NoSequentialUniform_num%s.csv" %--% num
  df <- read.csv(filename)
  bart_times_10[i] <- df$runtimeBART
  i <- i+1
}

par(mfrow=c(1,1), pty="s")
plot(nums, gp_times, xlab="n", ylab="Time (s)", 
     main="", cex.main=1.9, cex.lab=2,lwd=6.8, cex.axis=2,
     ylim = c(0, 750), col="dodgerblue")
df <- data.frame(cbind(gp_times, nums^3))
colnames(df) <- c("gp_times", "num_cubed")
m <- lm(gp_times ~ num_cubed, data=df)
pred_df <- data.frame(seq(0, 10000, 1)^3)
colnames(pred_df) <- "num_cubed"
points(seq(0, 10000, 1), predict(m, pred_df), ty="l", col="dodgerblue")


points(nums, bart_times,lwd=6.8, col="orangered")
m <- lm(bart_times~ nums)
coeffs <-as.numeric(m$coefficients)
abline(a=coeffs[1], b=coeffs[2], col="orangered")

points(nums, bart_times_10,lwd=6.8, col="green")
m <- lm(bart_times_10~ nums)
coeffs <-as.numeric(m$coefficients)
abline(a=coeffs[1], b=coeffs[2], col="green")

legend("topleft", legend=c("GP-BQ d=1", "BART-Int d=1", "BART-Int d=10"),
       col=c("dodgerblue", "orangered", "green"), lty = c(1,1,1), cex=1.6, bty="n")
