library(bcf)
library(dbarts)

# # Toy example in the package manual
# p = 3 #two control variables and one moderator
# n = 250
# set.seed(1)

# x = matrix(rnorm(n*p), nrow=n)
# # create targeted selection
# q = -1*(x[,1]>(x[,2])) + 1*(x[,1]<(x[,2]))
# # generate treatment variable
# pi = pnorm(q)
# z = rbinom(n,1,pi)
# # tau is the true (homogeneous) treatment effect
# tau = (0.5*(x[,3] > -3/4) + 0.25*(x[,3] > 0) + 0.25*(x[,3]>3/4))
# # generate the response using q, tau and z
# mu = (q + tau*z)
# # set the noise level relative to the expected mean function of Y
# sigma = diff(range(q + tau*pi))/8
# # draw the response variable with additive error
# y = mu + sigma*rnorm(n)
# # If you didn't know pi, you would estimate it here
# pihat = pnorm(q)
# bcf_fit = bcf(y, z, x, x, pihat, nburn=2000, nsim=2000)
# # Get posterior of treatment effects
# tau_post = bcf_fit$tau
# tauhat = colMeans(tau_post)
# plot(tau, tauhat); abline(0,1)



# Experiment1
g <- function(x) {
  results <- length(x)
  results[x == 1] <- 2
  results[x == 2] <- -1
  results[x == 3] <- -4
  return(results)
}
set.seed(1)
p = 5
n = 250

x = matrix(rnorm(n*p), ncol = p)
x[, 4] = x[, 4] > 1
x[, 5] = sample(c(1, 2, 3), n, replace = TRUE)
x_input = makeModelMatrixFromDataFrame(data.frame(x[, 1:3], as.factor(x[, 4]), as.factor(x[, 5])))

tau = 3
mu = 1 + g(x[, 5]) + x[, 1] * x[, 3]

s = sd(mu)
pi = 0.8 * pnorm(3 * mu / s - 0.5 * x[, 1]) + 0.05 + runif(n) / 10
z = rbinom(n, rep(1, n), pi)

sigma = 1
y = mu + tau * z + sigma * rnorm(n)

pihat = pi

bcf_fit = bcf(y, z, x_input, x_input, pihat, nburn=2000, nsim=2000)
# Get posterior of treatment effects
tau_post = bcf_fit$tau
tauhat = colMeans(tau_post)
plot(tau, tauhat); abline(0,1)

