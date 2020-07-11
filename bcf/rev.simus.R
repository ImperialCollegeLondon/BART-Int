library(dbarts)
library(bcf)
library(grf)
library(bayeslm)
library(nnet)   # I added this


#set.seed(12193)
n = 1000
x1 = sample(c(-2,2),n,replace=TRUE) + 0.25*rnorm(n)
x2 = rbinom(n,1,0.7)
x3 = sample(1:3,n,replace=TRUE,prob = c(0.1,0.3,0.6))
x3 = as.factor(x3)

x4 = rnorm(n)
x5 = rbinom(n,10,0.5)
x = data.frame(x1,x2,x3,x4,x5)




## RIC

mu = function(x){
  lev = c(-1,2,0)
  
  # linear
 #  result = 1 + x[,1]*(2*x[,2] - 2*(1-x[,2])) + lev[x3]

  # nonlinear
  result = 10*pnorm(x[,1]/4+x[,5]/10) + lev[x[,3]]

  return(result)
}

#alpha = 0.5
alpha = 1 + x[,4]^2*(2*x[,2]-1)
alpha = alpha - min(alpha)
alpha = 1*sd(mu(x))*alpha/sd(alpha)


pi = 0.6*pnorm(mu(x)-mean(mu(x)),0,0.75*sd(mu(x))) + 0.2
#pi = 0.75
hist(pi,30)
z = rbinom(n,1,pi)

Ey = mu(x) + alpha*z

sig = .1*sd(Ey)

y = Ey + sig*rnorm(n)

x.mod = makeModelMatrixFromDataFrame(data.frame(x,z))
x.test = makeModelMatrixFromDataFrame(data.frame(rbind(x,x),z = c(rep(1,n),rep(0,n))))
fit.bart = bart(x.mod,y,x.test,verbose = FALSE)

alpha.hat.bart = fit.bart$yhat.test.mean[1:n] -  fit.bart$yhat.test.mean[1:n+n]

#fit.glm = glm(z~.,data = data.frame(x,z))
# without interactions it work bad

#fit.glm = glm(z~.^3,data = data.frame(x,z),family = binomial)
#pihat = fit.glm$fitted.values

x.mod = makeModelMatrixFromDataFrame(data.frame(x))

fitz = nnet(z~.,data = x.mod,size = 4,rang = 0.1, maxit = 1000,abstol = 1.0e-8, decay = 5e-2,trace=FALSE)
pihat = fitz$fitted.values

fit.bcf = bcf(y, z, x.mod, x.mod, pihat, 10000, 3000)
#fit.bcf = bcf(y, z, x.mod, x.mod, pihat, 3000, 3000, 1, include_pi="both")
#fit.bcf = bcf(y, z, x.mod, x.mod, pihat, 1000, 3000, 1, include_pi="both", use_mscale=FALSE, use_tauscale=FALSE)

alphas = colMeans(fit.bcf$tau)


tau.forest = causal_forest(x.mod, y, z,num.trees = 4000)
tau.hat = predict(tau.forest, x.mod,estimate.variance = TRUE)

# OLS
 fit.lm = lm(y~.^3,data = data.frame(y,x.mod,z))
 yhat = predict(fit.lm,newdata = data.frame(x.test))

# horseshoe
#fit.lm = bayeslm(y~.^3,data = data.frame(y,x.mod,z,pihat))
#X = model.matrix(y~.^3,data.frame(y=c(y,y),x.test,pihat = c(pihat,pihat)))
#yhat = X%*%colMeans(fit.lm$beta)
tau.lm = yhat[1:n]-yhat[1:n+n]

if (length(alpha) == 1){
hist(alpha.hat.bart,30)
hist(tau.hat$predictions,30, add=TRUE,col='purple') # grf
hist(alphas,30,add=TRUE,col='orange') # bcf
hist(tau.lm,30,add=TRUE,col='pink')
abline(v=alpha,col='red',lwd=3)

rmse = NULL
rmse$bcf = sqrt(mean((mean(alpha) - mean(alphas))^2))
rmse$bart = sqrt(mean((mean(alpha) - mean(alpha.hat.bart))^2))
rmse$grf = sqrt(mean((mean(alpha) - mean(tau.hat$prediction))^2))
rmse$lm = sqrt(mean((mean(alpha) - mean(tau.lm))^2))


print(rmse)

}else{
plot(alpha, alpha.hat.bart,pch=20)
points(alpha, tau.hat$predictions,pch=20,col='purple')
points(alpha, alphas,pch=20,col='orange') # bcf
points(alpha,tau.lm,pch=20,col='pink') # lm
rmse = NULL
rmse$bcf = sqrt(mean((alpha - alphas)^2))
rmse$bart = sqrt(mean((alpha - alpha.hat.bart)^2))
rmse$grf = sqrt(mean((alpha - tau.hat$prediction)^2))
rmse$lm = sqrt(mean((alpha - tau.lm)^2))

print(rmse)


}


