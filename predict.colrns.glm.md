# `predict.colrns.glm`

Make predictions from a "colrns.glm" object


## Description

This function makes predictions from a "colrns.glm" object


## Usage

```r
predict(object, newx, type = "response")
```


## Arguments

Argument      |Description
------------- |----------------
`object`     |     Fitted "colrns.glm" object.
`newx`     |     Matrix of new values for x at which predictions are to be made.
`type`     |     Type "link" gives the linear predictors. Type "response" gives the fitted probabilities. Default is "response".


## Value

fitted probabilities or linear predictors


## Examples

```r
### Simulated Data
library(MASS)
set.seed(1)
n=500
p=30
Cov.x <- matrix(0.2, ncol=p, nrow=p)
gp.cov <- matrix(0.6, ncol=10, nrow=10)
diag(gp.cov) <- 1
for(i in 1:3){
Cov.x[{10*(i-1)+1}:{10*i},{10*(i-1)+1}:{10*i}] <- gp.cov
}
x_mu = matrix(0,p,1)
x = mvrnorm(n, x_mu, Cov.x)
colnames(x)<-paste("x",1:p,sep="")
beta=rep(c(1,1,1,rep(0,7)),3)
prob=exp(x%*%beta)/(1+exp(x%*%beta))
y=rbinom(n, 1, prob)
### Fit the model
colrns.glm(x, y, lambda.min=FALSE)
### Predict
predict.colrns.glm(object=fit.colrns.glm, newx=x[1:10,])
predict.colrns.glm(object=fit.colrns.glm, newx=x[1:10,], type="link")
```


