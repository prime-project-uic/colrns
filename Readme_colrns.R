################################
### Install package 'colrns' ###
################################
### Install packages
#install.packages("usethis")
#install.packages("devtools")
#install.packages("glmnet")
#devtools::install_github("klutometis/roxygen")
library(usethis)
library(devtools)
library(glmnet)
library(roxygen2)

### Set your own working directory
setwd("C:/Users/jiyeo/Dropbox/PhD/Research/R/Package") 

### Put 'colrns' folder under your working directory 

### Create the documents for the functions
setwd("./colrns")
document()

### Install the package
setwd("..")
install("colrns")

### Check the help page for the functions
#In case of error, try after restarting R with the commend ".rs.restartR()"
?colrns
?colrns.glm
?predict.colrns
?predict.colrns.glm


#########################
### Function Examples ###
#########################
### Function "colrns"
### NANES Persistent Organic Pollutant (POPs) Dataset (available in the github 'colrns' repository)
nhs.data <- read.csv("nhanes.pops.data.csv",sep=",", header=TRUE)
y = as.vector(nhs.data[,1]) #log-transformed longer leukocyte telomere length
x = as.matrix(nhs.data[,2:37]) #18 log-transformed POPs exposures and 18 confounders
### Fit the model (not penalizing confounders)
set.seed(1)
fit.colrns = colrns(y=y, x=x, penalty.factor= c(rep(1,18),rep(0,18)))
fit.colrns 

### Function "predict.colrns"
predict.colrns <- function(object, newx){
  predictor <- cbind(rep(1, nrow(newx)),newx) %*% object$coefficients
  return(predictor)
}
predict.colrns(object=fit.colrns, newx=x[1:5,])

### Function "colrns.glm"
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
fit.colrns.glm = colrns.glm(x, y, lambda.min=FALSE)
fit.colrns.glm

### Function "predict.colrns.glm"
predict.colrns.glm <- function(object, newx, type="response"){
  predictor <- cbind(rep(1, nrow(newx)),newx) %*% object$coefficients
  p.hat <- exp(predictor)/(1+exp(predictor))
  y.hat <- rbinom(nrow(newx), 1, p.hat)
  if(type=="response"){
    return(y.hat)
  } else if(type=="link"){
    return(predictor)
  }
}
predict.colrns.glm(object=fit.colrns.glm, newx=x[1:10,])
predict.colrns.glm(object=fit.colrns.glm, newx=x[1:10,], type="link")

