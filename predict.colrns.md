# `predict.colrns`

Make predictions from a "colrns" object


## Description

This function makes predictions from a "colrns" object


## Usage

```r
list(list("predict"), list("colrns"))(object, newx)
```


## Arguments

Argument      |Description
------------- |----------------
`object`     |     Fitted "colrns" object.
`newx`     |     Matrix of new values for x at which predictions are to be made.


## Value

fitted values


## Examples

```r
### NANES Persistent Organic Pollutant (POPs) Dataset (available in the github 'colrns' repository)
nhs.data <- read.csv("nhanes.pops.data.csv",sep=",", header=TRUE)
y = as.vector(nhs.data[,1]) #log-transformed longer leukocyte telomere length
x = as.matrix(nhs.data[,2:37]) #18 log-transformed POPs exposures and 18 confounders
### Fit the model (not penalizing confounders)
fit.colrns = colrns(y=y, x=x, penalty.factor= c(rep(1,18),rep(0,18)))
### Predict
predict.colrns(object=fit.colrns, newx=x[1:5,])
```


