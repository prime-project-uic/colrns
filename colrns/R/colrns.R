#' Correlation Learning for Variable Selection (COLRNS)
#'
#' This function allows you to implement COLRNS for high dimensional data analysis with continuous responses.
#' This package is based on the 'glmnet' package (Friedman et al., 2019)
#'
#' @param y Response vector (continuous)
#' @param x Design matrix including predictors and confounders
#' @param penalty.factor Separate penalty factor can be applied to each coefficient. It can be 0 for some variables. If 0, corresponding variables are always included in the model. Default is 1 for all variables.
#' @param nfolds Number of folds that is used in cross-validation. Default is 10.
#' @param type.measure Measurement to be optimized in cross-validation. Default is "deviance", which uses squared-error for gaussian models (a.k.a type.measure="mse"). Mean absolute error "mae" is another option.
#' @param alpha Grid for the parameter alpha, the proportion between the lasso penalty (alpha=1) and the ridge penalty (alpha=0). Default is from 0.1 to 0.9 with the interval of 0.1.
#' @param lambda.min Whether to estimate the coefficients at lambda.min or lambda.1se. Using lambda.1se tends to have a more sparse solution. Default is lambda.min.
#' @param cutoff Cutoff to construct the lead cluster. Default is 0.5.
#' @param r.group Minimum ratio of predictors to be included in the lead cluster. The minimum number of predictors recruited in the lead cluster d can be n-1 or n/log(n) (Fan and Lv, 2008), or 2n/log(n) (Zhong and Zhu, 2015). The minimum ratio can be d/number of total predictors*100. When sample size n is high, we recommend to use between 0.2 and 0.4. Default is 0.3.
#' @param niter Number of iterative procedures. Default is 150.
#' @return coefficients: A named vector of coefficients.
#' @return select.vars: A list of selected variables, that is, the variables with non-zero coefficients.
#' @keywords Variable selection; Correlation learning; Collinearity; High dimensional data
#' @export
#' @examples
#' ### NANES Persistent Organic Pollutant (POPs) Dataset (available in the github 'colrns' repository)
#' nhs.data <- read.csv("nhanes.pops.data.csv",sep=",", header=TRUE)
#' y = as.vector(nhs.data[,1]) #log-transformed longer leukocyte telomere length
#' x = as.matrix(nhs.data[,2:37]) #18 log-transformed POPs exposures and 18 confounders
#' ### Fit the model (not penalizing confounders)
#' colrns(y=y, x=x, penalty.factor= c(rep(1,18),rep(0,18)))
#' @references Hastie, T., Tibshirani, R., and Friedman, J. (2009). The Elements of Statistical Learning: Data Mining, Inference, and Prediction, 2nd Edition. Springer, New York
#' @references Hastie, T., Tibshirani, R., and Wainwright, M. (2015). Statistical Learning with Sparsityy: the LASSO and Generalizations. Boca Raton, FL, USA: CRC Press.
#' @references Friedman, J., Hastie, T., Tibshirani, R., Narasimhan, B., Simon, N., and Qian, J. (2019). Lasso and elastic-net regularized generalized linear models.
#' @references Jang, J. (2020). Variable Selection in Presence of Strong Collinearity with Application to Environmental Mixtures. PhD Dissertation, University of Illinois at Chicago


colrns <- function(x, y, penalty.factor=rep(1, nvars), nfolds=10, type.measure="deviance", alpha=seq(0.1,0.9,by=0.1), lambda.min=TRUE, cutoff=0.5, r.group=0.3, niter=150){

  if(is.null(x)==TRUE){
    print("Missing x matrix")
  }
  if(is.null(y)==TRUE){
    print("Missing y vector")
  }

  stopifnot(length(y) > 0, is.numeric(y), anyNA(y) == FALSE)
  stopifnot(is.numeric(x), nrow(x) == length(y), anyNA(x) == FALSE)

  if(is.matrix(x)==FALSE){
    print("x should be a matrix")
  }
  if(is.vector(y)==FALSE){
    print("y should  be a vector")
  }

  nvars <- ncol(x)
  is_penalized = penalty.factor
  iter <- 0
  list.beta.colrns<-list()
  list.x.select<-list()
  list.fit.object <- list()
  vec.alpha <- rep(NA, length=niter)
  error <- rep(NA, length=niter)
  int.mt <- matrix(1,nrow=nrow(x),ncol=1)

  mydata.train <- cbind(y,x[,which(is_penalized!=0)])
  y.train <- y
  x.train <- x[,which(is_penalized!=0)]
  x.fixed <- x[,which(is_penalized==0)]
  y.cor0 <- abs(cor(mydata.train)[1,-1])
  x.cor0 <- abs(cor(mydata.train)[-1,-1])
  set.red <- colnames(x.train)
  lead.var.mat = matrix(0, nrow=1, ncol=niter)
  p <- ncol(x.train)

  lead.var <- names(which.max(y.cor0[set.red]))
  group.var <- group.var0 <- group.var0.def <- colnames(x.train)[(x.cor0[lead.var, ]>cutoff)]

  ## Check if the group size of xj's > r.group percent of total number of variables
  iter1=0
  if(length(group.var)/p<r.group){
    while(length(group.var)/p<r.group){
      iter1 = iter1+1
      group.var <- group.var0 <- group.var0.def <-colnames(x.cor0)[x.cor0[lead.var, ]>(cutoff-0.1*iter1)]
    }
  }

  if(ncol(x.fixed)==0 & length(group.var)==1){
    print("Increase the 'r.group' so that the lead cluster can include more than one varaible") #09122020 should p*r.group > 1
  }

  while (iter<niter) {

    lead.var.mat[,iter+1] <- lead.var
    iter = iter + 1

    is_cluster = c(rep(1,length(group.var)),rep(0,ncol(x.fixed)))
    cv.error <- NULL; cv.list <- NULL;  cv.elst1 <- NULL
    for (i in 1:length(alpha)){
      cv.list[[i]] <- cv.glmnet(y=y.train, x=cbind(x.train[,group.var],x.fixed), type.measure=type.measure, penalty.factor = is_cluster, nfolds=nfolds, family="gaussian",alpha=alpha[i])
      cv.error <- c(cv.error, min(cv.list[[i]]$cvm))
    }

    list.fit.object[[iter]] <- cv.elst1 <- cv.glmnet(y=y.train, x=cbind(x.train[,group.var],x.fixed), type.measure=type.measure, penalty.factor = is_cluster, nfolds=nfolds, family="gaussian",alpha=alpha[which.min(cv.error)])
    vec.alpha[iter] <- alpha[which.min(cv.error)]

    beta.group <- matrix(0,nrow=length(group.var)+ncol(x.fixed)+1,ncol=1)
    if(lambda.min==FALSE){
      beta.group[coef(cv.elst1)@i+1] <- coef(cv.elst1)@x
    } else if(lambda.min==TRUE){
      beta.group[coef(cv.elst1,s="lambda.min")@i+1] <- coef(cv.elst1,s="lambda.min")@x
    }
    rownames(beta.group) <- coef(cv.elst1)@Dimnames[[1]]

    ## Residual
    xbeta <- cbind(int.mt, cbind(x.train[,group.var],x.fixed)) %*% beta.group
    resd <- data.frame(r = y.train - xbeta)

    ## Set of selected x's
    x.select <- group.var[which(beta.group[-c(1, (length(group.var)+2):length(beta.group))] != 0)]

    ## Switching the pool for new candidate
    if(length(group.var)==p){
      if(length(x.select)==p){
        set.new1 <- x.select  #when all vars are selected in previous iteration, pick new lead var among all p vars
      } else {
        set.new1 <- set.red[set.red %in% x.select == "FALSE"] #when not all vars are selected in previous iteration, pick new lead var among (all p vars\x.select)
      }
    } else {
      set.new1 <- set.red[set.red %in% c(group.var0,x.select) == "FALSE"] #New pool to select var's: set.red \ (selected x's or candidate group just now)
    }

    ## Flow within the jth simulation
    error[iter] <- min(cv.elst1$cvm)
    list.x.select[[iter]] <- c(colnames(x.fixed), x.select)
    list.beta.colrns[[iter]] <- beta.group

    ### Update 'mydata.train', 'y.cor', 'lead.var', 'x.cor', 'group' var with new residuals
    mydata.train <- data.frame(resd, x.train[,set.new1])
    y.cor <- abs(cor(mydata.train)[1,-1])
    lead.var <- names(which.max(y.cor))
    if(length(set.new1)==1){
      lead.var <- set.new1
    }
    group.var <- group.var0 <- colnames(x.cor0)[x.cor0[lead.var, ]>cutoff]

    ## Check if the group size of xj's > 5 percent of total number of variables
    iter2=0
    if(length(group.var)/p<r.group){
      while(length(group.var)/p<r.group){
        iter2 = iter2+1
        group.var <- group.var0 <- colnames(x.cor0)[x.cor0[lead.var, ]>(cutoff-0.1*iter2)]
      }
    }

    ### Combining x_selected and group of xj's (candidates)
    if(length(x.select)==0){
      group.var <- group.var0 <- group.var0.def <- unique(c(group.var, group.var0.def)) # current candidate group + previous group with no variable selected
    } else {
      # Combine a selected predictors and a new set of x's correlated with the current lead var.
      group.var <- group.var0.def <- unique(c(x.select,group.var))
    }
  }
  beta.colrns <- matrix(0,nrow=p+ncol(x.fixed)+1,ncol=1)
  rownames(beta.colrns) <- c("(Intercept)", colnames(x))
  beta.colrns[dimnames(list.beta.colrns[[which.min(error)]])[[1]],] <- list.beta.colrns[[which.min(error)]]
  select.vars <- list.x.select[[which.min(error)]]
  select.vars <- select.vars[order(match(select.vars,set.red))]
  results <- list(coefficients=beta.colrns, select.vars=select.vars)
  return(results)
}


