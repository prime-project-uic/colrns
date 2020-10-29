#' Make predictions from a "colrns" object
#'
#' This function makes predictions from a "colrns" object
#'
#' @param object Fitted "colrns" object.
#' @param newx Matrix of new values for x at which predictions are to be made.
#' @return fitted values
#' @keywords prediction
#' @export
#' @examples
#' ### NANES Persistent Organic Pollutant (POPs) Dataset (available in the github 'colrns' repository)
#' nhs.data <- read.csv("nhanes.pops.data.csv",sep=",", header=TRUE)
#' y = as.vector(nhs.data[,1]) #log-transformed longer leukocyte telomere length
#' x = as.matrix(nhs.data[,2:37]) #18 log-transformed POPs exposures and 18 confounders
#' ### Fit the model (not penalizing confounders)
#' fit.colrns = colrns(y=y, x=x, penalty.factor= c(rep(1,18),rep(0,18)))
#' ### Predict
#' predict.colrns(object=fit.colrns, newx=x[1:5,])


predict.colrns <- function(object, newx){
  predictor <- cbind(rep(1, nrow(newx)),newx) %*% object$coefficients
  return(predictor)
}
