#' Predict Method for Functional Nonparametric Regression
#'
#' Predicts response values for new functional data using a trained FuNopaRe model
#'
#' @param object Object of class FuNopaRe from FuNopaRe function
#' @param newdata Matrix of new functional data to predict
#' @param method.params List of parameters for bootstrap (if Bootstrapping = TRUE):
#'   \itemize{
#'     \item Resampling.Method: "wild.continuous", "wild.twopoint", or "homoscedatic"
#'     \item NB: Number of bootstrap samples
#'     \item neighbours: Number of neighbours for bandwidth selection
#'     \item alpha: Significance level for confidence intervals
#'   }
#' @param Bootstrapping Logical indicating whether to compute bootstrap confidence intervals
#' @param ... Additional arguments (not currently used)
#' @return Updated FuNopaRe object with additional field:
#'   \itemize{
#'     \item Prediction: Vector of predicted values
#'     \item loConfInt: Lower confidence interval (if Bootstrapping = TRUE)
#'     \item upConfInt: Upper confidence interval (if Bootstrapping = TRUE)
#'     \item method.params.bootstrap: Bootstrap parameters (if Bootstrapping = TRUE)
#'   }
#' @export
#' @method predict FuNopaRe
#' @examples
#' \dontrun{
#'   # Example usage
#'   X <- matrix(rnorm(100 * 50), 100, 50)
#'   Y <- rnorm(100)
#'   params <- list(q = 2, nknot = 10, range.grid = c(0, 1))
#'   model <- FuNopaRe(X, Y, "Deriv", params, "kNNgCV")
#'   
#'   X_new <- matrix(rnorm(10 * 50), 10, 50)
#'   predictions <- predict(model, X_new)
#' }
predict.FuNopaRe <- function(object, 
                             newdata, 
                             method.params = NULL, 
                             Bootstrapping = FALSE, 
                             ...) {
  
  # Compute distance matrix between learning and new data
  Dist <- Semimetric(object$X.learn, 
                     newdata, 
                     object$Semimetric, 
                     object$semimetric.params)
  DistMat <- Dist$semimetric
  
  if (object$Method == "KernelPredictionCV") {
    # Use fixed bandwidth
    Y <- KernelPrediction(DistMat, 
                         object$Y.learn, 
                         object$h.opt)
    object$Prediction <- Y
    
  } else if (object$Method == "KernelPredictionkNNgCV") {
    # Use k-NN with global bandwidth
    Y <- KernelPredictionkNN(DistMat, 
                             object$Y.learn, 
                             object$k.opt,
                             FALSE)
    object$Prediction <- Y
    
  } else if (object$Method == "KernelPredictionkNNlCV") {
    # Use k-NN with local bandwidth
    Y <- KernelPredictionkNN(DistMat, 
                             object$Y.learn, 
                             object$k.opt,
                             TRUE)
    object$Prediction <- Y
  }
  
  if (Bootstrapping) {
    if (is.null(method.params)) {
      stop("method.params must be provided when Bootstrapping = TRUE")
    }
    
    # Generate bootstrap samples
    W <- BootstrapData(object$Y.learn, 
                      object$Y.hat, 
                      method.params$Resampling.Method, 
                      method.params$NB)
    
    # Compute bootstrap predictions
    R <- KernelPredictionBoot(DistMat,
                             object$Y.learn,
                             Y,
                             W,
                             method.params$neighbours)
    
    object$Prediction <- R$pred
    z <- qnorm(method.params$alpha / 2)
    object$loConfInt <- R$mu - z * sqrt(R$sigma / length(object$Prediction))
    object$upConfInt <- R$mu + z * sqrt(R$sigma / length(object$Prediction))
    object$method.params.bootstrap <- method.params
  }
  
  return(object)
}
