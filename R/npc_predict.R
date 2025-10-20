#' Predict Method for Functional Nonparametric Classification
#'
#' Predicts class labels for new functional data using a trained FuNopaCl model
#'
#' @param object Object of class FuNopaCl from FuNopaCl function
#' @param newdata Matrix of new functional data to predict
#' @param ... Additional arguments (not currently used)
#' @return Updated FuNopaCl object with additional fields:
#'   \itemize{
#'     \item classes.pred: Predicted class labels
#'     \item Prob.pred: Probability matrix for predictions
#'   }
#' @export
#' @method predict FuNopaCl
#' @examples
#' \dontrun{
#'   # Example usage
#'   X <- matrix(rnorm(100 * 50), 100, 50)
#'   classes <- rep(1:2, each = 50)
#'   params <- list(q = 2, nknot = 10, range.grid = c(0, 1))
#'   model <- FuNopaCl(X, classes, "Deriv", params)
#'   
#'   X_new <- matrix(rnorm(10 * 50), 10, 50)
#'   predictions <- predict(model, X_new)
#' }
predict.FuNopaCl <- function(object, newdata, ...) {
  
  # Compute distance matrix between learning and new data
  Dist <- Semimetric(object$X.learn, 
                     newdata, 
                     object$Semimetric, 
                     object$semimetric.params)
  DistMat <- Dist$semimetric
  
  # Perform k-NN classification
  Y <- KernelClassificationkNN(DistMat, 
                               object$classes.learn, 
                               object$k.opt)
  
  object$classes.pred <- Y$classes.predicted
  object$Prob.pred <- Y$Prob.predicted
  
  return(object)
}
