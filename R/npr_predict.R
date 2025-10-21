#' Predict Method for Functional Nonparametric Regression
#'
#' Predicts response values for new functional data using a trained FuNopaRe model
#'
#' @param object Object of class FuNopaRe from FuNopaRe function
#' @param newdata Matrix of new functional data to predict
#' @param bootstrap Logical indicating whether to compute bootstrap confidence intervals
#' @param bootstrap_resampling_method Character string specifying resampling method:
#'   "wild.continuous", "wild.twopoint", or "homoscedatic"
#' @param bootstrap_samples Integer number of bootstrap samples
#' @param bootstrap_neighbours Integer number of neighbours for bandwidth selection
#' @param bootstrap_alpha Significance level for confidence intervals
#' @param ... Additional arguments (not currently used)
#' @return Vector of predicted values
#' @export
#' @method predict FuNopaRe
#' @examples
#' \dontrun{
#'   # Example usage
#'   X <- matrix(rnorm(100 * 50), 100, 50)
#'   Y <- rnorm(100)
#'   semimetric <- semiemetric_deriv$new(q = 2, nknot = 10, range.grid = c(0, 1))
#'   model <- FuNopaRe(X, Y, semimetric, k_nearest_neighbors = 20L)
#'   
#'   X_new <- matrix(rnorm(10 * 50), 10, 50)
#'   predictions <- predict(model, X_new)
#' }
predict.FuNopaRe <- function(object, 
                             newdata, 
                             bootstrap = FALSE, 
                             bootstrap_resampling_method = c("wild.continuous", "wild.twopoint", "homoscedatic"),
                             bootstrap_samples = 200,
                             bootstrap_neighbours = 20,
                             bootstrap_alpha = 0.05,
                             ...) {
  # Check if object is a valid FuNopaRe object
  if (!inherits(object, "FuNopaRe")) {
    stop("object must be a valid FuNopaRe object")
  }
  # Check if newdata is a matrix
  if (!is.matrix(newdata)) {
    stop("newdata must be a matrix")
  }
  # Check if bootstrap is a logical
  if (!is.logical(bootstrap)) {
    stop("bootstrap must be a logical")
  }

  # Compute distance matrix between learning and new data
  distance_matrix <- object$semimetric$calculate(object$X, newdata)
  if (object$bandwidth_method == "CV") {
    # Use fixed bandwidth
    Y <- KernelPrediction(distance_matrix, 
                         object$Y, 
                         object$h_opt)
  } else if (object$bandwidth_method == "kNNgCV") {
    # Use k-NN with global selection
    Y <- KernelPredictionkNN(distance_matrix, 
                             object$Y, 
                             object$k_opt,
                             FALSE)
  } else if (object$bandwidth_method == "kNNlCV") {
    # Use k-NN with local bandwidth
    Y <- KernelPredictionkNN(distance_matrix, 
                             object$Y, 
                             object$k_opt,
                             TRUE)
  }
  
  if (bootstrap) {
    if (is.null(bootstrap_parameters)) {
      stop("bootstrap_parameters must be provided when bootstrap = TRUE")
    }
    
    # Generate bootstrap samples
    W <- BootstrapData(object$Y, 
                      object$Yhat, 
                      bootstrap_resampling_method, 
                      bootstrap_samples)
    
    # Compute bootstrap predictions
    R <- KernelPredictionBoot(distance_matrix,
                             object$Y,
                             Y,
                             W,
                             bootstrap_neighbours)
    
    Y <- R$pred
    z <- qnorm(bootstrap_alpha / 2)
    lower <- R$mu - z * sqrt(R$sigma / length(object$Yhat))
    upper <- R$mu + z * sqrt(R$sigma / length(object$Yhat))
  }
  
  return(list(yHat = Y, lower = lower, upper = upper))
}
