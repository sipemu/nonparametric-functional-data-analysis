#' Functional Nonparametric Regression Learning
#'
#' Performs nonparametric regression for functional data using kernel methods
#'
#' @param X Matrix of functional data (each row is an observation)
#' @param Y Numeric vector of response values
#' @param semimetric Semimetric object
#' @param bandwidth_method Character string specifying bandwidth selection method:
#'   "CV" for cross-validation, "kNNgCV" for k-NN global CV, "kNNlCV" for k-NN local CV
#' @param grid_of_bandwidths Numeric vector of bandwidth values for cross-validation
#' @param k_nearest_neighbors Integer specifying the number of nearest neighbors for k-NN cross-validation
#' @return Object of class FuNopaRe containing:
#'   \itemize{
#'     \item semimetric: Semimetric object
#'     \item bandwidth_method: Bandwidth selection method used
#'     \item grid_of_bandwidths: Grid of bandwidth values used for cross-validation
#'     \item k_opt: Optimal k values from cross-validation
#'     \item mse: Learning error rate
#'     \item Yhat: Estimated response values for learning set
#'     \item X: Learning data
#'     \item Y: Learning response values
#'   }
#' @export
#' @examples
#' \dontrun{
#'   # Example usage
#'   X <- matrix(rnorm(100 * 50), 100, 50)
#'   Y <- rnorm(100)
#'   semimetric <- semiemetric_deriv$new(q = 2, nknot = 10, range.grid = c(0, 1))
#'   k_nearest_neighbors <- 20L
#'   model <- FuNopaRe(X, Y, semimetric, k_nearest_neighbors)
#' }
FuNopaRe <- function(X, 
                     Y, 
                     semimetric, 
                     bandwidth_method = c("CV", "kNNgCV", "kNNlCV"),
                     grid_of_bandwidths = NULL,
                     k_nearest_neighbors = 20L) {
  bandwidth_method <- match.arg(bandwidth_method, c("CV", "kNNgCV", "kNNlCV"))
  
  # Check if semimetric is a valid semimetric object
  if (!inherits(semimetric, "semimetric")) {
    stop("semimetric must be a valid semimetric object")
  }
  
  # Initialize result object as a list
  z <- list()
  z$bandwidth_method <- bandwidth_method
  
  # Compute distance matrix
  distance_matrix <- semimetric$calculate(X, X)
  
  if (bandwidth_method == "CV") {
    # Cross-validation bandwidth selection
    if (is.null(grid_of_bandwidths)) {
      grid_of_bandwidths <- quantile(distance_matrix[row(distance_matrix) > col(distance_matrix)], 0.05)
    }
    h <- KernelPredictionCV(distance_matrix, Y, grid_of_bandwidths)
    
    z$h_opt <- h$hopt
    z$mse <- h$mse
    z$hseq <- h$hseq
    z$Yhat <- h$yhat
    
  } else if (bandwidth_method == "kNNgCV") {
    # k-NN global cross-validation
    k <- KernelPredictionkNNgCV(distance_matrix, Y, k_nearest_neighbors)
    
    z$k_opt <- k$k.opt + 1
    z$mse <- k$mse
    z$Yhat <- k$yhat
    
  } else if (bandwidth_method == "kNNlCV") {
    # k-NN local cross-validation
    k <- KernelPredictionkNNlCV(distance_matrix, Y, k_nearest_neighbors)
    
    z$k_opt <- as.vector(k$kopt) + 1
    z$mse <- k$mse
    z$Yhat <- k$yhat
    
  } else {
    stop(paste("Unknown bandwidth selection method:", bandwidth_method))
  }
  
  class(z) <- "FuNopaRe"
  z$semimetric <- semimetric
  z$X <- X
  z$Y <- Y
  
  return(z)
}
