#' Functional Nonparametric Classification Learning
#'
#' Performs nonparametric classification for functional data using k-NN
#' with local cross-validation
#'
#' @param X Matrix of functional data (each row is an observation)
#' @param Y Integer vector of class labels
#' @param semimetric Semimetric object
#' @param k_nearest_neighbors Integer specifying the number of nearest neighbors
#' @return Object of class FuNopaCl containing:
#'   \itemize{
#'     \item semimetric: Semimetric object
#'     \item k_opt: Optimal k values from cross-validation
#'     \item mse: Learning error rate
#'     \item Yhat: Estimated response values for learning set
#'     \item prob_estimated: Probability matrix for learning set
#'     \item X: Learning data
#'     \item Y: Learning response values
#'   }
#' @export
#' @examples
#' \dontrun{
#'   # Example usage
#'   X <- matrix(rnorm(100 * 50), 100, 50)
#'   Y <- rep(1:2, each = 50)
#'   semimetric <- semiemetric_deriv$new(q = 2, nknot = 10, range.grid = c(0, 1))
#'   model <- FuNopaCl(X, Y, semimetric, k_nearest_neighbors)
#' }
FuNopaCl <- function(X, Y, semimetric, k_nearest_neighbors = 20L) {
  
  # Initialize result object as a list
  z <- list()
  z$semimetric <- semimetric
  
  # Compute distance matrix
  distance_matrix <- semimetric$calculate(X, X)
  
  # Perform k-NN classification with local cross-validation
  k <- KernelClassificationkNNlCV(distance_matrix, Y, k_nearest_neighbors)
  
  z$k_opt <- k$kopt + 1
  z$mse <- mean(k$classes.estimated != Y)
  z$Yhat <- k$classes.estimated
  z$prob_estimated <- k$Prob.estimated
  z$X <- X
  z$Y <- Y
  
  class(z) <- "FuNopaCl"
  
  return(z)
}
