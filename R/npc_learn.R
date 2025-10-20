#' Functional Nonparametric Classification Learning
#'
#' Performs nonparametric classification for functional data using k-NN
#' with local cross-validation
#'
#' @param X Matrix of functional data (each row is an observation)
#' @param classes Integer vector of class labels
#' @param semimetric Character string specifying semimetric type (default: "Deriv")
#' @param semimetric.params List of parameters for the semimetric
#' @return Object of class FuNopaCl containing:
#'   \itemize{
#'     \item Semimetric: Name of the semimetric used
#'     \item semimetric.params: Parameters used for the semimetric
#'     \item k.opt: Optimal k values from cross-validation
#'     \item mse.learn: Learning error rate
#'     \item classes.estimated: Estimated classes for learning set
#'     \item Prob.estimated: Probability matrix for learning set
#'     \item X.learn: Learning data
#'     \item classes.learn: Learning class labels
#'   }
#' @export
#' @examples
#' \dontrun{
#'   # Example usage
#'   X <- matrix(rnorm(100 * 50), 100, 50)
#'   classes <- rep(1:2, each = 50)
#'   params <- list(q = 2, nknot = 10, range.grid = c(0, 1))
#'   model <- FuNopaCl(X, classes, "Deriv", params)
#' }
FuNopaCl <- function(X, classes, semimetric = "Deriv", semimetric.params) {
  
  # Initialize result object as a list
  z <- list()
  sm <- paste("Semimetric", semimetric, sep = "")
  z$Semimetric <- sm
  z$semimetric.params <- semimetric.params
  
  # Compute distance matrix
  Dist <- Semimetric(X, X, z$Semimetric, z$semimetric.params)
  DistMat <- Dist$semimetric
  
  # Perform k-NN classification with local cross-validation
  k <- KernelClassificationkNNlCV(DistMat, classes, 20L)
  
  z$k.opt <- k$kopt
  z$mse.learn <- mean(k$classes.estimated != classes)
  z$classes.estimated <- k$classes.estimated
  z$Prob.estimated <- k$Prob.estimated
  z$X.learn <- X
  z$classes.learn <- classes
  
  class(z) <- "FuNopaCl"
  
  return(z)
}
