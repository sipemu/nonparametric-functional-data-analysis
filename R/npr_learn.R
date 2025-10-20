#' Functional Nonparametric Regression Learning
#'
#' Performs nonparametric regression for functional data using kernel methods
#'
#' @param X Matrix of functional data (each row is an observation)
#' @param Y Numeric vector of response values
#' @param semimetric Character string specifying semimetric type (default: "Deriv")
#' @param semimetric.params List of parameters for the semimetric
#' @param bandwidth Character string specifying bandwidth selection method:
#'   "CV" for cross-validation, "kNNgCV" for k-NN global CV, "kNNlCV" for k-NN local CV
#' @return Object of class FuNopaRe containing:
#'   \itemize{
#'     \item Semimetric: Name of the semimetric used
#'     \item Method: Method used for prediction
#'     \item semimetric.params: Parameters used for the semimetric
#'     \item h.opt or k.opt: Optimal bandwidth or k value
#'     \item mse.learn: Learning mean squared error
#'     \item Y.hat: Fitted values
#'     \item X.learn: Learning data
#'     \item Y.learn: Learning response values
#'   }
#' @export
#' @examples
#' \dontrun{
#'   # Example usage
#'   X <- matrix(rnorm(100 * 50), 100, 50)
#'   Y <- rnorm(100)
#'   params <- list(q = 2, nknot = 10, range.grid = c(0, 1))
#'   model <- FuNopaRe(X, Y, "Deriv", params, "kNNgCV")
#' }
FuNopaRe <- function(X, 
                     Y, 
                     semimetric = "Deriv", 
                     semimetric.params, 
                     bandwidth = "CV") {
  
  # Initialize result object as a list
  z <- list()
  sm <- paste("Semimetric", semimetric, sep = "")
  method <- paste("KernelPrediction", bandwidth, sep = "")
  z$Semimetric <- sm
  z$Method <- method
  z$semimetric.params <- semimetric.params
  
  # Compute distance matrix
  Dist <- Semimetric(X, X, z$Semimetric, z$semimetric.params)
  DistMat <- Dist$semimetric
  
  if (method == "KernelPredictionCV") {
    # Cross-validation bandwidth selection
    band <- quantile(DistMat[row(DistMat) > col(DistMat)], 0.05)
    h <- KernelPredictionCV(DistMat, Y, band)
    
    z$h.opt <- h$hopt
    z$mse.learn <- h$mse
    z$hseq <- h$hseq
    z$Y.hat <- h$yhat
    
  } else if (method == "KernelPredictionkNNgCV") {
    # k-NN global cross-validation
    k <- KernelPredictionkNNgCV(DistMat, Y, 20L)
    
    z$k.opt <- k$kopt
    z$mse.learn <- k$mse
    z$Y.hat <- k$yhat
    
  } else if (method == "KernelPredictionkNNlCV") {
    # k-NN local cross-validation
    k <- KernelPredictionkNNlCV(DistMat, Y, 20L)
    
    z$k.opt <- as.vector(k$kopt)
    z$mse.learn <- k$mse
    z$Y.hat <- k$yhat
    
  } else {
    stop(paste("Unknown bandwidth selection method:", bandwidth))
  }
  
  class(z) <- "FuNopaRe"
  z$X.learn <- X
  z$Y.learn <- Y
  
  return(z)
}
