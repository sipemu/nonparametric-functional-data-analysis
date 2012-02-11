FuNopaRe <- function(X, 
                     Y, 
                     semimetric = "Deriv", 
                     semimetric.params, 
                     bandwidth = "CV") {
  
  z <- c()
  sm <- paste("Semimetric", semimetric, sep = "")
  method <- paste("KernelPrediction", bandwidth, sep = "")
  z$Semimetric <- sm
  z$Method <- method
  z$semimetric.params <- semimetric.params
  
  Dist <- Semimetric (X, X, z$Semimetric, z$semimetric.params)
  DistMat <- Dist$semimetric
  
  if (method == "KernelPredictionCV") {
    band <- quantile(DistMat[row(DistMat) > col(DistMat)], 0.05)
    h <- .Call ("KernelPredictionCV", 
                DistMat, 
                Y, 
                band, 
                PACKAGE = "nfda")  
  
    z$h.opt <- h$hopt
    z$mse.learn <- h$mse
    z$hseq <- h$hseq
    z$Y.hat <- h$yhat
  } else if (method == "KernelPredictionkNNgCV") {
    k <- .Call ("KernelPredictionkNNgCV", 
                DistMat, 
                Y, 
                20, 
                PACKAGE = "nfda")
    z$k.opt <- k$kopt
    z$mse.learn <- k$mse
    z$Y.hat <- k$yhat
  } else if (method == "KernelPredictionkNNlCV") {
    k <- .Call ("KernelPredictionkNNlCV", 
                DistMat, 
                Y, 
                20, 
                PACKAGE = "nfda")
    z$k.opt <- as.vector(k$kopt)
    z$mse.learn <- k$mse
    z$Y.hat <- k$yhat
  }
  class(z) <- "FuNopaRe"
  z$X.learn <- X
  z$Y.learn <- Y
  z
}