FuNopaCl <- function(X, classes, semimetric = "Deriv", semimetric.params) {
  z <- c()
  sm <- paste("Semimetric", semimetric, sep = "")
  z$Semimetric <- sm
  z$semimetric.params <- semimetric.params
  
  Dist <- Semimetric (X, X, z$Semimetric, z$semimetric.params)
  DistMat <- Dist$semimetric
  
  k <- .Call ("KernelClassificationkNNlCV", 
              DistMat, 
              classes, 
              20, 
              PACKAGE = "nfda")
  z$k.opt <- k$kopt
  z$mse.learn <- mean(k$classes.estimated != classes)
  z$classes.estimated <- k$classes.estimated
  z$Prob.estimated <- k$Prob.estimated
  class(z) <- "FuNopaCl"
  z$X.learn <- X
  z$classes.learn <- classes
  z
}