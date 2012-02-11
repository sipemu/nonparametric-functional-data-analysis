predict.FuNopaRe <- function(object, 
                             newdata, 
                             method.params, 
                             Bootstrapping = FALSE, ...) {
  
  Dist <- Semimetric (object$X.learn, 
                      newdata, 
                      object$Semimetric, 
                      object$semimetric.params)
  DistMat <- Dist$semimetric
  if (object$Method == "KernelPredictionCV") {   
    Y <- .Call ("KernelPrediction", 
                DistMat, 
                object$Y.learn, 
                object$h.opt, 
                PACKAGE = "nfda")  
    object$Prediction <- Y
  } else if (object$Method == "KernelPredictionkNNgCV") {   
      Y <- .Call ("KernelPredictionkNN", 
                    DistMat, 
                    object$Y.learn, 
                    object$k.opt,
                    FALSE,
                    PACKAGE = "nfda")
      object$Prediction <- Y
  } else if (object$Method == "KernelPredictionkNNlCV") {   
      Y <- .Call ("KernelPredictionkNN", 
                    DistMat, 
                    object$Y.learn, 
                    object$k.opt,
                    TRUE,
                    PACKAGE = "nfda")
      object$Prediction <- Y
  }
  if (Bootstrapping == TRUE) {
    W <- BootstrapData (object$Y.learn, 
                        object$Y.hat, 
                        method.params$Resampling.Method, 
                        method.params$NB)
    R <- .Call ("KernelPredictionBoot", 
                DistMat,
                object$Y.learn,
                Y,
                W,
                method.params$neighbours,
                PACKAGE = "nfda")
    object$Prediction <- R$pred
    z <- qnorm(method.params$alpha / 2)
    object$loConfInt <- R$mu - z * sqrt(R$sigma / length(object$Prediction))
    object$upConfInt <- R$mu + z * sqrt(R$sigma / length(object$Prediction))
    object$method.params.bootstrap <- method.params
  }
  object
}