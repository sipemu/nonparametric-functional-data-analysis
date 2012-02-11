predict.FuNopaCl <- function(object, 
                             newdata, ...) {
  
  Dist <- Semimetric (object$X.learn, 
                      newdata, 
                      object$Semimetric, 
                      object$semimetric.params)
  DistMat <- Dist$semimetric  
  Y <- .Call ("KernelClassificationkNN", 
              DistMat, 
              object$classes.learn, 
              object$k.opt, 
              PACKAGE = "nfda")  
  object$classes.pred <- Y$classes.predicted
  object$Prob.pred <- Y$Prob.predicted
  object
}