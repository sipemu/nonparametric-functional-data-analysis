Semimetric <- function(Data1, Data2, semimetric, semimetric.params) {
  
  if (semimetric == "SemimetricDeriv") {
    Dist <- SemimetricDeriv (Data1, 
                             Data2, 
                             q = semimetric.params$q, 
                             nknot = semimetric.params$nknot, 
                             range.grid = semimetric.params$range.grid,
                             Hhalf = semimetric.params$HhalfDeriv)
  } else if (semimetric == "SemimetricPCA") {
    Dist <- SemimetricPCA (Data1, 
                           Data2, 
                           q = semimetric.params$q,
                           EigenVec = semimetric.params$EigenVec)
  } else if (semimetric == "SemimetricPLS") {
    Dist <- SemimetricPLS (Data1, 
                           Data2, 
                           q = semimetric.params$q)
  } else {
    
  }
  return(Dist)
}

################################################################################
#
# Semimetric based on the L_2-norm of the qth derivatives
#
################################################################################
SemimetricDeriv <- function (Data1, Data2, 
                             q = 2, nknot = 20, range.grid = c(0,1), 
                             Hhalf = NULL) {
    
    library (splines)
    if (is.vector(Data1)) Data1 <- as.matrix(t(Data1))
    if (is.vector(Data2)) Data2 <- as.matrix(t(Data2))
    testfordim <- sum (dim (Data1) == dim (Data2)) == 2
    twodatasets <- 1
    if (testfordim) twodatasets <- sum (Data1 == Data2) != prod (dim (Data1))
    #####################################################################
    # B-spline approximation of the curves containing in DATASET
    #####################################################################
    p <- ncol (Data1)
    a <- range.grid[1]
    b <- range.grid[2]
    x <- seq (a, b, length = p)
    order.Bspline <- q + 3
    nknotmax <- (p - order.Bspline - 1) %/% 2
    if (nknot > nknotmax) {
      stop (paste ("give a number nknot smaller than ", 
                   nknotmax, " for avoiding ill-conditioned matrix"))
    }
    Knot <- seq (a, b, length = nknot + 2)[ - c(1, nknot + 2)]
    delta <- sort (c(rep(c(a, b), order.Bspline), Knot))
    Bspline <- splineDesign (delta, x, order.Bspline)
    #######################################################################
    # Numerical integration by the Gauss method 
    #######################################################################   
    point.gauss <- c(-0.9324695142, -0.6612093865, -0.2386191861, 
                     0.2386191861, 0.6612093865, 0.9324695142)
    weight.gauss <- c(0.1713244924, 0.360761573, 0.4679139346, 
                      0.4679139346, 0.360761573, 0.1713244924)
    x.gauss <- 0.5 * ((b + a) + (b - a) * point.gauss)
    lx.gauss <- length(x.gauss)
    
    Bspline.deriv <- splineDesign (delta, 
                                   x.gauss, 
                                   order.Bspline, 
                                   rep(q, lx.gauss))
    if (is.atomic(Hhalf)) {
    Hhalf <- .Call ("SemimetricDerivDesign", 
                    Bspline.deriv, 
                    weight.gauss, 
                    b - a, 
                    PACKAGE = "nfda")
    }
    if (twodatasets) {      
      semimetric <- .Call ("SemimetricDeriv", 
                           Hhalf, Bspline, 
                           Data1, Data2,
                           twodatasets,
                           PACKAGE = "nfda")
    } else {    
      semimetric <- .Call ("SemimetricDeriv", 
                           Hhalf, Bspline, 
                           Data1, Data1,
                           twodatasets,
                           PACKAGE = "nfda")
    }
    return(list(semimetric = semimetric, Hhalf = Hhalf))
}
################################################################################
#
# Semimetric based on Fourier coefficients
#
################################################################################



################################################################################
#
# Semimetric based on PCA
#
################################################################################
SemimetricPCA <- function (Data1, Data2, q, EigenVec = NULL) {
  
  if (is.vector(Data1)) Data1 <- as.matrix(t(Data1))
  if (is.vector(Data2)) Data2 <- as.matrix(t(Data2))
	testfordim <- sum(dim(Data1)==dim(Data2))==2
	twodatasets <- 1
	if (testfordim) twodatasets <- sum(Data1==Data2)!=prod(dim(Data1))
	qmax <- ncol(Data1)
	if (q > qmax) stop(paste("give a integer q smaller than ", qmax))
  if (is.atomic(EigenVec)) {
      EigenVec <- .Call ("SemimetricPCAEV", 
                          Data1,
                          q, 
                          PACKAGE = "nfda")
  }
  if (twodatasets) {
    semimetric <- .Call ("SemimetricPCA",  
                         EigenVec,
                         Data1, Data2,
                         twodatasets,
                         PACKAGE = "nfda")
  } else {
    semimetric <- .Call ("SemimetricPCA",  
                         EigenVec,
                         Data1, Data1,
                         twodatasets,
                         PACKAGE = "nfda")
  }
  return(list(semimetric = semimetric, EigenVec = EigenVec))
}
################################################################################
#
# Semimetric based on PLS
#
################################################################################
SemimetricPLS <- function (Response, Data1, Data2, q) {
  
  library(pls)
  if (is.vector(Data1)) Data1 <- as.matrix(t(Data1))
  if (is.vector(Data2)) Data2 <- as.matrix(t(Data2))
  testfordim <- sum (dim (Data1) == dim (Data2)) == 2
  twodatasets <- 1
  if (testfordim) twodatasets <- sum (Data1 == Data2) != prod (dim (Data1))
  qmax <- ncol(Data1)
	if(q > qmax) stop(paste("give a integer q smaller than ", qmax))
  
  n <- nrow(Data1)
  m <- ncol(Data1)
  mplsr.res <- plsr (Response ~ Data1, ncomp = q)
  Coef <- matrix (mplsr.res$coefficients, m, n)
  b0 <- mplsr.res$Ymeans - t(Coef) %*% mplsr.res$Xmeans
  
  Comp1 <- Data1 %*% Coef
  Comp1 <- outer (rep(1, n), as.vector(b0)) + Comp1
  if(twodatasets) {
    n2 <- nrow (Data2)
    Comp2 <- Data2 %*% Coef
    Comp2 <- outer (rep(1, n2), as.vector(b0)) + Comp2
  }
  else {
    Comp2 <- Comp1
  }
  Semimetric <- 0
  for(g in 1:q)
    Semimetric <- Semimetric + outer (Comp1[ , g], Comp2[ , g], "-")^2
  return(sqrt(Semimetric))
}