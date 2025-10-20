#' Semimetric Wrapper Function
#'
#' Main wrapper function for computing semimetric distances between functional data
#'
#' @param Data1 Matrix of functional data (first dataset)
#' @param Data2 Matrix of functional data (second dataset)
#' @param semimetric Character string specifying the semimetric type
#' @param semimetric.params List of parameters for the semimetric
#' @return List containing semimetric distance matrix and additional outputs
#' @export
Semimetric <- function(Data1, Data2, semimetric, semimetric.params) {
  
  if (semimetric == "SemimetricDeriv") {
    Dist <- semimetric_deriv_wrapper(Data1, 
                                     Data2, 
                                     q = semimetric.params$q, 
                                     nknot = semimetric.params$nknot, 
                                     range.grid = semimetric.params$range.grid,
                                     Hhalf = semimetric.params$HhalfDeriv)
  } else if (semimetric == "SemimetricPCA") {
    Dist <- semimetric_pca_wrapper(Data1, 
                                   Data2, 
                                   q = semimetric.params$q,
                                   EigenVec = semimetric.params$EigenVec)
  } else if (semimetric == "SemimetricPLS") {
    Dist <- SemimetricPLS(Data1, 
                           Data2, 
                           q = semimetric.params$q)
  } else {
    stop("Unknown semimetric type")
  }
  return(Dist)
}

#' Semimetric Based on Derivatives
#'
#' Computes semimetric distance based on L2-norm of qth derivatives
#'
#' @param Data1 Matrix of functional data (first dataset)
#' @param Data2 Matrix of functional data (second dataset)
#' @param q Integer order of derivative (default: 2)
#' @param nknot Integer number of knots for B-spline (default: 20)
#' @param range.grid Numeric vector of length 2 specifying grid range (default: c(0,1))
#' @param Hhalf Matrix H-half from previous computation (default: NULL)
#' @return List containing semimetric distance matrix and Hhalf matrix
#' @export
semimetric_deriv_wrapper <- function(Data1, Data2, 
                            q = 2, nknot = 20, range.grid = c(0, 1), 
                             Hhalf = NULL) {
    
    if (is.vector(Data1)) Data1 <- as.matrix(t(Data1))
    if (is.vector(Data2)) Data2 <- as.matrix(t(Data2))
  testfordim <- sum(dim(Data1) == dim(Data2)) == 2
  twodatasets <- TRUE
  if (testfordim) twodatasets <- sum(Data1 == Data2) != prod(dim(Data1))
  
  # B-spline approximation of the curves
  p <- ncol(Data1)
  a <- range.grid[1]
  b <- range.grid[2]
  x <- seq(a, b, length = p)
  order.Bspline <- q + 3
  nknotmax <- (p - order.Bspline - 1) %/% 2
  
  if (nknot > nknotmax) {
    stop(paste("Give a number nknot smaller than", 
               nknotmax, "for avoiding ill-conditioned matrix"))
  }
  
  Knot <- seq(a, b, length = nknot + 2)[-c(1, nknot + 2)]
  delta <- sort(c(rep(c(a, b), order.Bspline), Knot))
  Bspline <- splines::splineDesign(delta, x, order.Bspline)
  
  # Numerical integration by the Gauss method 
  point.gauss <- c(-0.9324695142, -0.6612093865, -0.2386191861, 
                   0.2386191861, 0.6612093865, 0.9324695142)
  weight.gauss <- c(0.1713244924, 0.360761573, 0.4679139346, 
                    0.4679139346, 0.360761573, 0.1713244924)
  x.gauss <- 0.5 * ((b + a) + (b - a) * point.gauss)
  lx.gauss <- length(x.gauss)
  
  Bspline.deriv <- splines::splineDesign(delta, 
                                         x.gauss, 
                                         order.Bspline, 
                                         rep(q, lx.gauss))
  
  if (is.null(Hhalf)) {
    Hhalf <- SemimetricDerivDesign(Bspline.deriv, 
                    weight.gauss, 
                                   b - a)
  }
  
  semimetric <- SemimetricDeriv(Hhalf, 
                                Bspline, 
                                Data1, 
                                Data2,
                                twodatasets)
  
    return(list(semimetric = semimetric, Hhalf = Hhalf))
}

#' Semimetric Based on PCA
#'
#' Computes semimetric distance based on principal component analysis
#'
#' @param Data1 Matrix of functional data (first dataset)
#' @param Data2 Matrix of functional data (second dataset)
#' @param q Integer number of principal components
#' @param EigenVec Matrix of eigenvectors from previous computation (default: NULL)
#' @return List containing semimetric distance matrix and eigenvector matrix
#' @export
semimetric_pca_wrapper <- function(Data1, Data2, q, EigenVec = NULL) {
  
  if (is.vector(Data1)) Data1 <- as.matrix(t(Data1))
  if (is.vector(Data2)) Data2 <- as.matrix(t(Data2))
  testfordim <- sum(dim(Data1) == dim(Data2)) == 2
  twodatasets <- TRUE
  if (testfordim) twodatasets <- sum(Data1 == Data2) != prod(dim(Data1))
  
	qmax <- ncol(Data1)
  if (q > qmax) stop(paste("Give an integer q smaller than", qmax))
  
  if (is.null(EigenVec)) {
    EigenVec <- SemimetricPCAEV(Data1, q)
  }
  
  semimetric <- SemimetricPCA(EigenVec, Data1, Data2, twodatasets)
  
  return(list(semimetric = semimetric, EigenVec = EigenVec))
}

#' Semimetric Based on PLS
#'
#' Computes semimetric distance based on partial least squares
#'
#' @param Response Numeric vector of responses
#' @param Data1 Matrix of functional data (first dataset)
#' @param Data2 Matrix of functional data (second dataset)
#' @param q Integer number of PLS components
#' @return Matrix of semimetric distances
#' @export
SemimetricPLS <- function(Response, Data1, Data2, q) {
  
  if (!requireNamespace("pls", quietly = TRUE)) {
    stop("Package 'pls' is required for SemimetricPLS")
  }
  
  if (is.vector(Data1)) Data1 <- as.matrix(t(Data1))
  if (is.vector(Data2)) Data2 <- as.matrix(t(Data2))
  testfordim <- sum(dim(Data1) == dim(Data2)) == 2
  twodatasets <- TRUE
  if (testfordim) twodatasets <- sum(Data1 == Data2) != prod(dim(Data1))
  
  qmax <- ncol(Data1)
  if (q > qmax) stop(paste("Give an integer q smaller than", qmax))
  
  n <- nrow(Data1)
  m <- ncol(Data1)
  mplsr.res <- pls::plsr(Response ~ Data1, ncomp = q)
  Coef <- matrix(mplsr.res$coefficients, m, n)
  b0 <- mplsr.res$Ymeans - t(Coef) %*% mplsr.res$Xmeans
  
  Comp1 <- Data1 %*% Coef
  Comp1 <- outer(rep(1, n), as.vector(b0)) + Comp1
  
  if (twodatasets) {
    n2 <- nrow(Data2)
    Comp2 <- Data2 %*% Coef
    Comp2 <- outer(rep(1, n2), as.vector(b0)) + Comp2
  } else {
    Comp2 <- Comp1
  }
  
  Semimetric <- 0
  for (g in 1:q) {
    Semimetric <- Semimetric + outer(Comp1[, g], Comp2[, g], "-")^2
  }
  
  return(sqrt(Semimetric))
}
