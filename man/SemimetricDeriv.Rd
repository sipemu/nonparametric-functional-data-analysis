\name{SemimetricDeriv}
\alias{SemimetricDeriv}
\title{SemimetricDeriv}
\description{
  \code{SemimetricDeriv} is a function that calculates a distance matrix based 
  on the L_2-norm of the derivatives of the functional data. A Rcpp port of the 
  semimetric function based on derivatives of Ferraty and Vieu.
}
\usage{
SemimetricDeriv(Data1, 
                Data2, 
                q = 2, nknot = 20, range.grid = c(0,1), 
                Hhalf = NULL)
}
\arguments{
  \item{Data1}{ Matrix with functional data (curves) each row one; n rows}
  \item{Data2}{ Matrix with a second set of functional data (curves) 
                each row one; m rows}
  \item{q}{ Order of derivation}
  \item{nknot}{ Number of the interior knots for the B-spline basis}
  \item{range.grid}{ Vector of length two, containing the domain of the curves}
  \item{Hhalf}{ A design matrix that is based on the B-spline basis design}
}
\value{
  \code{Semimetric} returns a list with: 
                    \code{semimetric} as the n x m - distance-matrix 
                    \code{Hhalf} A design matrix that is based on the B-spline 
                    basis design, which can be reused in following semimetric
                    calculations
}
\references{ 
  http://www.math.univ-toulouse.fr/staph/npfda/

  Ferraty, F. and Vieu, P. \emph{Nonparametric Functional Data Analysis}. 
  Springer 2006.
 }
\author{Simon Mueller \email{simon.mueller@mathematik.uni-stuttgart.de}}
\seealso{\link{Semimetric}}
\examples{
# functional data sets
require (fds)

# fat spectrum dataset
X <- t(Fatspectrum$y)

# setup semimetric parameters
q <- 2
nknot <- 20
range.grid <- c (min (Fatspectrum$x), 
                 max (Fatspectrum$x))

# calculate distance matrix
Dist1 <- SemimetricDeriv (X[1:100, ], 
                          X[101:150, ], 
                          q, nknot, range.grid)

# calculate distance matrix using Hhalf
Dist2 <- SemimetricDeriv (X[1:100, ], 
                          X[151:215, ], 
                          q, nknot, range.grid, 
                          Dist1$Hhalf)
}
