\name{Semimetric}
\alias{Semimetric}
\title{Semimetric}
\description{
  \code{Semimetric} is a function that calculates a distance matrix for 
  functional data. A Rcpp port of the semimetric functions of Ferraty and Vieu.
}
\usage{
Semimetric(Data1, Data2, semimetric, semimetric.params)
}
\arguments{
  \item{Data1}{ Matrix with functional data (curves) each row one; n rows}
  \item{Data2}{ Matrix with a second set of functional data (curves) 
                each row one; m rows}
  \item{semimetric}{ A string of choosing the semimetric; allowed are:
                     "Deriv", "PCA", and "PLS"}
  \item{semimetric.params}{ Parameters for the semimetric function.}
}
\value{
  \code{Semimetric} returns a list with: 
                    \code{semimetric} as the n x m - distance-matrix 
                    \code{...} reusable calculations, 
                    e.g. PCA: the q eigen vectors 
}
\references{ 
  http://www.math.univ-toulouse.fr/staph/npfda/

  Ferraty, F. and Vieu, P. \emph{Nonparametric Functional Data Analysis}. 
  Springer 2006.
 }
\author{Simon Mueller \email{simon.mueller@mathematik.uni-stuttgart.de}}
\seealso{\link{SemimetricDeriv}, \link{SemimetricPCA}, \link{SemimetricPLS}}
\examples{
# functional data sets
require (fds)

# fat spectrum dataset
X <- t(Fatspectrum$y)

# setup semimetric parameters
semimetric.params <- c ()
semimetric.params$q <- 2
semimetric.params$nknot <- 20
semimetric.params$range.grid <- c (min (Fatspectrum$x), 
                                   max (Fatspectrum$x))

# calculate distance matrix
Dist <- Semimetric (X, 
                    X, 
                    semimetric = "SemimetricDeriv", 
                    semimetric.params)
}
