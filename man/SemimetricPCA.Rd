\name{SemimetricPCA}
\alias{SemimetricPCA}
\title{SemimetricPCA}
\description{
  \code{SemimetricPCA} is a function that calculates a distance matrix based 
  on the PCA of the functional data. A Rcpp port of the 
  PCA-type semimetric function of Ferraty and Vieu.
}
\usage{
SemimetricPCA(Data1, 
              Data2, 
              q, 
              EigenVec = NULL)
}
\arguments{
  \item{Data1}{ Matrix with functional data (curves) each row one; n rows}
  \item{Data2}{ Matrix with a second set of functional data (curves) 
                each row one; m rows}
  \item{q}{ Number of principal components}
  \item{EigenVec}{ The Eigenvectors based on DATA1}
}
\value{
  \code{Semimetric} returns a list with: 
                    \code{semimetric} as the n x m - distance-matrix 
                    \code{EigenVec} The Eigenvectors based on DATA1, which can
                    be reused in following semimetric calculations
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

# setup semimetric parameter
q <- 4

# calculate distance matrix
Dist1 <- SemimetricPCA (X[1:100, ], 
                        X[101:150, ], 
                        q)

# calculate distance matrix using Hhalf
Dist2 <- SemimetricPCA (X[1:100, ], 
                        X[151:215, ], 
                        q, 
                        Dist1$EigenVec)
}
