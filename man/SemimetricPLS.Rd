\name{SemimetricPLS}
\alias{SemimetricPLS}
\title{SemimetricPLS}
\description{
  \code{SemimetricPLS} is a function that calculates a distance matrix based 
  on the PLS components of the functional data. Response need not to be classes
  as by Ferraty and Vieu.
}
\usage{
SemimetricPLS(Response,
              Data1, 
              Data2, 
              q)
}
\arguments{
  \item{Response}{The response of Data1}
  \item{Data1}{ Matrix with functional data (curves) each row one; n rows}
  \item{Data2}{ Matrix with a second set of functional data (curves) 
                each row one; m rows}
  \item{q}{ Number of pls components}
}
\value{
  \code{Semimetric} returns a list with: 
                    \code{semimetric} as the n x m - distance-matrix 
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
Y <- Fatvalues

# setup semimetric parameter
q <- 4

# calculate distance matrix
Dist1 <- SemimetricPLS (X[1:100],
                        X[1:100, ], 
                        X[101:150, ], 
                        q)
}
