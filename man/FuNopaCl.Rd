\name{FuNopaCl}
\alias{FuNopaCl}
\title{Nonparametric Classification for Functional Data}
\description{
  \code{FuNopaCl} is a function that estimates optimal bandwidth by 
  k nearest neighbour local cross-validation for the kernel estimate 
  based on the learning data.
}
\usage{
FuNopaCl(X, classes, semimetric, semimetric.params)
}
\arguments{
  \item{X}{ Matrix with the functional data (curves) each row one}
  \item{classes}{ Vector of the classes}
  \item{semimetric}{ A string of choosing the semimetric; allowed are:
                     "Deriv" and "PCA"}
  \item{semimetric.params}{ Parameters for the semimetric function.}
}
\value{
  \code{FuNopaCl} returns an object of the \code{class} \code{FuNopaCl};
}
\references{ 
  Ferraty, F. and Vieu, P. \emph{Nonparametric Functional Data Analysis}. 
  Springer 2006.
 }
\author{Simon Mueller \email{simon.mueller@mathematik.uni-stuttgart.de}}
\seealso{\link{Semimetric}, \link{predict.FuNopaCl}}
\examples{
# functional datasets
require (fds)

# fat spectrum dataset
Y <- Fatvalues
X <- t(Fatspectrum$y)
Y[Y < 20] <- 1
Y[Y >= 20] <- 2

# setup semimetric params
semimetric.params <- c()
semimetric.params$q <- 2
semimetric.params$nknot  <- 20
semimetric.params$range.grid  <- c (min (Fatspectrum$x), 
                                    max (Fatspectrum$x))
# learn and testsample
learn <- 1:160
test <- 161:215

# parameter estimation and prediction by cross-validation
Learn.Fat.CV <- FuNopaCl (X[learn, ], 
                          Y[learn], 
                          semimetric = "Deriv", 
                          semimetric.params)
}
