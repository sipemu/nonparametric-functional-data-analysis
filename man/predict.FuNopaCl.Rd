\name{predict.FuNopaCl}
\alias{predict.FuNopaCl}
\title{Nonparametric Classification for Functional Data}
\description{
  \code{predict.FuNopaCl} is a function for predictions from the result of
  nonparametric modell fitting.
}
\usage{
\method{predict}{FuNopaCl}(object, newdata,  ...)
}
\arguments{
  \item{object}{ Afitted object of class inheriting from \code{FuNopaCl}}
  \item{newdata}{ Matrix with the functional data (curves) each row one for prediction}
  \item{...}{further arguments passed to or from other methods}
}
\value{
  \code{FuNopaCl} returns an object of the \code{class} \code{FuNopaCl}; 
  additional with the predictions of the classes in the vector \code{Prediction}.
}
\references{ 
  Ferraty, F. and Vieu, P. \emph{Nonparametric Functional Data Analysis}. 
  Springer 2006.
 }
\author{Simon Mueller \email{simon.mueller@mathematik.uni-stuttgart.de}}
\seealso{\link{Semimetric}}
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
Learn.Fat <- FuNopaCl (X[learn, ], 
                       Y[learn], 
                       semimetric = "Deriv", 
                       semimetric.params)
                          
Predict.Fat <- predict (Learn.Fat, 
                        X[test, ])

mean (Predict.Fat$classes.pred != Y[161:215])                          
}
