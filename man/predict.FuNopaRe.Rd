\name{predict.FuNopaRe}
\alias{predict.FuNopaRe}
\title{Nonparametric Regression for Functional Data}
\description{
  \code{predict.FuNopaRe} is a function for predictions from the result of
  nonparametric modell fitting.
}
\usage{
\method{predict}{FuNopaRe}(object, newdata, method.params, Bootstrapping = FALSE, ...)
}
\arguments{
  \item{object}{ Afitted object of class inheriting from \code{FuNopaRe}}
  \item{newdata}{ Matrix with the functional data (curves) each row one for prediction}
  \item{method.params}{ Parameters for bootstrapping}
  \item{Bootstrapping}{ Using bootstrapping for local adaptive bandwidth
                        selection}
  \item{...}{further arguments passed to or from other methods}
}
\value{
  \code{FuNopaRe} returns an object of the \code{class} \code{FuNopaRe}; 
  additional with the predictions in the vector \code{Prediction}.
}
\references{ 
  Ferraty, F. and Vieu, P. \emph{Nonparametric Functional Data Analysis}. 
  Springer 2006.
  
  Rachdi, M. and Vieu, P. \emph{Nonparametric regression for functional data: 
  automatic smoothing parameter selection}. Journal of Statistical Planning 
  and Inference 137, 9 (2007), 2784-2801.
  
  Benhenni, K., Ferraty, F., Rachdi, M., and Vieu, P. \emph{Local smoothing 
  regression with functional data}. Computational Statistics 22, 3 (2007) 353???369.
 }
\author{Simon Mueller \email{simon.mueller@mathematik.uni-stuttgart.de}}
\seealso{\link{Semimetric}}
\examples{
# functional data sets
library (fds)

# fat spectrum dataset
Y <- Fatvalues
X <- t(Fatspectrum$y)

# setup semimetric parameters
semimetric.params <- c()
semimetric.params$q <- 2
semimetric.params$nknot  <- 20
semimetric.params$range.grid  <- c (min (Fatspectrum$x), 
                                    max (Fatspectrum$x))

# learn and testsample
learn <- 1:160
test <- 161:215

# parameter estimation and prediction by cross-validation
Learn.Fat.CV <- FuNopaRe (X[learn, ], 
                          Y[learn], 
                          semimetric = "Deriv", 
                          semimetric.params, 
                          bandwidth = "CV")
                          
Predict.Fat.CV <- predict (Learn.Fat.CV, 
                           X[test, ], 
                           method.params = NULL)

plot (Predict.Fat.CV$Prediction, Y[161:215])                          
}
