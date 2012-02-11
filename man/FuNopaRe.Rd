\name{FuNopaRe}
\alias{FuNopaRe}
\title{Nonparametric Regression for Functional Data}
\description{
  \code{FuNopaRe} is a function that estimates optimal bandwidth of the kernel 
  estimate based on the learning data.
}
\usage{
FuNopaRe(X, Y, semimetric, semimetric.params, bandwidth = "CV")
}
\arguments{
  \item{X}{ Matrix with the functional data (curves) each row one}
  \item{Y}{ Vector of the scalar responses}
  \item{semimetric}{ A string of choosing the semimetric; allowed are:
                     "Deriv" and "PCA"}
  \item{semimetric.params}{ Parameters for the semimetric function.}
  \item{bandwidth}{ Method for choosing the bandwidth; allowed are:
                    "CV"(default), "kNNgCV", and "kNNlCV"}
}
\value{
  \code{FuNopaRe} returns an object of the \code{class} \code{FuNopaRe};
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
\seealso{\link{Semimetric}, \link{predict.FuNopaRe}}
\examples{
# functional datasets
require (fds)

# fat spectrum dataset
Y <- Fatvalues
X <- t(Fatspectrum$y)

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
Learn.Fat.CV <- FuNopaRe (X[learn, ], 
                          Y[learn], 
                          semimetric = "Deriv", 
                          semimetric.params, 
                          bandwidth = "CV")
}
