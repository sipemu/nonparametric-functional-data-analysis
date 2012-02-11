\name{nfda-package}
\alias{nfda-package}
\alias{nfda}
\docType{package}
\title{
Nonparametric Functional Data Analysis in R
}
\description{
C++-functions for different nonparametric kernel estimates. Originally written 
in R by Ferraty and Vieu (2006) Nonparametric Functional Data Analysis. Original 
code can be downloaded on their website 
http://www.math.univ-toulouse.fr/staph/npfda/
}
\details{
\tabular{ll}{
Package: \tab nfda\cr
Type: \tab Package\cr
Version: \tab 0.2-1\cr
Date: \tab 2011-12-15\cr
License: \tab GPL-2\cr
LazyLoad: \tab yes\cr
}
}
\author{
Maintainer:
Simon Mueller <Simon.Mueller@mathematik.uni-stuttgart.de>
}
\references{
  http://www.math.univ-toulouse.fr/staph/npfda/

  Ferraty, F. and Vieu, P. \emph{Nonparametric Functional Data Analysis}. 
  Springer 2006.
  
  Rachdi, M. and Vieu, P. \emph{Nonparametric regression for functional data: 
  automatic smoothing parameter selection}. Journal of Statistical Planning 
  and Inference 137, 9 (2007), 2784-2801.
  
  Benhenni, K., Ferraty, F., Rachdi, M., and Vieu, P. \emph{Local smoothing 
  regression with functional data}. Computational Statistics 22, 3 (2007) 353-369.
}
\examples{
##
#
# an example on nonparametric kernel regression
#
##
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
