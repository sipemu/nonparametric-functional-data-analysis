\name{BootstrapData}
\alias{BootstrapData}
\title{BootstrapData}
\description{
  \code{BootstrapData} is a function that calculates the bootstrap response 
  based on Yhat
}
\usage{
BootstrapData (Y.learn, 
               Y.learn.hat, 
               Resampling.Method = "homoscedatic", 
               NB = 100)
}
\arguments{
  \item{Y.learn}{ n response values of the learning set }
  \item{Y.learn.hat}{ n prediction of the response values of the learning set}
  \item{Resampling.Method}{ Bootstrap method:
                     "homoscedastic" (default), "wild.continuous", 
                     and "wild.twopoint"}
  \item{NB}{ Number of bootstrap samples}
}
\value{
  \code{W} n x NB-Matrix of bootstrapped response.
}
\references{ 
  Ferraty, F., Van Keillegom, I., and Vieu, P. \emph{On the Validity of the 
  Bootstrap in Non-Parametric Functional Regression}. Scandinavian Journal of 
  Statistics, Vol. 37 (2010), pp. 286-306.
  
  Davidson, J., Monticini, A., and Peel, D. \emph{Implementing the wild 
  bootstrap using a two-point distribution}. Economics Letters, Vol. 96, 
  No. 3. (2007), pp. 309-315.
  
  Mammen, E. \emph{Resampling Methods for Nonparametric Regression}. 
  In Smoothing and Regres- sion: Approaches, Computation, and Application, 
  M. G. Schimek, Ed. John Wiley, 2000.
 }
\author{Simon Mueller \email{simon.mueller@mathematik.uni-stuttgart.de}}
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

method.params <- c()
method.params$NB <- 100
method.params$Resampling.Method <- "homoscedatic"
method.params$neighbours <- 20
method.params$alpha <- 0.05
Predict.Fat.CV <- predict (Learn.Fat.CV, 
                           X[test, ], 
                           method.params = method.params,
                           Bootstrapping = TRUE)

plot (Predict.Fat.CV$Prediction, Y[161:215])                          
}