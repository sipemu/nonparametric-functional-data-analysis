\name{BootstrapKStest}
\alias{BootstrapKStest}
\title{BootstrapKStest}
\description{
  \code{BootstrapKStest} is a function that calculates the optimal parameter of
  the two-point distribution based on the Kolmogorov-Smirnov-Test.
}
\usage{
BootstrapKStest(Residuals, NB = 100, grid = 10)
}
\arguments{
  \item{Residuals}{ Y.learn - Y.lern.hat}
  \item{NB}{ Number of bootstrap samples}
  \item{grid}{ Size of the grid}
}
\value{
  \code{W} n x NB-Matrix of bootstrapped residuals.
}
\references{ 
  Davidson, J., Monticini, A., and Peel, D. \emph{Implementing the wild 
  bootstrap using a two-point distribution}. Economics Letters, Vol. 96, 
  No. 3. (2007), pp. 309-315.
 }
\author{Simon Mueller \email{simon.mueller@mathematik.uni-stuttgart.de}}