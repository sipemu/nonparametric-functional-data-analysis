#' Generate Bootstrap Data
#'
#' Generates bootstrap samples using various resampling methods
#'
#' @param Y.learn Numeric vector of learning responses
#' @param Y.learn.hat Numeric vector of fitted values for learning data
#' @param Resampling.Method Character string specifying resampling method:
#'   "homoscedatic", "wild.continuous", or "wild.twopoint"
#' @param NB Integer number of bootstrap samples (default: 100)
#' @return Matrix of bootstrap samples (n x NB)
#' @export
#' @examples
#' \dontrun{
#'   Y <- rnorm(100)
#'   Y.hat <- rnorm(100)
#'   W <- BootstrapData(Y, Y.hat, "homoscedatic", 100)
#' }
BootstrapData <- function(Y.learn, Y.learn.hat, 
                          Resampling.Method = "homoscedatic", NB = 100) {
  
  Residuals <- Y.learn - Y.learn.hat
  n <- length(Residuals)
  
  if (Resampling.Method == "wild.continuous") {
    # Wild bootstrap with continuous distribution
    a <- 0.8535  # sqrt(1 - 2 * 20^(-2/3))
    b <- 0.3684  # 20^(-1/3)
    W <- matrix(rnorm(n * NB, 0, 1), n, NB)
    return((a * W + b * (W * W - 1)) * Residuals + matrix(Y.learn.hat, n, NB))
    
  } else if (Resampling.Method == "wild.twopoint") {
    # Wild bootstrap with two-point distribution
    W <- BootstrapKStest(Residuals, NB)
    return(W * Residuals)
    
  } else if (Resampling.Method == "homoscedatic") {
    # Homoscedatic bootstrap
    eps <- Residuals - mean(Residuals)
    return(matrix(sample(eps, n * NB, replace = TRUE), n, NB) + 
           matrix(Y.learn.hat, n, NB))
    
  } else {
    stop(paste("Unknown resampling method:", Resampling.Method,
               "\nChoose one of: wild.continuous, wild.twopoint, or homoscedatic"))
  }
}

#' Calculate Optimal Two-Point Distribution Parameter
#'
#' Uses Kolmogorov-Smirnov test to find optimal parameter for two-point distribution
#'
#' @param Residuals Numeric vector of residuals
#' @param NB Integer number of bootstrap samples (default: 100)
#' @param grid Integer number of grid points for parameter search (default: 10)
#' @return Matrix of bootstrap samples (n x NB)
#' @export
#' @examples
#' \dontrun{
#'   Residuals <- rnorm(100)
#'   W <- BootstrapKStest(Residuals, 100, 10)
#' }
BootstrapKStest <- function(Residuals, NB = 100, grid = 10) {
  
  a.seq <- seq(1, 0.5 * (1 + sqrt(5)), length = grid)
  n <- length(Residuals)
  ord <- order(Residuals)
  Residuals.ordered <- Residuals[ord]
  t <- seq(0, 1, length = n)
  KN <- numeric(grid)
  
  for (i in 1:grid) {
    # Generate two-point distribution samples
    W <- matrix(
      data = sample(
        x = c(a.seq[i], -1 / a.seq[i]), 
        size = n * NB, 
        replace = TRUE, 
        prob = c(1 / (1 + a.seq[i]^2), a.seq[i]^2 / (1 + a.seq[i]^2))
      ), 
      nrow = n, 
      ncol = NB
    )
    W <- Residuals * W
    
    # Calculate empirical distribution
    G <- numeric(n)
    for (j in 1:n) {
      G[j] <- sum(W < Residuals.ordered[j]) / (n * NB)
    }
    
    # Kolmogorov-Smirnov test
    k <- ks.test(G, "punif", 0, 1)
    KN[i] <- k$p.value
  }
  
  # Select best parameter
  k <- order(KN)[1]
  
  # Generate final bootstrap samples with optimal parameter
  result <- Residuals * matrix(
    x = sample(
      x = c(a.seq[k], -1 / a.seq[k]), 
      size = n * NB, 
      replace = TRUE, 
      prob = c(1 / (1 + a.seq[k]^2), a.seq[k]^2 / (1 + a.seq[k]^2))
    ), 
    nrow = n, 
    ncol = NB
  )
  
  return(result)
}
