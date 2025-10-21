#' Generate Bootstrap Data
#'
#' Generates bootstrap samples using various resampling methods
#'
#' @param Y Numeric vector of learning responses
#' @param Yhat Numeric vector of fitted values for learning data
#' @param Resampling.Method Character string specifying resampling method:
#'   "homoscedatic", "wild.continuous", or "wild.twopoint"
#' @param NB Integer number of bootstrap samples (default: 100)
#' @return Matrix of bootstrap samples (n x NB)
#' @export
#' @examples
#' \dontrun{
#'   Y <- rnorm(100)
#'   Yhat <- rnorm(100)
#'   W <- BootstrapData(Y, Yhat, "homoscedatic", 100)
#' }
BootstrapData <- function(Y, Yhat, 
                          bootstrap_resampling_method = c("homoscedatic", "wild.continuous", "wild.twopoint"),
                          bootstrap_samples = 100) {
  bootstrap_resampling_method <- match.arg(bootstrap_resampling_method, c("homoscedatic", "wild.continuous", "wild.twopoint"))
  
  residuals <- Y - Yhat
  n <- length(residuals)
  
  if (bootstrap_resampling_method == "wild.continuous") {
    # Wild bootstrap with continuous distribution
    a <- 0.8535  # sqrt(1 - 2 * 20^(-2/3))
    b <- 0.3684  # 20^(-1/3)
    W <- matrix(rnorm(n * bootstrap_samples, 0, 1), n, bootstrap_samples)
    return((a * W + b * (W * W - 1)) * residuals + matrix(Yhat, n, bootstrap_samples))
    
  } else if (bootstrap_resampling_method == "wild.twopoint") {
    # Wild bootstrap with two-point distribution
    W <- BootstrapKStest(residuals, bootstrap_samples)
    return(W * residuals)
    
  } else if (bootstrap_resampling_method == "homoscedatic") {
    # Homoscedatic bootstrap
    eps <- residuals - mean(residuals)
    return(matrix(sample(eps, n * bootstrap_samples, replace = TRUE), n, bootstrap_samples) + 
           matrix(Yhat, n, bootstrap_samples))
    
  }
}

#' Calculate Optimal Two-Point Distribution Parameter
#'
#' Uses Kolmogorov-Smirnov test to find optimal parameter for two-point distribution
#'
#' @param residuals Numeric vector of residuals
#' @param bootstrap_samples Integer number of bootstrap samples (default: 100)
#' @param n_grid_points Integer number of grid points for parameter search (default: 10)
#' @return Matrix of bootstrap samples (n x bootstrap_samples)
#' @export
#' @examples
#' \dontrun{
#'   residuals <- rnorm(100)
#'   W <- BootstrapKStest(residuals, bootstrap_samples = 100, n_grid_points = 10)
#' }
BootstrapKStest <- function(residuals, bootstrap_samples = 100, n_grid_points = 10) {
  
  a_seq <- seq(1, 0.5 * (1 + sqrt(5)), length = n_grid_points)
  n <- length(residuals)
  ord <- order(residuals)
  residuals_ordered <- residuals[ord]
  t <- seq(0, 1, length = n)
  KS_test_p_values <- numeric(n_grid_points)
  
  for (i in 1:n_grid_points) {
    # Generate two-point distribution samples
    W <- matrix(
      sample(
        c(a_seq[i], -1 / a_seq[i]), 
        size = n * bootstrap_samples, 
        replace = TRUE, 
        prob = c(1 / (1 + a_seq[i]^2), a_seq[i]^2 / (1 + a_seq[i]^2))
      ), 
      nrow = n, 
      ncol = bootstrap_samples
    )
    W <- residuals * W
    
    # Calculate empirical distribution
    G <- numeric(n)
    for (j in 1:n) {
      G[j] <- sum(W < residuals_ordered[j]) / (n * bootstrap_samples)
    }
    
    # Kolmogorov-Smirnov test
    k <- ks.test(G, "punif", 0, 1)
    KS_test_p_values[i] <- k$p.value
  }
  
  # Select best parameter
  best_parameter_index <- order(KS_test_p_values)[1]
  
  # Generate final bootstrap samples with optimal parameter
  result <- residuals * matrix(
    sample(
      c(a_seq[best_parameter_index], -1 / a_seq[best_parameter_index]), 
      size = n * bootstrap_samples, 
      replace = TRUE, 
      prob = c(1 / (1 + a_seq[best_parameter_index]^2), a_seq[best_parameter_index]^2 / (1 + a_seq[best_parameter_index]^2))
    ), 
    nrow = n, 
    ncol = bootstrap_samples
  )
  
  return(result)
}
