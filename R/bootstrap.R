################################################################################
#
# Get bootstrap data
#
################################################################################
BootstrapData <- function (Y.learn, 
                           Y.learn.hat, 
                           Resampling.Method = "homoscedatic", 
                           NB = 100) {
  Residuals <- Y.learn - Y.learn.hat
  n <- length(Residuals)
  if ( Resampling.Method == "wild.continuous" ) {
    a <- 0.8535 # sqrt(1 - 2 * 20^(-2/3))
    b <- 0.3684 # 20^(-1/3)
    W <- matrix(.Internal(rnorm(n * NB, 0, 1)), n, NB)
    return((a * W + b * (W * W - 1)) * Residuals + matrix(Y.learn.hat, n, NB))
  }
  else if ( Resampling.Method == "wild.twopoint" ) {
    W <- BootstrapKStest(Residuals)
    return(W * Residuals)
  }
  else if ( Resampling.Method == "homoscedatic" ) {
    eps <- Residuals - .Internal(mean(Residuals))
    return(matrix(sample(eps, 
                         n * NB, 
                         replace = T), 
           n, NB) + matrix(Y.learn.hat, n, NB))
  }
  else {
    print ("Choose one of the following resampling methods: 
           wild.continuous, wild.twopoint or homoscedatic!")
  }
}
################################################################################
#
# Calculate the optimal parameter of the the two point distribution 
#
################################################################################
BootstrapKStest <- function (Residuals, 
                             NB = 100, 
                             grid = 10) {
  a.seq <- seq(1, 0.5 * (1 + sqrt(5)), length = grid)
	n <- length(Residuals)
	ord <- order(Residuals)
	Residuals.ordered <- Residuals[ord]
	t <- seq(0,1,length=n)
	KN <- c()
	for (i in 1:grid) {
		W <- matrix(sample(c(a.seq[i], -1 / a.seq[i]), 
                       n * NB, 
                       replace = T, 
                       prob = c(1 / (1 + a.seq[i]^2), 
                                a.seq[i]^2 / (1 + a.seq[i]^2))
                       ), 
                n, NB)
		W <- Residuals * W
    G <- rep(0, n)
		for(j in 1:n) {
			G[j] <- sum(W < Residuals.ordered[j]) / (n * NB)
		} 
    k <- ks.test(G, "punif", 0, 1) # 1 / sqrt(n) * max(abs(G - t / n))
    KN[i] <- k$p.value
	}
	k <- order(KN)[1]
	return (Residuals * matrix(sample(c(a.seq[k], -1 / a.seq[k]), 
                                    n * NB, replace = T, 
                                    prob = c(1 / (1 + a.seq[k]^2), 
                                             a.seq[k]^2 / (1 + a.seq[k]^2))
                                    ), 
                             n, NB))
}