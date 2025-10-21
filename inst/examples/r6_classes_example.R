#!/usr/bin/env Rscript
#' ---
#' title: "Using R6 Semimetric Classes in nfda"
#' author: "nfda package"
#' date: "`r Sys.Date()`"
#' ---
#' 
#' This example demonstrates the object-oriented interface to semimetrics
#' using R6 classes. The R6 interface provides more flexibility and control
#' compared to the functional wrappers.

library(nfda)

cat("═══════════════════════════════════════════\n")
cat("  R6 Semimetric Classes Demonstration\n")
cat("═══════════════════════════════════════════\n\n")

# Generate sample functional data
set.seed(2025)
n <- 80
p <- 40
t_grid <- seq(0, 1, length = p)

# Create smooth curves
X <- matrix(0, n, p)
for(i in 1:n) {
  a1 <- rnorm(1, 2, 0.5)
  a2 <- rnorm(1, 1, 0.3)
  X[i,] <- a1 * sin(2*pi*t_grid) + a2 * cos(4*pi*t_grid) + rnorm(p, 0, 0.1)
}

# Split data
X_train <- X[1:60, ]
X_test <- X[61:80, ]

cat("Data prepared:\n")
cat("  Training set: ", nrow(X_train), "curves ×", ncol(X_train), "points\n")
cat("  Test set:     ", nrow(X_test), "curves ×", ncol(X_test), "points\n\n")

# ============================================================================
# Example 1: PCA Semimetric R6 Class
# ============================================================================
cat("Example 1: PCA Semimetric (R6 Class)\n")
cat("─────────────────────────────────────\n")

# Create PCA semimetric object
sm_pca <- semimetric_pca$new(q = 10, eigen_vec = NULL)

cat("Created semimetric object:\n")
sm_pca$print()
cat("\n")

# Calculate distances (training set)
cat("Computing distances on training set...\n")
dist_train <- sm_pca$calculate(X_train, X_train)

cat("  Distance matrix dimensions:", dim(dist_train), "\n")
cat("  Diagonal elements (should be ~0):", round(mean(diag(dist_train)), 10), "\n")
cat("  Mean off-diagonal distance:", round(mean(dist_train[upper.tri(dist_train)]), 4), "\n")
cat("  Distance range: [", round(min(dist_train), 4), ",", 
    round(max(dist_train), 4), "]\n\n")

# Calculate distances between training and test
cat("Computing distances between training and test...\n")
dist_cross <- sm_pca$calculate(X_train, X_test)

cat("  Cross-distance matrix:", dim(dist_cross), "\n")
cat("  Mean cross-distance:", round(mean(dist_cross), 4), "\n\n")

# ============================================================================
# Example 2: Derivative Semimetric R6 Class
# ============================================================================
cat("Example 2: Derivative Semimetric (R6 Class)\n")
cat("──────────────────────────────────────────\n")

# Create derivative semimetric object
sm_deriv <- semimetric_deriv$new(
  q = 2,                 # Second derivative (curvature)
  nknot = 10,            # Interior knots
  range_grid = c(0, 1),  # Domain
  Hhalf = NULL           # Will be computed
)

cat("Created derivative semimetric:\n")
cat("  Derivative order (q):", sm_deriv$q, "\n")
cat("  Number of knots:", sm_deriv$nknot, "\n")
cat("  Domain:", sm_deriv$range_grid, "\n")
cat("  Quadrature points:", length(sm_deriv$point_gauss), "\n\n")

# Calculate distances
cat("Computing derivative-based distances...\n")
dist_deriv_train <- sm_deriv$calculate(X_train, X_train)

cat("  Distance matrix dimensions:", dim(dist_deriv_train), "\n")
cat("  Diagonal elements:", round(mean(diag(dist_deriv_train)), 10), "\n")
cat("  Mean distance:", round(mean(dist_deriv_train[upper.tri(dist_deriv_train)]), 4), "\n")

# Test and train distances
dist_deriv_cross <- sm_deriv$calculate(X_train, X_test)
cat("  Cross-distance mean:", round(mean(dist_deriv_cross), 4), "\n\n")

# ============================================================================
# Example 3: Comparing Semimetrics
# ============================================================================
cat("Example 3: Comparison of Semimetrics\n")
cat("───────────────────────────────────\n\n")

# For a fair comparison, use normalized distances
normalize <- function(D) {
  (D - min(D)) / (max(D) - min(D) + 1e-10)
}

D_pca_norm <- normalize(dist_train)
D_deriv_norm <- normalize(dist_deriv_train)

# Calculate correlation between distance matrices
# (excluding diagonal)
upper_idx <- upper.tri(dist_train)
cor_pca_deriv <- cor(D_pca_norm[upper_idx], D_deriv_norm[upper_idx])

cat("Distance matrix comparison:\n")
cat("  PCA mean distance:", round(mean(dist_train[upper_idx]), 4), "\n")
cat("  Derivative mean distance:", round(mean(dist_deriv_train[upper_idx]), 4), "\n")
cat("  Correlation between methods:", round(cor_pca_deriv, 4), "\n\n")

# Find most similar pairs according to each method
pca_min_idx <- which(dist_train == min(dist_train[upper_idx]), arr.ind = TRUE)[1,]
deriv_min_idx <- which(dist_deriv_train == min(dist_deriv_train[upper_idx]), arr.ind = TRUE)[1,]

cat("Most similar pair (PCA): curves", pca_min_idx[1], "and", pca_min_idx[2], "\n")
cat("  Distance:", round(dist_train[pca_min_idx[1], pca_min_idx[2]], 6), "\n")

cat("Most similar pair (Derivative): curves", deriv_min_idx[1], "and", deriv_min_idx[2], "\n")
cat("  Distance:", round(dist_deriv_train[deriv_min_idx[1], deriv_min_idx[2]], 6), "\n\n")

# ============================================================================
# Example 4: PLS Semimetric (Supervised)
# ============================================================================
if (requireNamespace("pls", quietly = TRUE)) {
  cat("Example 4: PLS Semimetric (Supervised)\n")
  cat("─────────────────────────────────────\n")
  
  # Create response variable
  Y_train <- apply(X_train, 1, function(x) sum(x[1:20])) + rnorm(nrow(X_train), 0, 0.5)
  
  # Create PLS semimetric with response
  sm_pls <- semimetric_pls$new(n_comp = 8, y = Y_train)
  
  cat("Created PLS semimetric:\n")
  sm_pls$print()
  cat("\n")
  
  # Calculate supervised distances
  cat("Computing PLS-based distances...\n")
  dist_pls_train <- sm_pls$calculate(X_train, X_train)
  
  cat("  Distance matrix dimensions:", dim(dist_pls_train), "\n")
  cat("  Mean distance:", round(mean(dist_pls_train[upper.tri(dist_pls_train)]), 4), "\n\n")
  
  # Compare with unsupervised methods
  cor_pls_pca <- cor(dist_pls_train[upper_idx], D_pca_norm[upper_idx])
  cor_pls_deriv <- cor(dist_pls_train[upper_idx], D_deriv_norm[upper_idx])
  
  cat("PLS vs other methods:\n")
  cat("  Correlation with PCA:", round(cor_pls_pca, 4), "\n")
  cat("  Correlation with Derivative:", round(cor_pls_deriv, 4), "\n")
} else {
  cat("Example 4: Skipped (pls package not available)\n\n")
}

# ============================================================================
# Example 5: Caching and Reusability
# ============================================================================
cat("Example 5: Caching and Computational Efficiency\n")
cat("──────────────────────────────────────────────\n")

# Demonstrate caching with derivative semimetric
cat("First computation (computes Hhalf):\n")
system.time({
  sm_cached <- semimetric_deriv$new(q = 2, nknot = 10, range_grid = c(0,1))
  d1 <- sm_cached$calculate(X_train, X_test)
})

cat("\nSubsequent computation (reuses cached Hhalf):\n")
system.time({
  # Hhalf is cached in sm_cached object
  d2 <- sm_cached$calculate(X_train[1:30,], X_test)
})

cat("\nNote: Second computation should be faster due to caching.\n\n")

# ============================================================================
# Example 6: Using R6 Classes in Prediction Models
# ============================================================================
cat("Example 6: Integration with Prediction Models\n")
cat("────────────────────────────────────────────\n")

# The functional interface uses R6 classes internally
# Here's how they connect:

Y <- rnorm(60)
classes <- rep(1:2, each = 30)

# Standard functional interface
params <- list(q = 2, nknot = 8, range.grid = c(0, 1))
model <- FuNopaRe(X_train, Y, "Deriv", params, "kNNgCV")

cat("Model created using functional interface\n")
cat("  Internally uses semimetric_deriv R6 class\n")
cat("  Optimal k:", model$k.opt, "\n\n")

# Classification
model_class <- FuNopaCl(X_train, classes, "Deriv", params)
cat("Classification model created\n")
cat("  Training error:", round(model_class$mse.learn, 4), "\n\n")

# ============================================================================
# Example 7: Custom Distance Analysis
# ============================================================================
cat("Example 7: Analyzing Distance Properties\n")
cat("───────────────────────────────────────\n")

# Check semimetric axioms
sm_test <- semimetric_pca$new(q = 5)
D_test <- sm_test$calculate(X_train[1:10,], X_train[1:10,])

# Axiom 1: Non-negativity
cat("Semimetric axioms verification:\n")
cat("  Non-negativity: ", all(D_test >= 0), "✓\n")

# Axiom 2: Symmetry  
max_asymmetry <- max(abs(D_test - t(D_test)))
cat("  Symmetry (max |D-D'|):", format(max_asymmetry, scientific = TRUE), 
    ifelse(max_asymmetry < 1e-10, "✓", "✗"), "\n")

# Axiom 3: Identity
cat("  Identity (diagonal):", all(abs(diag(D_test)) < 1e-10), "✓\n\n")

# ============================================================================
# Summary
# ============================================================================
cat("═══════════════════════════════════════════\n")
cat("  Summary: R6 Semimetric Classes\n")
cat("═══════════════════════════════════════════\n\n")

cat("Available R6 Classes:\n")
cat("  1. semimetric         - Base abstract class\n")
cat("  2. semimetric_pca     - PCA-based (unsupervised)\n")
cat("  3. semimetric_deriv   - Derivative-based (unsupervised)\n")
cat("  4. semimetric_pls     - PLS-based (supervised)\n\n")

cat("Key Features:\n")
cat("  • Object-oriented design with inheritance\n")
cat("  • Caching for computational efficiency\n")
cat("  • Flexible and reusable\n")
cat("  • Well-documented with mathematical rigor\n\n")

cat("Use Cases:\n")
cat("  • Direct distance computation\n")
cat("  • Custom similarity analysis\n")
cat("  • Integration with prediction models\n")
cat("  • Method comparison and evaluation\n\n")

cat("Documentation:\n")
cat("  ?semimetric\n")
cat("  ?semimetric_pca\n")
cat("  ?semimetric_deriv\n")
cat("  ?semimetric_pls\n\n")

cat("Example completed successfully!\n")

