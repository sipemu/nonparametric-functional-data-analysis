#!/usr/bin/env Rscript
#' ---
#' title: "Fat Content Prediction using Nonparametric Regression"
#' author: "nfda package example"
#' date: "`r Sys.Date()`"
#' ---
#' 
#' This example demonstrates how to use the nfda package for nonparametric
#' functional regression using the Fat spectrometric dataset from the fds package.
#' 
#' The dataset contains near-infrared absorbance spectra (functional data) and
#' corresponding fat content measurements (response variable) for meat samples.

# Load required packages
library(nfda)
library(fds)

cat("====================================\n")
cat("Fat Content Prediction Example\n")
cat("====================================\n\n")

# Load the Fat dataset
data(Fatspectrum, package = "fds")
data(Fatvalues, package = "fds")

cat("Dataset information:\n")
cat("  - Number of curves:", nrow(Fatspectrum$y), "\n")
cat("  - Number of observations:", ncol(Fatspectrum$y), "\n")
cat("  - Fat content range:", range(Fatvalues), "\n\n")

# Prepare data: transpose to get samples x features format
X <- t(Fatspectrum$y)
Y <- Fatvalues

cat("Prepared data dimensions:\n")
cat("  - X (functional data):", dim(X), "\n")
cat("  - Y (response):", length(Y), "\n\n")

# Split into training (70%) and test (30%) sets
set.seed(2025)
n <- nrow(X)
train_size <- floor(0.7 * n)
train_idx <- sample(1:n, size = train_size)
test_idx <- setdiff(1:n, train_idx)

X_train <- X[train_idx, ]
Y_train <- Y[train_idx]
X_test <- X[test_idx, ]
Y_test <- Y[test_idx]

cat("Training set size:", length(Y_train), "\n")
cat("Test set size:", length(Y_test), "\n\n")

# ============================================================================
# Example 1: Derivative-based semimetric with k-NN global cross-validation
# ============================================================================
cat("Example 1: Derivative-based semimetric\n")
cat("---------------------------------------\n")

params_deriv <- list(
  q = 2,              # Second derivative
  nknot = 10,         # Number of B-spline knots
  range.grid = c(0, 1) # Grid range
)

cat("Fitting model with kNN global CV...\n")
model_deriv <- FuNopaRe(X_train, Y_train, 
                        semimetric = "Deriv", 
                        semimetric.params = params_deriv,
                        bandwidth = "kNNgCV")

cat("  Optimal k:", model_deriv$k.opt, "\n")
cat("  Training MSE:", round(model_deriv$mse.learn, 4), "\n\n")

# Make predictions on test set
cat("Making predictions on test set...\n")
pred_deriv <- predict(model_deriv, X_test)

# Calculate test MSE and R-squared
test_mse <- mean((pred_deriv$Prediction - Y_test)^2)
test_mae <- mean(abs(pred_deriv$Prediction - Y_test))
ss_tot <- sum((Y_test - mean(Y_test))^2)
ss_res <- sum((Y_test - pred_deriv$Prediction)^2)
r_squared <- 1 - ss_res / ss_tot

cat("Test set performance:\n")
cat("  MSE:", round(test_mse, 4), "\n")
cat("  MAE:", round(test_mae, 4), "\n")
cat("  R-squared:", round(r_squared, 4), "\n\n")

# ============================================================================
# Example 2: PCA-based semimetric
# ============================================================================
cat("Example 2: PCA-based semimetric\n")
cat("---------------------------------------\n")

params_pca <- list(
  q = 15,             # Number of principal components
  EigenVec = NULL     # Will be computed automatically
)

cat("Fitting model with PCA semimetric...\n")
model_pca <- FuNopaRe(X_train, Y_train, 
                      semimetric = "PCA", 
                      semimetric.params = params_pca,
                      bandwidth = "kNNgCV")

cat("  Optimal k:", model_pca$k.opt, "\n")
cat("  Training MSE:", round(model_pca$mse.learn, 4), "\n\n")

# Make predictions
pred_pca <- predict(model_pca, X_test)

test_mse_pca <- mean((pred_pca$Prediction - Y_test)^2)
test_mae_pca <- mean(abs(pred_pca$Prediction - Y_test))
ss_res_pca <- sum((Y_test - pred_pca$Prediction)^2)
r_squared_pca <- 1 - ss_res_pca / ss_tot

cat("Test set performance:\n")
cat("  MSE:", round(test_mse_pca, 4), "\n")
cat("  MAE:", round(test_mae_pca, 4), "\n")
cat("  R-squared:", round(r_squared_pca, 4), "\n\n")

# ============================================================================
# Example 3: Local cross-validation
# ============================================================================
cat("Example 3: Local cross-validation (adaptive bandwidth)\n")
cat("-------------------------------------------------------\n")

cat("Fitting model with kNN local CV...\n")
model_local <- FuNopaRe(X_train, Y_train, 
                        semimetric = "Deriv", 
                        semimetric.params = params_deriv,
                        bandwidth = "kNNlCV")

cat("  Number of different k values used:", length(unique(model_local$k.opt)), "\n")
cat("  k range: [", min(model_local$k.opt), ",", max(model_local$k.opt), "]\n")
cat("  Training MSE:", round(model_local$mse.learn, 4), "\n\n")

# Make predictions
pred_local <- predict(model_local, X_test)

test_mse_local <- mean((pred_local$Prediction - Y_test)^2)
r_squared_local <- 1 - sum((Y_test - pred_local$Prediction)^2) / ss_tot

cat("Test set performance:\n")
cat("  MSE:", round(test_mse_local, 4), "\n")
cat("  R-squared:", round(r_squared_local, 4), "\n\n")

# ============================================================================
# Summary and Comparison
# ============================================================================
cat("====================================\n")
cat("Summary of Methods\n")
cat("====================================\n\n")

results <- data.frame(
  Method = c("Deriv + kNN Global", "PCA + kNN Global", "Deriv + kNN Local"),
  Test_MSE = c(test_mse, test_mse_pca, test_mse_local),
  R_squared = c(r_squared, r_squared_pca, r_squared_local)
)

print(results, row.names = FALSE)

cat("\nBest method:", results$Method[which.min(results$Test_MSE)], "\n")
cat("\nExample completed successfully!\n")

