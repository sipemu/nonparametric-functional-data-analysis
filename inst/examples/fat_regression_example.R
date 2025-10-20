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
bootstrap.params <- list("Resampling.Method" = "homoscedatic", 
                         "NB" = 500, 
                         "neighbours" = 10, 
                         "alpha" = 0.2)
pred_deriv <- predict(model_deriv, X_test, method.params = bootstrap.params, Bootstrapping = T)

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
# Example 2: PLS-based semimetric
# ============================================================================
cat("Example 2: PLS-based semimetric\n")
cat("---------------------------------------\n")

params_pls <- list(
  q = 25, 
  y = Y_train
)

cat("Fitting model with PLS semimetric...\n")
model_pls <- FuNopaRe(X_train, Y_train, 
                      semimetric = "PLS", 
                      semimetric.params = params_pls,
                      bandwidth = "kNNgCV")

cat("  Optimal k:", model_pls$k.opt, "\n")
cat("  Training MSE:", round(model_pls$mse.learn, 4), "\n\n")

# Make predictions
pred_pls <- predict(model_pls, X_test)

test_mse_pls <- mean((pred_pls$Prediction - Y_test)^2)
test_mae_pls <- mean(abs(pred_pls$Prediction - Y_test))
ss_res_pls <- sum((Y_test - pred_pls$Prediction)^2)
r_squared_pls <- 1 - ss_res_pls / ss_tot

cat("Test set performance:\n")
cat("  MSE:", round(test_mse_pls, 4), "\n")
cat("  MAE:", round(test_mae_pls, 4), "\n")
cat("  R-squared:", round(r_squared_pls, 4), "\n\n")

# ============================================================================
# Example 3: Local cross-validation
# ============================================================================
cat("Example 3: Local cross-validation (adaptive bandwidth)\n")
cat("-------------------------------------------------------\n")

cat("Fitting model with kNN local CV...\n")
model_local <- FuNopaRe(X_train, Y_train, 
                        semimetric = "PLS", 
                        semimetric.params = params_pls,
                        bandwidth = "kNNlCV")

cat("  Number of different k values used:", length(unique(model_local$k.opt)), "\n")
cat("  k range: [", min(model_local$k.opt), ",", max(model_local$k.opt), "]\n")
cat("  Training MSE:", round(model_local$mse.learn, 4), "\n\n")

# Make predictions
bootstrap.params <- list("Resampling.Method" = "homoscedatic", "NB" = 200, "neighbours" = 20, "alpha" = 0.05)
pred_local <- predict(model_local, X_test, method.params = bootstrap.params, Bootstrapping = T)

test_mse_local <- mean((pred_local$Prediction - Y_test)^2)
test_mae_local <- mean(abs(pred_local$Prediction - Y_test))
r_squared_local <- 1 - sum((Y_test - pred_local$Prediction)^2) / ss_tot

cat("Test set performance:\n")
cat("  MSE:", round(test_mse_local, 4), "\n")
cat("  MAE:", round(test_mae_local, 4), "\n")
cat("  R-squared:", round(r_squared_local, 4), "\n\n")

# ============================================================================
# Summary and Comparison
# ============================================================================
cat("====================================\n")
cat("Summary of Methods\n")
cat("====================================\n\n")

results <- data.frame(
  Method = c("Deriv + kNN Global", "PLS + kNN Global", "Deriv + kNN Local"),
  Test_MSE = c(test_mse, test_mse_pls, test_mse_local),
  R_squared = c(r_squared, r_squared_pls, r_squared_local)
)

print(results, row.names = FALSE)

cat("\nBest method:", results$Method[which.min(results$Test_MSE)], "\n")
cat("\nExample completed successfully!\n")

