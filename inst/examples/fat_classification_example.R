#!/usr/bin/env Rscript
#' ---
#' title: "Fat Content Classification using Nonparametric Methods"
#' author: "nfda package example"
#' date: "`r Sys.Date()`"
#' ---
#' 
#' This example demonstrates how to use the nfda package for nonparametric
#' functional classification. We convert the fat content prediction problem
#' into classification by categorizing samples into low, medium, and high
#' fat content groups.

# Load required packages
library(nfda)
library(fds)

cat("====================================\n")
cat("Fat Content Classification Example\n")
cat("====================================\n\n")

# Load the Fat dataset
data(Fatspectrum, package = "fds")
data(Fatvalues, package = "fds")

# Prepare data
X <- t(Fatspectrum$y)
Y <- Fatvalues

cat("Dataset information:\n")
cat("  - Number of observations:", nrow(X), "\n")
cat("  - Number of wavelengths:", ncol(X), "\n")
cat("  - Fat content range:", range(Y), "\n\n")

# ============================================================================
# Example 1: Binary Classification (Low vs High Fat)
# ============================================================================
cat("Example 1: Binary Classification\n")
cat("---------------------------------\n")

# Create binary classes based on median
classes_binary <- as.integer(Y >= median(Y)) + 1

cat("Class distribution:\n")
cat("  Class 1 (Low fat):", sum(classes_binary == 1), "\n")
cat("  Class 2 (High fat):", sum(classes_binary == 2), "\n\n")

# Split into training and test sets
set.seed(2025)
n <- nrow(X)
train_idx <- sample(1:n, size = floor(0.7 * n))
test_idx <- setdiff(1:n, train_idx)

X_train <- X[train_idx, ]
classes_train <- classes_binary[train_idx]
X_test <- X[test_idx, ]
classes_test <- classes_binary[test_idx]

# Fit classification model
params_deriv <- list(
  q = 2,
  nknot = 10,
  range.grid = c(0, 1)
)

cat("Fitting binary classification model...\n")
model_binary <- FuNopaCl(X_train, classes_train, 
                         semimetric = "Deriv",
                         semimetric.params = params_deriv)

cat("  Training error rate:", round(model_binary$mse.learn, 4), "\n")
cat("  Optimal k range: [", min(model_binary$k.opt), ",", 
    max(model_binary$k.opt), "]\n\n")

# Make predictions
pred_binary <- predict(model_binary, X_test)

# Calculate test accuracy
accuracy <- mean(pred_binary$classes.pred == classes_test)
confusion_matrix <- table(Predicted = pred_binary$classes.pred, 
                         Actual = classes_test)

cat("Test set performance:\n")
cat("  Accuracy:", round(accuracy, 4), "\n")
cat("  Error rate:", round(1 - accuracy, 4), "\n\n")

cat("Confusion Matrix:\n")
print(confusion_matrix)
cat("\n")

# Calculate sensitivity and specificity
sensitivity <- confusion_matrix[1, 1] / sum(confusion_matrix[, 1])
specificity <- confusion_matrix[2, 2] / sum(confusion_matrix[, 2])

cat("  Sensitivity (Class 1):", round(sensitivity, 4), "\n")
cat("  Specificity (Class 2):", round(specificity, 4), "\n\n")

# ============================================================================
# Example 2: Multi-class Classification (Low, Medium, High)
# ============================================================================
cat("Example 2: Three-Class Classification\n")
cat("--------------------------------------\n")

# Create three classes based on tertiles
quantiles <- quantile(Y, probs = c(1/3, 2/3))
classes_multi <- cut(Y, 
                     breaks = c(-Inf, quantiles[1], quantiles[2], Inf),
                     labels = c(1, 2, 3))
classes_multi <- as.integer(as.character(classes_multi))

cat("Class distribution:\n")
cat("  Class 1 (Low fat):", sum(classes_multi == 1), "\n")
cat("  Class 2 (Medium fat):", sum(classes_multi == 2), "\n")
cat("  Class 3 (High fat):", sum(classes_multi == 3), "\n\n")

# Use same train/test split
classes_train_multi <- classes_multi[train_idx]
classes_test_multi <- classes_multi[test_idx]

cat("Fitting three-class classification model...\n")
model_multi <- FuNopaCl(X_train, classes_train_multi, 
                        semimetric = "Deriv",
                        semimetric.params = params_deriv)

cat("  Training error rate:", round(model_multi$mse.learn, 4), "\n\n")

# Make predictions
pred_multi <- predict(model_multi, X_test)

# Calculate test accuracy
accuracy_multi <- mean(pred_multi$classes.pred == classes_test_multi)
confusion_matrix_multi <- table(Predicted = pred_multi$classes.pred, 
                                Actual = classes_test_multi)

cat("Test set performance:\n")
cat("  Accuracy:", round(accuracy_multi, 4), "\n")
cat("  Error rate:", round(1 - accuracy_multi, 4), "\n\n")

cat("Confusion Matrix:\n")
print(confusion_matrix_multi)
cat("\n")

# Calculate per-class accuracy
for (i in 1:3) {
  if (sum(confusion_matrix_multi[, i]) > 0) {
    class_acc <- confusion_matrix_multi[i, i] / sum(confusion_matrix_multi[, i])
    cat("  Class", i, "accuracy:", round(class_acc, 4), "\n")
  }
}
cat("\n")

# ============================================================================
# Example 3: PCA-based Semimetric for Classification
# ============================================================================
cat("Example 3: PCA-based Semimetric\n")
cat("--------------------------------\n")

params_pca <- list(
  q = 20,
  EigenVec = NULL
)

cat("Fitting model with PCA semimetric...\n")
model_pca <- FuNopaCl(X_train, classes_train, 
                      semimetric = "PCA",
                      semimetric.params = params_pca)

cat("  Training error rate:", round(model_pca$mse.learn, 4), "\n\n")

# Make predictions
pred_pca <- predict(model_pca, X_test)

accuracy_pca <- mean(pred_pca$classes.pred == classes_test)

cat("Test set performance:\n")
cat("  Accuracy:", round(accuracy_pca, 4), "\n\n")

# ============================================================================
# Analyzing Prediction Probabilities
# ============================================================================
cat("Example 4: Prediction Probabilities Analysis\n")
cat("---------------------------------------------\n")

# Look at probability distribution for a few test samples
cat("Probability distributions for selected test samples:\n\n")

for (i in 1:min(5, nrow(X_test))) {
  cat("Sample", i, "(True class:", classes_test[i], "):\n")
  cat("  Predicted class:", pred_binary$classes.pred[i], "\n")
  cat("  Class 1 prob:", round(pred_binary$Prob.pred[i, 1], 4), "\n")
  cat("  Class 2 prob:", round(pred_binary$Prob.pred[i, 2], 4), "\n\n")
}

# ============================================================================
# Summary
# ============================================================================
cat("====================================\n")
cat("Summary of Classification Results\n")
cat("====================================\n\n")

results <- data.frame(
  Problem = c("Binary (Deriv)", "Multi-class (Deriv)", "Binary (PCA)"),
  Classes = c(2, 3, 2),
  Accuracy = c(accuracy, accuracy_multi, accuracy_pca),
  Error_Rate = c(1 - accuracy, 1 - accuracy_multi, 1 - accuracy_pca)
)

print(results, row.names = FALSE)

cat("\nBest binary classifier:", 
    ifelse(accuracy > accuracy_pca, "Derivative-based", "PCA-based"), "\n")
cat("\nExample completed successfully!\n")

