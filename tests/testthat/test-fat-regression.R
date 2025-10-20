# Test nonparametric regression with Fat dataset
library(testthat)
library(nfda)
library(fds)

test_that("Regression with Fat dataset works", {
  skip_if_not_installed("fds")
  
  # Load data
  data(Fatspectrum, package = "fds")
  data(Fatvalues, package = "fds")
  
  # Prepare data: transpose to get samples x features
  X <- t(Fatspectrum$y)
  Y <- Fatvalues
  
  # Check data dimensions
  expect_equal(nrow(X), length(Y))
  expect_equal(nrow(X), 215)
  expect_equal(ncol(X), 100)
  
  # Split into training and test sets
  set.seed(123)
  n <- nrow(X)
  train_idx <- sample(1:n, size = floor(0.7 * n))
  test_idx <- setdiff(1:n, train_idx)
  
  X_train <- X[train_idx, ]
  Y_train <- Y[train_idx]
  X_test <- X[test_idx, ]
  Y_test <- Y[test_idx]
  
  # Test with derivative-based semimetric
  params_deriv <- list(
    q = 2,
    nknot = 10,
    range.grid = c(0, 1)
  )
  
  # Fit model with k-NN global CV
  model_knn <- FuNopaRe(X_train, Y_train, "Deriv", params_deriv, "kNNgCV")
  
  expect_s3_class(model_knn, "FuNopaRe")
  expect_true(!is.null(model_knn$k.opt))
  expect_true(!is.null(model_knn$Y.hat))
  expect_true(length(model_knn$Y.hat) == length(Y_train))
  expect_true(model_knn$k.opt >= 0 && model_knn$k.opt <= 20)
  
  # Make predictions
  predictions <- predict(model_knn, X_test)
  
  expect_true(!is.null(predictions$Prediction))
  expect_equal(length(predictions$Prediction), nrow(X_test))
  
  # Check prediction quality (MSE should be reasonable)
  mse <- mean((predictions$Prediction - Y_test)^2)
  expect_true(is.finite(mse))
  expect_true(mse > 0)
  
  # Test with PCA-based semimetric
  params_pca <- list(
    q = 10,
    EigenVec = NULL
  )
  
  model_pca <- FuNopaRe(X_train, Y_train, "PCA", params_pca, "kNNgCV")
  
  expect_s3_class(model_pca, "FuNopaRe")
  expect_true(!is.null(model_pca$k.opt))
  
  predictions_pca <- predict(model_pca, X_test)
  expect_equal(length(predictions_pca$Prediction), nrow(X_test))
})

test_that("Different bandwidth selection methods work with Fat dataset", {
  skip_if_not_installed("fds")
  
  data(Fatspectrum, package = "fds")
  data(Fatvalues, package = "fds")
  
  X <- t(Fatspectrum$y)[1:50, ]  # Use subset for faster testing
  Y <- Fatvalues[1:50]
  
  params <- list(q = 2, nknot = 10, range.grid = c(0, 1))
  
  # Test kNNgCV
  model1 <- FuNopaRe(X, Y, "Deriv", params, "kNNgCV")
  expect_s3_class(model1, "FuNopaRe")
  expect_equal(model1$Method, "KernelPredictionkNNgCV")
  
  # Test kNNlCV
  model2 <- FuNopaRe(X, Y, "Deriv", params, "kNNlCV")
  expect_s3_class(model2, "FuNopaRe")
  expect_equal(model2$Method, "KernelPredictionkNNlCV")
  expect_true(is.vector(model2$k.opt))
  expect_equal(length(model2$k.opt), length(Y))
})

