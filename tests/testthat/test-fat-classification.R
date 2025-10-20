# Test nonparametric classification with Fat dataset
library(testthat)
library(nfda)
library(fds)

test_that("Classification with Fat dataset works", {
  skip_if_not_installed("fds")
  
  # Load data
  data(Fatspectrum, package = "fds")
  data(Fatvalues, package = "fds")
  
  # Prepare data
  X <- t(Fatspectrum$y)
  Y <- Fatvalues
  
  # Create binary classification: low fat (< median) vs high fat (>= median)
  classes <- as.integer(Y >= median(Y)) + 1
  
  # Split into training and test sets
  set.seed(456)
  n <- nrow(X)
  train_idx <- sample(1:n, size = floor(0.7 * n))
  test_idx <- setdiff(1:n, train_idx)
  
  X_train <- X[train_idx, ]
  classes_train <- classes[train_idx]
  X_test <- X[test_idx, ]
  classes_test <- classes[test_idx]
  
  # Test with derivative-based semimetric
  params_deriv <- list(
    q = 2,
    nknot = 10,
    range.grid = c(0, 1)
  )
  
  # Fit classification model
  model <- FuNopaCl(X_train, classes_train, "Deriv", params_deriv)
  
  expect_s3_class(model, "FuNopaCl")
  expect_true(!is.null(model$k.opt))
  expect_true(!is.null(model$classes.estimated))
  expect_true(!is.null(model$Prob.estimated))
  expect_equal(length(model$classes.estimated), length(classes_train))
  expect_equal(nrow(model$Prob.estimated), length(classes_train))
  expect_equal(ncol(model$Prob.estimated), 2)  # Two classes
  
  # Check that error rate is reasonable
  expect_true(model$mse.learn >= 0 && model$mse.learn <= 1)
  
  # Make predictions
  predictions <- predict(model, X_test)
  
  expect_true(!is.null(predictions$classes.pred))
  expect_true(!is.null(predictions$Prob.pred))
  expect_equal(length(predictions$classes.pred), nrow(X_test))
  expect_equal(nrow(predictions$Prob.pred), nrow(X_test))
  
  # Check prediction accuracy
  accuracy <- mean(predictions$classes.pred == classes_test)
  expect_true(accuracy > 0.5)  # Should be better than random
  
  # Check that probability matrix has correct dimensions
  expect_true(is.matrix(predictions$Prob.pred))
  expect_equal(dim(predictions$Prob.pred), c(nrow(X_test), 2))
})

test_that("Multi-class classification works with Fat dataset", {
  skip_if_not_installed("fds")
  
  data(Fatspectrum, package = "fds")
  data(Fatvalues, package = "fds")
  
  X <- t(Fatspectrum$y)[1:100, ]  # Use subset for faster testing
  Y <- Fatvalues[1:100]
  
  # Create 3-class problem based on tertiles
  quantiles <- quantile(Y, probs = c(1/3, 2/3))
  classes <- cut(Y, 
                 breaks = c(-Inf, quantiles[1], quantiles[2], Inf),
                 labels = FALSE)
  
  params <- list(q = 2, nknot = 10, range.grid = c(0, 1))
  
  # Fit model
  model <- FuNopaCl(X, classes, "Deriv", params)
  
  expect_s3_class(model, "FuNopaCl")
  expect_equal(ncol(model$Prob.estimated), 3)  # Three classes
  expect_true(all(model$classes.estimated %in% 1:3))
  
  # Test predictions
  X_test <- X[1:10, ]
  predictions <- predict(model, X_test)
  
  expect_equal(ncol(predictions$Prob.pred), 3)
  expect_true(all(predictions$classes.pred %in% 1:3))
})

test_that("PCA semimetric works for classification", {
  skip_if_not_installed("fds")
  
  data(Fatspectrum, package = "fds")
  data(Fatvalues, package = "fds")
  
  X <- t(Fatspectrum$y)[1:80, ]
  Y <- Fatvalues[1:80]
  
  classes <- as.integer(Y >= median(Y)) + 1
  
  params_pca <- list(
    q = 15,
    EigenVec = NULL
  )
  
  model <- FuNopaCl(X, classes, "PCA", params_pca)
  
  expect_s3_class(model, "FuNopaCl")
  expect_equal(ncol(model$Prob.estimated), 2)
  
  # Predictions
  predictions <- predict(model, X[1:10, ])
  expect_equal(length(predictions$classes.pred), 10)
})

