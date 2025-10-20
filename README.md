# nfda: Nonparametric Functional Data Analysis

[![License: GPL-2](https://img.shields.io/badge/License-GPL%202-blue.svg)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)

## Overview

The `nfda` package provides nonparametric methods for functional data analysis in R. It implements kernel-based regression and classification techniques for functional data, with various semimetric options and bandwidth selection methods.

This package is based on the seminal work of Ferraty and Vieu on nonparametric functional data analysis. The C++ implementation using Rcpp and RcppArmadillo provides efficient computation for large datasets.

## Features

### Nonparametric Regression
- Kernel prediction with fixed bandwidth
- Cross-validation for bandwidth selection
- k-Nearest Neighbors (k-NN) with global and local cross-validation
- Bootstrap confidence intervals

### Nonparametric Classification
- k-NN classification with local cross-validation
- Probability estimates for class membership

### Semimetric Options
- **PCA-based**: Distances based on principal component scores
- **Derivative-based**: Distances based on functional derivatives using B-splines
- **PLS-based**: Distances based on partial least squares components

## Installation

### From Source

```r
# Install dependencies
install.packages(c("Rcpp", "RcppArmadillo", "splines", "pls"))

# Install from source
install.packages("nfda", type = "source")
```

### Development Version

```r
# Install devtools if needed
install.packages("devtools")

# Install from GitHub (when available)
devtools::install_github("simonm/Nonparametric-Functional-Data-Analysis")
```

## Requirements

- R >= 3.5.0
- C++11 compiler
- Rcpp >= 1.0.0
- RcppArmadillo >= 0.9.0

## Quick Start

### Nonparametric Regression

```r
library(nfda)

# Generate example functional data
n <- 100
p <- 50
X <- matrix(rnorm(n * p), n, p)
Y <- rnorm(n)

# Set up semimetric parameters for derivative-based semimetric
params <- list(
  q = 2,           # Second derivative
  nknot = 10,      # Number of B-spline knots
  range.grid = c(0, 1)
)

# Fit nonparametric regression model
model <- FuNopaRe(X, Y, 
                  semimetric = "Deriv", 
                  semimetric.params = params,
                  bandwidth = "kNNgCV")

# Predict on new data
X_new <- matrix(rnorm(10 * p), 10, p)
predictions <- predict(model, X_new)
```

### Nonparametric Classification

```r
library(nfda)

# Generate example functional data
n <- 100
p <- 50
X <- matrix(rnorm(n * p), n, p)
classes <- rep(1:2, each = 50)

# Set up semimetric parameters
params <- list(
  q = 5,           # Number of principal components
  EigenVec = NULL  # Will be computed automatically
)

# Fit nonparametric classification model
model <- FuNopaCl(X, classes,
                  semimetric = "PCA",
                  semimetric.params = params)

# Predict on new data
X_new <- matrix(rnorm(10 * p), 10, p)
predictions <- predict(model, X_new)

# View predicted classes
print(predictions$classes.pred)

# View class probabilities
print(predictions$Prob.pred)
```

## Bandwidth Selection Methods

The package offers several bandwidth selection methods:

1. **CV**: Cross-validation for continuous bandwidth
2. **kNNgCV**: k-Nearest Neighbors with global cross-validation
3. **kNNlCV**: k-Nearest Neighbors with local cross-validation (adaptive bandwidth)

## Semimetric Details

### PCA-based Semimetric

Distance based on principal component scores:
```r
d(f, g) = sqrt(sum_{j=1}^q (s_j(f) - s_j(g))^2)
```
where `s_j` are the principal component scores.

### Derivative-based Semimetric

Distance based on L2-norm of functional derivatives:
```r
d(f, g) = sqrt(integral (f^(q)(t) - g^(q)(t))^2 dt)
```
where `f^(q)` denotes the q-th derivative.

### PLS-based Semimetric

Distance based on partial least squares components for supervised learning.

## Performance

The package uses C++11 with Rcpp and RcppArmadillo for efficient computation:
- Vectorized operations using Armadillo
- Const-correctness for better optimization
- Modern C++ features for improved performance

## References

Ferraty, F. and Vieu, P. (2006). *Nonparametric Functional Data Analysis: Theory and Practice*. Springer Series in Statistics.

Additional methods available at: http://www.math.univ-toulouse.fr/staph/npfda/

## Original Author

Simon Müller

Original code created in 2011, modernized in 2025.

## License

MIT License - see LICENSE file for details

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## Issues

Please report any issues at: https://github.com/simonm/Nonparametric-Functional-Data-Analysis/issues

