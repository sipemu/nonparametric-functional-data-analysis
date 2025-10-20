# nfda Package Examples

This directory contains comprehensive examples demonstrating the use of the nfda package for nonparametric functional data analysis.

## Available Examples

### 1. Fat Regression Example (`fat_regression_example.R`)

Demonstrates nonparametric regression using the Fat spectrometric dataset.

**What it shows:**
- How to prepare functional data for regression
- Multiple bandwidth selection methods (kNN global, kNN local)
- Different semimetric types (Derivative, PCA)
- Model evaluation with MSE and R-squared
- Comparison of different approaches

**Run it:**
```bash
Rscript examples/fat_regression_example.R
```

**Expected output:**
- Model fitting results
- Test set performance metrics
- Method comparison table
- Best method identification

**Typical results:**
- Derivative-based: R² ≈ 0.98
- PCA-based: R² ≈ 0.38

### 2. Fat Classification Example (`fat_classification_example.R`)

Demonstrates nonparametric classification using the Fat spectrometric dataset.

**What it shows:**
- Binary classification (2 classes)
- Multi-class classification (3 classes)
- Different semimetric types (Derivative, PCA)
- Confusion matrices and accuracy metrics
- Probability distribution analysis

**Run it:**
```bash
Rscript examples/fat_classification_example.R
```

**Expected output:**
- Classification accuracy for different problems
- Confusion matrices
- Per-class performance metrics
- Prediction probability examples

**Typical results:**
- Binary classification: ~91% accuracy
- Three-class classification: ~77% accuracy

## Dataset Information

Both examples use the **Fat spectrometric dataset** from the fds package:

- **Source**: Near-infrared absorbance spectra of meat samples
- **Samples**: 215 meat samples
- **Features**: 100 wavelength measurements per sample
- **Response**: Fat content (0.9% to 49.1%)

This is a real-world dataset commonly used in chemometrics for demonstrating functional data analysis methods.

## Requirements

To run these examples, you need:
```r
install.packages(c("nfda", "fds"))
```

## How to Use

### Running from R

```r
# Source the example file
source("examples/fat_regression_example.R")

# Or run specific parts interactively
library(nfda)
library(fds)
data(Fatspectrum)
data(Fatvalues)
# ... your analysis
```

### Running from Command Line

```bash
# Make executable (optional)
chmod +x examples/fat_regression_example.R

# Run
./examples/fat_regression_example.R
# or
Rscript examples/fat_regression_example.R
```

## Understanding the Output

### Regression Example Output

```
====================================
Fat Content Prediction Example
====================================

Dataset information:
  - Number of curves: 100 
  - Number of observations: 215 
  - Fat content range: 0.9 49.1 

Example 1: Derivative-based semimetric
---------------------------------------
Fitting model with kNN global CV...
  Optimal k: 7 
  Training MSE: ...

Test set performance:
  MSE: 3.44
  MAE: 1.40
  R-squared: 0.976
```

### Classification Example Output

```
====================================
Fat Content Classification Example
====================================

Example 1: Binary Classification
---------------------------------
Test set performance:
  Accuracy: 0.9077 

Confusion Matrix:
         Actual
Predicted  1  2
        1 34  6
        2  0 25
```

## Tips for Adaptation

Want to use these examples with your own data? Here's how:

1. **Replace the data loading**:
   ```r
   # Instead of:
   data(Fatspectrum)
   X <- t(Fatspectrum$y)
   
   # Use your own:
   X <- your_functional_data_matrix  # samples × features
   Y <- your_response_vector
   ```

2. **Adjust parameters**:
   ```r
   params <- list(
     q = 2,              # Adjust derivative order or PC count
     nknot = 10,         # Adjust knot count for your data density
     range.grid = c(0, 1) # Match your data domain
   )
   ```

3. **Choose appropriate method**:
   - Use `"Deriv"` for smooth functional data
   - Use `"PCA"` for more general functional data
   - Try both and compare!

4. **Select bandwidth method**:
   - `"kNNgCV"`: Global bandwidth (same for all)
   - `"kNNlCV"`: Local bandwidth (adaptive, may be slower)
   - `"CV"`: Traditional cross-validation

## Common Issues

### Package not found
```r
# Install required packages
install.packages("fds")
```

### Data format issues
Make sure your data is in the right format:
- X: matrix with samples in rows, features in columns
- Y: numeric vector (regression) or integer vector (classification)

### Memory issues with large datasets
For very large datasets, consider:
- Using a subset for initial testing
- Reducing the number of knots (`nknot` parameter)
- Using fewer principal components (`q` parameter)

## Additional Resources

- Package documentation: `?FuNopaRe`, `?FuNopaCl`
- Original paper: Ferraty and Vieu (2006)
- fds package: `help(package="fds")`

## Contributing

Have a cool example to share? Feel free to contribute!

