//
//  semimetric.cpp
//  Modern implementation with Rcpp attributes
//
//  Originally created by Simon Müller on 18.08.11.
//  Modernized for current C++ and Rcpp standards
//
//  Original written in R by Ferraty and Vieu. More semimetrics are available
//  on their website http://www.math.univ-toulouse.fr/staph/npfda/ 

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "semimetric.h"

//' Semimetric PCA Eigenvector Calculation
//' 
//' Calculates eigenvectors for PCA-based semimetric
//' 
//' @param Data1 Numeric matrix of functional data
//' @param q Integer number of principal components
//' @return Numeric matrix of eigenvectors
//' @keywords internal
// [[Rcpp::export]]
arma::mat SemimetricPCAEV(const arma::mat& Data1, int q) {
    
    const int n = Data1.n_rows;
    
    // Calculate covariance matrix and eigendecomposition
    arma::mat Cov = Data1.t() * Data1 / n;
    arma::vec eigval;
    arma::mat Eigvec;
    
    arma::eig_sym(eigval, Eigvec, Cov);
    
    // Select q largest eigenvectors (they are in ascending order)
    arma::mat Eigvecnq(Eigvec.n_rows, q);
    for (int i = 0; i < q; i++) {
        Eigvecnq.col(i) = Eigvec.col(Eigvec.n_cols - i - 1);
    }
    
    return Eigvecnq;
}

//' Semimetric Based on PCA
//' 
//' Computes semimetric distance based on principal component analysis
//' 
//' @param Eigenvec Numeric matrix of eigenvectors
//' @param Data1 Numeric matrix of functional data (first dataset)
//' @param Data2 Numeric matrix of functional data (second dataset)
//' @param twodatasets Logical indicating if two datasets are used
//' @return Numeric matrix of semimetric distances
//' @keywords internal
// [[Rcpp::export]]
arma::mat SemimetricPCA(const arma::mat& Eigenvec, 
                        const arma::mat& Data1, 
                        const arma::mat& Data2, 
                        bool twodatasets) {
    
    const int n = Data1.n_rows;
    const int q = Eigenvec.n_cols;
    
    // Calculate principal components for first dataset
    arma::mat Comp1 = Data1 * Eigenvec;
    
    if (!twodatasets) {
        // Single dataset case
        arma::mat Semimetric = arma::zeros<arma::mat>(n, n);
        
        for (int i = 0; i < q; i++) {
            arma::vec tempvec = Comp1.col(i);
            // Outer difference: outer(x, y, "-") creates matrix where M[i,j] = x[i] - y[j]
            for (arma::uword row = 0; row < n; row++) {
                for (arma::uword col = 0; col < n; col++) {
                    double diff = tempvec(row) - tempvec(col);
                    Semimetric(row, col) += diff * diff;
                }
            }
        }
        return arma::sqrt(Semimetric);
    } else {
        // Two datasets case
        const int m = Data2.n_rows;
        arma::mat Comp2 = Data2 * Eigenvec;
        arma::mat Semimetric = arma::zeros<arma::mat>(n, m);
        
        for (int i = 0; i < q; i++) {
            arma::vec tempvec1 = Comp1.col(i);
            arma::vec tempvec2 = Comp2.col(i);
            
            for (arma::uword row = 0; row < n; row++) {
                for (arma::uword col = 0; col < m; col++) {
                    double diff = tempvec1(row) - tempvec2(col);
                    Semimetric(row, col) += diff * diff;
                }
            }
        }
        return arma::sqrt(Semimetric);
    }
}

//' Semimetric Derivative Design Matrix
//' 
//' Calculates design matrix for derivative-based semimetric
//' 
//' @param BsplineDeriv Numeric matrix of B-spline derivatives
//' @param weightgauss Numeric vector of Gaussian quadrature weights
//' @param span Numeric span of the grid
//' @return Numeric matrix H-half matrix
//' @keywords internal
// [[Rcpp::export]]
arma::mat SemimetricDerivDesign(const arma::mat& BsplineDeriv, 
                                const arma::vec& weightgauss, 
                                double span) {
    
    const int n = BsplineDeriv.n_rows;
    const int m = BsplineDeriv.n_cols;
    
    // Apply Gaussian quadrature weights
    arma::vec wg = 0.5 * span * weightgauss;
    arma::mat Temp1(n, m);
    
    for (int i = 0; i < m; i++) {
        Temp1.col(i) = BsplineDeriv.col(i) % wg;
    }
    
    arma::mat H = BsplineDeriv.t() * Temp1;
    
    // Eigendecomposition
    arma::vec eigval;
    arma::mat eigvec;
    arma::eig_sym(eigval, eigvec, H);
    
    // Take square root of positive eigenvalues
    eigval = arma::sqrt(eigval % (eigval > 0));
    
    // Compute H^{1/2}
    for (arma::uword i = 0; i < H.n_rows; i++) {
        H.col(i) = eigval(i) * eigvec.col(i);
    }
    H = (eigvec * H.t()).t();
    
    return H;
}

//' Semimetric Based on Derivatives
//' 
//' Computes semimetric distance based on functional derivatives
//' 
//' @param Hhalf Numeric matrix H-half from SemimetricDerivDesign
//' @param Bspline Numeric matrix of B-spline basis
//' @param Data1 Numeric matrix of functional data (first dataset)
//' @param Data2 Numeric matrix of functional data (second dataset)
//' @param twodatasets Logical indicating if two datasets are used
//' @return Numeric matrix of semimetric distances
//' @keywords internal
// [[Rcpp::export]]
arma::mat SemimetricDeriv(const arma::mat& Hhalf, 
                          const arma::mat& Bspline, 
                          const arma::mat& Data1, 
                          const arma::mat& Data2, 
                          bool twodatasets) {
    
    const int n = Data1.n_rows;
    const int nbasis = Hhalf.n_rows;
    
    // Calculate coefficient matrices
    arma::mat Cmat = Bspline.t() * Bspline;
    arma::mat Dmat1 = Bspline.t() * Data1.t();
    
    arma::mat Coef1 = (Hhalf * arma::solve(Cmat, Dmat1)).t();
    
    if (!twodatasets) {
        // Single dataset case
        arma::mat Semimetric = arma::zeros<arma::mat>(n, n);
        
        for (int i = 0; i < nbasis; i++) {
            arma::vec tempvec = Coef1.col(i);
            
            for (arma::uword row = 0; row < n; row++) {
                for (arma::uword col = 0; col < n; col++) {
                    double diff = tempvec(row) - tempvec(col);
                    Semimetric(row, col) += diff * diff;
                }
            }
        }
        return arma::sqrt(Semimetric);
    } else {
        // Two datasets case
        const int m = Data2.n_rows;
        arma::mat Dmat2 = Bspline.t() * Data2.t();
        arma::mat Coef2 = (Hhalf * arma::solve(Cmat, Dmat2)).t();
        
        arma::mat Semimetric = arma::zeros<arma::mat>(n, m);
        
        for (int i = 0; i < nbasis; i++) {
            arma::vec tempvec1 = Coef1.col(i);
            arma::vec tempvec2 = Coef2.col(i);
            
            for (arma::uword row = 0; row < n; row++) {
                for (arma::uword col = 0; col < m; col++) {
                    double diff = tempvec1(row) - tempvec2(col);
                    Semimetric(row, col) += diff * diff;
                }
            }
        }
        return arma::sqrt(Semimetric);
    }
}
