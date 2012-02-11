//
//  semimetric.cpp
//  
//
//  Created by Simon Müller on 18.08.11.
//
//  Original written in R by Ferraty and Vieu. More semimetrics are available
//  on their website http://www.math.univ-toulouse.fr/staph/npfda/ 

#include "semimetric.h"

/*
 
        SEMIMETRIC BASED ON PCA
 
*/
RcppExport SEXP SemimetricPCAEV(SEXP Data1, SEXP q) {

    Rcpp::NumericMatrix D1r(Data1);
    int qq = Rcpp::as<int>(q);
    int n = D1r.nrow();
    arma::mat D1(D1r.begin(), D1r.nrow(), D1r.ncol(), false);
    arma::vec eigval;
    arma::mat Eigvec;
    arma::mat Eigvecnq;
    
    arma::mat Cov = trans(D1) * D1 / n;
    arma::eig_sym(eigval, Eigvec, Cov);
    Eigvecnq.zeros(Eigvec.n_rows, qq);
    for(int i = 0; i < qq; i++) {
        Eigvecnq.col(i) = Eigvec.col(Eigvec.n_cols - i - 1);
    }
    return Rcpp::wrap(Eigvecnq);
}

RcppExport SEXP SemimetricPCA(SEXP Eigenvec, SEXP Data1, SEXP Data2, SEXP twodatasets) {
    
    Rcpp::NumericMatrix D1r(Data1);
    Rcpp::NumericMatrix Evr(Eigenvec);
    bool tds = Rcpp::as<bool>(twodatasets);
    int n = D1r.nrow();
    int q = Evr.ncol();
    
    Rcpp::NumericMatrix Comp1r(n, q);
    arma::mat Comp1(Comp1r.begin(), Comp1r.nrow(), Comp1r.ncol(), false);
    arma::mat D1(D1r.begin(), D1r.nrow(), D1r.ncol(), false);
    arma::mat Ev(Evr.begin(), Evr.nrow(), Evr.ncol(), false);
    
    Comp1 = D1 * Ev;
    
    if (!tds) {
        arma::mat Semimetric;
        Semimetric.zeros(n, n);
        Rcpp::NumericMatrix Tempr;
        Rcpp::NumericVector tempvec;
        for (int i = 0; i < q; i++) {
            tempvec = Comp1r(Rcpp::_, i);
            Tempr = Rcpp::outer(tempvec, tempvec, std::minus<double>());
            Semimetric = Semimetric + arma::square(Rcpp::as<arma::mat>(Tempr));
        }
        return Rcpp::wrap(arma::sqrt(Semimetric));
    } else {
        Rcpp::NumericMatrix D2r(Data2);
        int m = D2r.nrow();
        Rcpp::NumericMatrix Comp2r(m, q);
        arma::mat Comp2(Comp2r.begin(), Comp2r.nrow(), Comp2r.ncol(), false);
        arma::mat D2(D2r.begin(), D2r.nrow(), D2r.ncol(), false);
        
        arma::mat Semimetric;
        Semimetric.zeros(n, m);
        Rcpp::NumericMatrix Tempr;
        
        Comp2 = D2 * Ev;
        
        Rcpp::NumericVector tempvec1;
        Rcpp::NumericVector tempvec2;
        for (int i = 0; i < q; i++) {
            tempvec1 = Comp1r(Rcpp::_, i);
            tempvec2 = Comp2r(Rcpp::_, i);
            Tempr = Rcpp::outer(tempvec1, tempvec2, std::minus<double>());
            Semimetric = Semimetric + arma::square(Rcpp::as<arma::mat>(Tempr));
        }
        return Rcpp::wrap(arma::sqrt(Semimetric)); 
    }
}


/*
 
 SEMIMETRIC BASED ON DERIVATIVES
 
*/
RcppExport SEXP SemimetricDerivDesign(SEXP BsplineDeriv, SEXP weightgauss, SEXP span) {
    
    Rcpp::NumericMatrix BDr(BsplineDeriv);
    Rcpp::NumericVector wgr(weightgauss);
    double sp = Rcpp::as<double>(span);
    int n = BDr.nrow();
    int m = BDr.ncol();
    arma::mat BD(BDr.begin(), n, m, false);
    arma::mat Temp1(n, m);
    arma::colvec wg(wgr.begin(), wgr.size(), true);
    
    wg = 0.5 * sp * wg;
    for (int i = 0; i < m; i++) {
        Temp1.col(i) = BD.col(i) % wg;
    }
    arma::mat H = trans(BD) * Temp1;
    arma::vec eigval;
    arma::mat eigvec;
    
    arma::eig_sym(eigval, eigvec, H);
    eigval = arma::sqrt(eigval % (eigval > 0));
    
    for (int i = 0; i < H.n_rows; i++) {
        H.col(i) = eigval(i) * eigvec.col(i);
    }
    H = arma::trans(eigvec * arma::trans(H));
    
    return Rcpp::wrap(H);
}

RcppExport SEXP SemimetricDeriv(SEXP Hhalf, SEXP Bspline, SEXP Data1, SEXP Data2, SEXP twodatasets) {
    
    bool tds = Rcpp::as<bool>(twodatasets);
    if (!tds) {
        Rcpp::NumericMatrix Hr(Hhalf);
        Rcpp::NumericMatrix Bsr(Bspline);
        Rcpp::NumericMatrix Dr(Data1);
        int n = Dr.nrow();
        int nbasis = Hr.nrow();
        
        arma::mat H(Hr.begin(), Hr.nrow(), Hr.ncol(), false);
        arma::mat Bs(Bsr.begin(), Bsr.nrow(), Bsr.ncol(), false);
        arma::mat D(Dr.begin(), Dr.nrow(), Dr.ncol(), false);
        
        arma::mat Cmat = arma::trans(Bs) * Bs;
        arma::mat Dmat1 = arma::trans(Bs) * arma::trans(D);
        
        Rcpp::NumericMatrix C1r(n, nbasis);
        arma::mat Coef1(C1r.begin(), C1r.nrow(), C1r.ncol(), false);
        Coef1 = arma::trans(H * arma::solve(Cmat, Dmat1));
        
        Rcpp::NumericMatrix tempr;
        arma::mat Semimetric;
        Semimetric.zeros(n, n);
        
        Rcpp::NumericVector tempvec = C1r(Rcpp::_, 0);
        for (int i = 0; i < nbasis; i++) {
            tempvec = C1r(Rcpp::_, i);
            tempr = Rcpp::outer(tempvec, tempvec, std::minus<double>());
            Semimetric = Semimetric + arma::square(Rcpp::as<arma::mat>(tempr));
        }
        return Rcpp::wrap(arma::sqrt(Semimetric)); 
    } else {
        Rcpp::NumericMatrix Hr(Hhalf);
        Rcpp::NumericMatrix Bsr(Bspline);
        Rcpp::NumericMatrix D1r(Data1);
        Rcpp::NumericMatrix D2r(Data2);
        int n = D1r.nrow();
        int m = D2r.nrow();
        int nbasis = Hr.nrow();
        
        arma::mat H(Hr.begin(), Hr.nrow(), Hr.ncol(), false);
        arma::mat Bs(Bsr.begin(), Bsr.nrow(), Bsr.ncol(), false);
        arma::mat D1(D1r.begin(), D1r.nrow(), D1r.ncol(), false);
        arma::mat D2(D2r.begin(), D2r.nrow(), D2r.ncol(), false);
        
        arma::mat Cmat = arma::trans(Bs) * Bs;
        arma::mat Dmat1 = arma::trans(Bs) * arma::trans(D1);
        arma::mat Dmat2 = arma::trans(Bs) * arma::trans(D2);
        
        Rcpp::NumericMatrix C1r(n, nbasis);
        Rcpp::NumericMatrix C2r(m, nbasis);
        arma::mat Coef1(C1r.begin(), C1r.nrow(), C1r.ncol(), false);
        arma::mat Coef2(C2r.begin(), C2r.nrow(), C2r.ncol(), false);
        Coef1 = arma::trans(H * arma::solve(Cmat, Dmat1));
        Coef2 = arma::trans(H * arma::solve(Cmat, Dmat2));
        
        
        Rcpp::NumericMatrix tempr(n, m);
        arma::mat temp(tempr.begin(), n, m, false);
        arma::mat Semimetric;
        Semimetric.zeros(n, m);
        
        Rcpp::NumericVector tempvec1 = C1r(Rcpp::_, 0);
        Rcpp::NumericVector tempvec2 = C2r(Rcpp::_, 0);
        for (int i = 0; i < nbasis; i++) {
            tempvec1 = C1r( Rcpp::_, i);
            tempvec2 = C2r( Rcpp::_, i);
            tempr = Rcpp::outer(tempvec1, tempvec2, std::minus<double>());
            Semimetric = Semimetric + arma::square(Rcpp::as<arma::mat>(tempr));
        }
        return Rcpp::wrap(arma::sqrt(Semimetric)); 
    }
}
