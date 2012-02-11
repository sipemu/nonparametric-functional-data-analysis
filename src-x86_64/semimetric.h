//
//  semimetric.h
//  
//
//  Created by Simon Müller on 18.08.11.
//

#ifndef _semimetric_h
#define _semimetric_h

#include <RcppArmadillo.h>

// Semimetric based on PCA
RcppExport SEXP SemimetricPCAEV(SEXP Data1, SEXP q);
RcppExport SEXP SemimetricPCA(SEXP Eigenvec, SEXP Data1, SEXP Data2, SEXP twodatasets);

// Semimetric based on Derivatives
RcppExport SEXP SemimetricDerivDesign(SEXP BsplineDeriv, SEXP weightgauss, SEXP span);
RcppExport SEXP SemimetricDeriv(SEXP Hhalf, SEXP Bspline, SEXP Data1, SEXP Data2, SEXP twodatasets);

#endif
