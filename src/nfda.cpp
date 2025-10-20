//
//  nfda.cpp
//  Modern implementation with Rcpp attributes
//
//  Originally created by Simon Müller on 18.08.11.
//  Modernized for current C++ and Rcpp standards
//
//  Ported from the R-code written by Ferraty and Vieu. More methods are available
//  on their website http://www.math.univ-toulouse.fr/staph/npfda/ 

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "nfda.h"

//' Kernel Prediction with Fixed Bandwidth
//' 
//' Performs kernel prediction using a fixed bandwidth
//' 
//' @param DistanceMatrix Numeric matrix of distances
//' @param response Numeric vector of responses
//' @param bandwidth Numeric bandwidth parameter
//' @return Numeric vector of predictions
//' @export
// [[Rcpp::export]]
arma::vec KernelPrediction(const arma::mat& DistanceMatrix, 
                           const arma::vec& response, 
                           double bandwidth) {
    
    const int n = DistanceMatrix.n_rows;
    const int m = DistanceMatrix.n_cols;
    
    // Correct bandwidth if necessary
    arma::mat KNN = arma::sort(DistanceMatrix);
    double mm = KNN.row(1).max();
    double h = (bandwidth < mm) ? mm : bandwidth;
    
    // Calculate kernel weights
    arma::mat K = DistanceMatrix / h;
    K = arma::ones<arma::mat>(n, m) - (K % K);
    K = K % (K >= 0) % (K <= 1);
    
    // Calculate predictions
    arma::vec resp = (K.t() * response) / arma::trans(arma::sum(K, 0));
    
    return resp;
}

//' Kernel Prediction with Cross-Validation
//' 
//' Performs kernel prediction with cross-validation for bandwidth selection
//' 
//' @param DistanceMatrix Numeric matrix of distances (square, n x n)
//' @param response Numeric vector of responses
//' @param bandwidth Initial bandwidth value
//' @return List containing hseq, mse, hopt, and yhat
//' @export
// [[Rcpp::export]]
Rcpp::List KernelPredictionCV(const arma::mat& DistanceMatrix, 
                              const arma::vec& response, 
                              double bandwidth) {
    
    const int n = DistanceMatrix.n_rows;
    const int hlen = 20;  // Number of bandwidth values to try
    
    arma::vec b(hlen);
    arma::vec mse(hlen);
    arma::mat Kones = arma::ones<arma::mat>(n, n);
    arma::vec tnzeros = arma::zeros<arma::vec>(n);
    
    double h = bandwidth;
    
    for (int i = 0; i < hlen; i++) {
        h *= 1.1;
        arma::mat Kern = DistanceMatrix / h;
        Kern = Kones - (Kern % Kern);
        Kern.diag() = tnzeros;
        Kern = Kern % (Kern >= 0) % (Kern <= 1);
        
        arma::rowvec tmp = arma::sum(Kern, 0);
        
        // Ensure all denominators are non-zero
        bool logic = arma::any(tmp == 0);
        while (logic) {
            h *= 1.1;
            Kern = DistanceMatrix / h;
            Kern = Kones - (Kern % Kern);
            Kern.diag() = tnzeros;
            Kern = Kern % (Kern >= 0) % (Kern <= 1);
            tmp = arma::sum(Kern, 0);
            logic = arma::any(tmp == 0);
        }
        
        arma::vec resp = Kern.t() * response / tmp.t();
        resp = arma::square(resp - response);
        mse(i) = arma::sum(resp);
        b(i) = h;
    }
    
    arma::uword index = mse.index_min();
    
    // Calculate final prediction with optimal bandwidth
    arma::mat Kern = DistanceMatrix / b(index);
    Kern = Kones - (Kern % Kern);
    Kern.diag() = tnzeros;
    Kern = Kern % (Kern >= 0) % (Kern <= 1);
    arma::vec resp = Kern.t() * response / arma::trans(arma::sum(Kern, 0));
    
    return Rcpp::List::create(
        Rcpp::Named("hseq") = b,
        Rcpp::Named("mse") = mse / n,
        Rcpp::Named("hopt") = b(index),
        Rcpp::Named("yhat") = resp
    );
}

//' Kernel Prediction with Bootstrap
//' 
//' Performs kernel prediction with bootstrap confidence intervals
//' 
//' @param DistanceMatrixPred Numeric matrix of distances (n x m)
//' @param YLearn Numeric vector of learning responses
//' @param predresponse Numeric vector of predicted responses
//' @param BootMat Numeric matrix of bootstrap samples
//' @param neighbours Integer number of neighbours
//' @return List containing pred, mu, sigma, and Mse
//' @export
// [[Rcpp::export]]
Rcpp::List KernelPredictionBoot(const arma::mat& DistanceMatrixPred, 
                                const arma::vec& YLearn, 
                                const arma::vec& predresponse, 
                                const arma::mat& BootMat, 
                                int neighbours) {
    
    const int n = DistanceMatrixPred.n_rows;
    const int m = DistanceMatrixPred.n_cols;
    const int nb = BootMat.n_cols;
    const int hlen = neighbours;
    
    arma::mat TempMat(n, hlen);
    arma::vec tnones = arma::ones<arma::vec>(n);
    arma::vec tnbones = arma::ones<arma::vec>(nb);
    
    arma::vec bootpred(m);
    arma::vec mu(m);
    arma::vec sigma(m);
    arma::mat Mse(m, hlen);
    
    // Setup bandwidth matrix
    arma::mat KNNBand = arma::sort(DistanceMatrixPred, "ascend", 0);
    KNNBand.shed_row(0);
    
    // Loop over elements
    for (int i = 0; i < m; i++) {
        arma::vec tempvec2(hlen);
        
        // Loop over bandwidths
        for (int j = 0; j < hlen; j++) {
            double h = 0.5 * (KNNBand(j, i) + KNNBand(j + 1, i));
            arma::vec Kern = DistanceMatrixPred.col(i) / h;
            Kern = tnones - Kern % Kern;
            Kern = Kern % (Kern >= 0) % (Kern <= 1);
            TempMat.col(j) = Kern;
            
            arma::vec tempvec = (BootMat.t() * Kern) / arma::sum(Kern);
            tempvec = tempvec - (predresponse(i) * tnbones);
            tempvec = arma::square(tempvec);
            Mse(i, j) = arma::sum(tempvec);
            tempvec2(j) = Mse(i, j);
        }
        
        arma::uword index = tempvec2.index_min();
        double h = arma::sum(TempMat.col(index));
        arma::vec Kern = (TempMat.col(index) % YLearn) / h;
        bootpred(i) = arma::sum(Kern);
        
        // Calculate parameters for naive bootstrapped confidence intervals
        arma::vec tempvec = (BootMat.t() * TempMat.col(index)) / h;
        mu(i) = arma::mean(tempvec);
        sigma(i) = arma::stddev(tempvec);
    }
    
    return Rcpp::List::create(
        Rcpp::Named("pred") = bootpred,
        Rcpp::Named("mu") = mu,
        Rcpp::Named("sigma") = sigma,
        Rcpp::Named("Mse") = Mse / nb
    );
}

//' Kernel Prediction with k-Nearest Neighbors
//' 
//' Performs kernel prediction using k-nearest neighbors
//' 
//' @param DistanceMatrix Numeric matrix of distances
//' @param response Numeric vector of responses
//' @param neighbours Integer or vector of neighbour counts
//' @param local Logical indicating local or global bandwidth selection
//' @return Numeric vector of predictions
//' @export
// [[Rcpp::export]]
arma::vec KernelPredictionkNN(const arma::mat& DistanceMatrix, 
                              const arma::vec& response, 
                              SEXP neighbours, 
                              bool local) {
    
    const int n = DistanceMatrix.n_rows;
    const int m = DistanceMatrix.n_cols;
    
    arma::mat K = DistanceMatrix;
    arma::rowvec h(m);
    
    if (local) {
        // Local bandwidth selection
        arma::ivec knn = Rcpp::as<arma::ivec>(neighbours);
        arma::mat KNN = arma::sort(K, "ascend", 0);
        
        for (int i = 0; i < m; i++) {
            arma::vec tmp = K.col(i);
            arma::uword index = tmp.index_min();
            int k = knn(index);
            h(i) = 0.5 * (KNN(k, i) + KNN(k + 1, i));
            K.col(i) = K.col(i) / h(i);
        }
    } else {
        // Global bandwidth selection
        int knn = Rcpp::as<int>(neighbours);
        arma::mat KNN = arma::sort(K, "ascend", 0);
        h = 0.5 * (KNN.row(knn) + KNN.row(knn + 1));
        
        for (int i = 0; i < m; i++) {
            K.col(i) = K.col(i) / h(i);
        }
    }
    
    // Calculate kernel estimate
    K = arma::ones<arma::mat>(n, m) - (K % K);
    K = K % (K >= 0) % (K <= 1);
    arma::vec resp = (K.t() * response) / arma::trans(arma::sum(K, 0));
    
    return resp;
}

//' Kernel Prediction with k-NN and Global Cross-Validation
//' 
//' Performs kernel prediction with global cross-validation for k selection
//' 
//' @param DistanceMatrix Numeric matrix of distances (square, n x n)
//' @param response Numeric vector of responses
//' @param knnlen Integer maximum number of neighbours to consider
//' @return List containing kopt, mse, and yhat
//' @export
// [[Rcpp::export]]
Rcpp::List KernelPredictionkNNgCV(const arma::mat& DistanceMatrix, 
                                  const arma::vec& response, 
                                  int knnlen) {
    
    const int n = DistanceMatrix.n_rows;
    
    arma::mat Kones = arma::ones<arma::mat>(n, n);
    arma::vec tnzeros = arma::zeros<arma::vec>(n);
    
    // Setup bandwidth matrix
    arma::mat KNN = arma::sort(DistanceMatrix, "ascend", 0);
    KNN.shed_row(0);
    
    arma::vec mse(knnlen);
    
    for (int i = 0; i < knnlen; i++) {
        arma::rowvec h = 0.5 * (KNN.row(i) + KNN.row(i + 1));
        arma::mat K(n, n);
        
        for (int j = 0; j < n; j++) {
            K.col(j) = DistanceMatrix.col(j) / h(j);
        }
        
        K = Kones - (K % K);
        K.diag() = tnzeros;
        K = K % (K >= 0) % (K <= 1);
        
        arma::vec resp = K.t() * response / arma::trans(arma::sum(K, 0));
        resp = arma::square(resp - response);
        mse(i) = arma::sum(resp);
    }
    
    arma::uword index = mse.index_min();
    
    // Calculate final prediction with optimal k
    arma::rowvec h = 0.5 * (KNN.row(index) + KNN.row(index + 1));
    arma::mat K(n, n);
    
    for (int j = 0; j < n; j++) {
        K.col(j) = DistanceMatrix.col(j) / h(j);
    }
    
    K = Kones - (K % K);
    K = K % (K >= 0) % (K <= 1);
    arma::vec resp = K.t() * response / arma::trans(arma::sum(K, 0));
    
    return Rcpp::List::create(
        Rcpp::Named("kopt") = static_cast<int>(index),
        Rcpp::Named("mse") = mse,
        Rcpp::Named("yhat") = resp
    );
}

//' Kernel Prediction with k-NN and Local Cross-Validation
//' 
//' Performs kernel prediction with local cross-validation for k selection
//' 
//' @param DistanceMatrix Numeric matrix of distances (square, n x n)
//' @param response Numeric vector of responses
//' @param knnlen Integer maximum number of neighbours to consider
//' @return List containing kopt, yhat, and mse
//' @export
// [[Rcpp::export]]
Rcpp::List KernelPredictionkNNlCV(const arma::mat& DistanceMatrix, 
                                  const arma::vec& response, 
                                  int knnlen) {
    
    const int n = DistanceMatrix.n_rows;
    
    arma::vec tknn_lenones = arma::ones<arma::vec>(knnlen);
    arma::uvec bandwidth_opt(n);
    arma::vec yhat(n);
    
    for (int i = 0; i < n; i++) {
        arma::uvec tmp_ind = arma::sort_index(DistanceMatrix.col(i));
        arma::vec tmp = arma::sort(DistanceMatrix.col(i));
        
        arma::vec z(knnlen);
        arma::vec ytmp(knnlen);
        arma::vec bandwidth(knnlen);
        
        for (int k = 0; k < knnlen; k++) {
            z(k) = tmp(k + 1);
            ytmp(k) = response(tmp_ind(k + 1));
            bandwidth(k) = 0.5 * (tmp(k + 1) + tmp(k + 2));
        }
        
        // Calculate yhat for each bandwidth
        arma::vec yhat_tmp(knnlen);
        for (int j = 0; j < knnlen; j++) {
            arma::vec ktemp = z / bandwidth(j);
            ktemp = tknn_lenones - (ktemp % ktemp);
            ktemp = ktemp % (ktemp >= 0) % (ktemp <= 1);
            
            double denom = arma::accu(ktemp);
            yhat_tmp(j) = arma::accu(ktemp % ytmp) / denom;
        }
        
        // Choose best one
        arma::vec criterium = arma::abs(response(i) * tknn_lenones - yhat_tmp);
        arma::uword index = criterium.index_min();
        yhat(i) = yhat_tmp(index);
        bandwidth_opt(i) = index;
    }
    
    // Calculate MSE
    double mse_val = arma::accu(arma::square(yhat - response));
    
    return Rcpp::List::create(
        Rcpp::Named("kopt") = bandwidth_opt,
        Rcpp::Named("yhat") = yhat,
        Rcpp::Named("mse") = mse_val
    );
}

//' Kernel Classification with k-NN and Local Cross-Validation
//' 
//' Performs kernel classification with local cross-validation for k selection
//' 
//' @param DistanceMatrix Numeric matrix of distances (square, n x n)
//' @param classes Integer vector of class labels
//' @param knnlen Integer maximum number of neighbours to consider
//' @return List containing kopt, classes.estimated, and Prob.estimated
//' @export
// [[Rcpp::export]]
Rcpp::List KernelClassificationkNNlCV(const arma::mat& DistanceMatrix, 
                                       const arma::ivec& classes, 
                                       int knnlen) {
    
    const int n = DistanceMatrix.n_rows;
    const int nbclass = arma::max(classes);
    
    arma::vec tknn_lenones = arma::ones<arma::vec>(knnlen);
    arma::uvec bandwidth_opt(n);
    arma::uvec classes_est(n);
    
    // Create binary matrix for classes
    arma::umat Binary(n, nbclass);
    for (int g = 1; g <= nbclass; g++) {
        Binary.col(g - 1) = (classes == g);
    }
    
    arma::mat Probhat(n, nbclass);
    
    for (int i = 0; i < n; i++) {
        arma::uvec tmp_ind = arma::sort_index(DistanceMatrix.col(i));
        arma::vec tmp = arma::sort(DistanceMatrix.col(i));
        
        arma::vec z(knnlen);
        arma::vec bandwidth(knnlen);
        
        for (int k = 0; k < knnlen; k++) {
            z(k) = tmp(k + 1);
            bandwidth(k) = 0.5 * (tmp(k + 1) + tmp(k + 2));
        }
        
        // Calculate class probabilities for each bandwidth
        arma::mat yhat_tmp(knnlen, nbclass);
        arma::vec criterium(knnlen);
        
        for (int j = 0; j < knnlen; j++) {
            arma::vec ktemp = z / bandwidth(j);
            ktemp = tknn_lenones - (ktemp % ktemp);
            ktemp = ktemp % (ktemp >= 0) % (ktemp <= 1);
            
            double denom = arma::accu(ktemp);
            
            for (int g = 0; g < nbclass; g++) {
                arma::vec ytmp(knnlen);
                for (int k = 0; k < knnlen; k++) {
                    ytmp(k) = Binary(tmp_ind(k + 1), g);
                }
                yhat_tmp(j, g) = arma::accu(ktemp % ytmp) / denom;
            }
            
            criterium(j) = arma::accu(arma::square(yhat_tmp.row(j) - Binary.row(i)));
        }
        
        arma::uword index = criterium.index_min();
        bandwidth_opt(i) = index;
        Probhat.row(i) = yhat_tmp.row(index);
        
        arma::rowvec tmp2 = yhat_tmp.row(index);
        arma::uword class_index = tmp2.index_max();
        classes_est(i) = class_index + 1;
    }
    
    return Rcpp::List::create(
        Rcpp::Named("kopt") = bandwidth_opt,
        Rcpp::Named("classes.estimated") = classes_est,
        Rcpp::Named("Prob.estimated") = Probhat
    );
}

//' Kernel Classification with k-NN
//' 
//' Performs kernel classification using k-nearest neighbors
//' 
//' @param DistanceMatrix Numeric matrix of distances (n x m)
//' @param classes Integer vector of class labels
//' @param neighbours Integer vector of neighbour counts
//' @return List containing classes.predicted and Prob.predicted
//' @export
// [[Rcpp::export]]
Rcpp::List KernelClassificationkNN(const arma::mat& DistanceMatrix, 
                                    const arma::ivec& classes, 
                                    const arma::ivec& neighbours) {
    
    const int n = DistanceMatrix.n_rows;
    const int m = DistanceMatrix.n_cols;
    const int nbclass = arma::max(classes);
    
    arma::mat K = DistanceMatrix;
    arma::rowvec h(m);
    
    // Create binary matrix for classes
    arma::umat Binary(n, nbclass);
    for (int g = 1; g <= nbclass; g++) {
        Binary.col(g - 1) = (classes == g);
    }
    
    arma::mat KNN = arma::sort(K, "ascend", 0);
    arma::uvec classes_pred(m);
    arma::mat Probhat(m, nbclass);
    
    for (int i = 0; i < m; i++) {
        arma::vec tmp = K.col(i);
        arma::uword index = tmp.index_min();
        int k = neighbours(index);
        h(i) = 0.5 * (KNN(k, i) + KNN(k + 1, i));
        K.col(i) = K.col(i) / h(i);
    }
    
    K = arma::ones<arma::mat>(n, m) - (K % K);
    K = K % (K >= 0) % (K <= 1);
    
    for (int g = 0; g < nbclass; g++) {
        Probhat.col(g) = K.t() * arma::conv_to<arma::vec>::from(Binary.col(g)) / arma::trans(arma::sum(K, 0));
    }
    
    for (int i = 0; i < m; i++) {
        arma::rowvec tmp2 = Probhat.row(i);
        arma::uword index = tmp2.index_max();
        classes_pred(i) = index + 1;
    }
    
    return Rcpp::List::create(
        Rcpp::Named("classes.predicted") = classes_pred,
        Rcpp::Named("Prob.predicted") = Probhat
    );
}
