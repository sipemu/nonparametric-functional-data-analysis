//
//  nfda.cpp
//  
//
//  Created by Simon Müller on 18.08.11.
//
//  Portet from the R-code written by Ferraty and Vieu. More methods are available
//  on their website http://www.math.univ-toulouse.fr/staph/npfda/ 

#include "nfda.h"

RcppExport SEXP KernelPrediction(SEXP DistanceMatrix, SEXP response, SEXP bandwidth) {
    
    Rcpp::NumericMatrix Kr(DistanceMatrix);
    Rcpp::NumericVector yr(response);
    double h = Rcpp::as<double>(bandwidth);
    int n = Kr.nrow();
    int m = Kr.ncol();
    
    arma::mat Kc(Kr.begin(), n, m, false); 
    arma::mat K;
    arma::colvec y(yr.begin(), yr.size(), false);
    arma::colvec resp;
    
    // correct bandwidth if necessary
    arma::mat KNN = sort(Kc);
    arma::rowvec hvec = KNN.row(1);
    double mm = hvec.max(); 
    if (h < mm) {
        h = mm;
    }
    K = Kc / h;
    K = arma::ones<arma::mat>(n, m) - (K % K);
    K = K % (K >= 0);
    K = K % (K <= 1);
    resp = (arma::trans(K) * y) / arma::trans(arma::sum(K));

    return Rcpp::wrap(resp);
}

/*
 Portet from the R-code written by Ferraty and Vieu. More methods are available
 on their website http://www.math.univ-toulouse.fr/staph/npfda/ 
 
 Comment:   I observed that I get in this global cv case always the smallest mse 
            be the smallest corrected bandwidth. Probably this is an effect of 
            heterogenity of the data sets I examined. Alternatively, one can 
            set m(X_i) = 0, if the denominator of the weight corresponding to 
            X_i is equal to zero.
 */

RcppExport SEXP KernelPredictionCV(SEXP DistanceMatrix, SEXP response, SEXP bandwidth) {

    Rcpp::NumericMatrix Kr(DistanceMatrix);
    Rcpp::NumericVector yr(response);
    Rcpp::NumericVector br(bandwidth);
    int hlen = br.size();
    int n = Kr.nrow();
    
    arma::mat K(Kr.begin(), n, n, false);
    arma::colvec y(yr.begin(), yr.size(), false);
    arma::colvec b(br.begin(), br.size(), false);
    arma::mat Kern(n, n);
    
    arma::colvec tnzeros = arma::zeros<arma::colvec>(n);
    arma::mat Kones = arma::ones<arma::mat>(n, n);
    arma::colvec resp(n);
    arma::colvec mse(hlen);
    arma::u32 index;
    arma::rowvec tmp(n);
    
    int logic;
    double h;
    h = b(0);
    for (int i = 0; i < hlen; i++) {
        h = 1.1 * h;
        Kern = K / h;
        Kern = Kones - (Kern % Kern);
        Kern.diag() = tnzeros;
        Kern = Kern % (Kern >= 0);
        Kern = Kern % (Kern <= 1);
        tmp = arma::sum(Kern);
        logic = 0;
        for (int j = 0; j < n; j++) {
            if (tmp(j) == 0) {
                logic = 1;
            } 
        }
        while (logic == 1) {
            h = 1.1 * h;
            Kern = K / h;
            Kern = Kones - (Kern % Kern);
            Kern.diag() = tnzeros;
            Kern = Kern % (Kern >= 0);
            Kern = Kern % (Kern <= 1);
            tmp = arma::sum(Kern);
            logic = 0;
            for (int j = 0; j < n; j++) {
                if (tmp(j) == 0) {
                    logic = 1;
                } 
            }
        }
        resp = arma::trans(Kern) * y / arma::trans(tmp);
        
        resp = (resp - y) % (resp - y);
        mse(i) = arma::sum(resp);
        b(i) = h;
    }
    
    mse.min(index);
    
    Kern = K / b(index);
    Kern = Kones - (Kern % Kern);
    Kern.diag() = tnzeros;
    Kern = Kern % (Kern >= 0);
    Kern = Kern % (Kern <= 1);
    resp = arma::trans(Kern) * y / arma::trans(arma::sum(Kern));
    
    return Rcpp::List::create(Rcpp::Named("hseq") = Rcpp::wrap(b), Rcpp::Named("mse") = Rcpp::wrap(mse / n), Rcpp::Named("hopt") = Rcpp::wrap(b(index)), Rcpp::Named("yhat") = Rcpp::wrap(resp));
}


RcppExport SEXP KernelPredictionBoot(SEXP DistanceMatrixPred, SEXP YLearn, SEXP predresponse, SEXP BootMat, SEXP neighbours) {
    
    Rcpp::NumericMatrix Kr(DistanceMatrixPred);
    Rcpp::NumericMatrix Wr(BootMat);
    Rcpp::NumericVector yr(YLearn);   
    Rcpp::NumericVector pyr(predresponse);
    
    int hlen = Rcpp::as<int>(neighbours);
    int n = Kr.nrow();
    int m = Kr.ncol();
    int nb = Wr.ncol();
    
    arma::mat K(Kr.begin(), n, m, false);
    arma::mat W(Wr.begin(), Wr.nrow(), Wr.ncol(), false);
    arma::colvec Kern(n);
    arma::colvec y(yr.begin(), yr.size(), false);
    arma::colvec py(pyr.begin(), pyr.size(), false);
      
    arma::mat TempMat(n, hlen);
    arma::vec tempvec(nb);
    arma::vec tempvec2(hlen);
    arma::colvec tnones = arma::ones<arma::colvec>(n);
    arma::colvec tnbones = arma::ones<arma::colvec>(nb);

    arma::colvec bootpred(m);
    arma::colvec mu(m);
    arma::colvec sigma(m);
    arma::mat Mse(m, hlen);
    double h;
    arma::u32 index;
    
    // setup bandwidth matrix
    arma::mat KNNBand = sort(K);
    KNNBand.shed_row(0);
    
    // loop over elements
    for (int i = 0; i < m; i++) {
        // loop over bandwidths
        for (int j = 0; j < hlen; j++) {
            h = 0.5 * (KNNBand.col(i)(j) + KNNBand.col(i)(j + 1));
            Kern = K.col(i) / h;
            Kern = tnones - Kern % Kern;
            Kern = Kern % (Kern >= 0);
            Kern = Kern % (Kern <= 1);
            TempMat.col(j) = Kern; 
            tempvec = (arma::trans(W) * Kern) / arma::sum(Kern);
            tempvec = tempvec - (py(i) * tnbones);
            tempvec = tempvec % tempvec;
            Mse(i, j) = arma::sum(tempvec);
            tempvec2(j) = Mse(i, j);
        }
        tempvec2.min(index);
        h = arma::sum(TempMat.col(index));
        Kern = (TempMat.col(index) % y) / h;
        bootpred(i) = arma::sum(Kern);
        
        //calculate parameters for naive bootstrapped confidence intervals 
        //TODO: improve that!!!
        tempvec = (arma::trans(W) * TempMat.col(index)) / h;
        mu(i) = arma::mean(tempvec);
        sigma(i) = arma::stddev(tempvec);
    }
    return Rcpp::List::create(Rcpp::Named("pred") = Rcpp::wrap(bootpred), Rcpp::Named("mu") = Rcpp::wrap(bootpred), Rcpp::Named("sigma") = Rcpp::wrap(bootpred), Rcpp::Named("Mse") = Rcpp::wrap(Mse / nb)); 
}


RcppExport SEXP KernelPredictionkNN(SEXP DistanceMatrix, SEXP response, SEXP neighbours, SEXP local) {
    
    Rcpp::NumericMatrix Kr(DistanceMatrix);
    Rcpp::NumericVector yr(response);
    int n = Kr.nrow();
    int m = Kr.ncol();
    
    // if local == true we have a vector of bandwidths else just an integer value
    bool loc = Rcpp::as<bool>(local);
     
    arma::mat Kc(Kr.begin(), n, m, false); 
    arma::mat K = Kc;
    arma::colvec y(yr.begin(), yr.size(), false);
    arma::colvec resp;
    
    // setup bandwidth matrix: global and local
    arma::rowvec h(m);
    if (loc) {
        Rcpp::IntegerVector knnr(neighbours);
        arma::ivec knn(knnr.begin(), knnr.size(), false);
        arma::colvec tmp(n);
        arma::u32 index;
        arma::mat KNN = sort(K);
        int k;
        for (int i = 0; i < m; i++) {
            tmp = K.col(i);
            tmp.min(index);
            k = knn(index);
            tmp = KNN.col(i);
            h(i) = 0.5 * (tmp(k) + tmp(k + 1));
            K.col(i) = K.col(i) / h(i);
        } 
    } else {
        int knn = Rcpp::as<int>(neighbours);
        arma::mat KNN = sort(K);
        h = 0.5 * (KNN.row(knn) + KNN.row(knn + 1));
        for (int i = 0; i < m; i++) {
            K.col(i) = K.col(i) / h(i);
        }
    }    
    // calculate kernel estimate
    K = arma::ones<arma::mat>(n, m) - (K % K);
    K = K % (K >= 0);
    K = K % (K <= 1);
    resp = (arma::trans(K) * y) / arma::trans(arma::sum(K));
    
    return Rcpp::wrap(resp);
}




RcppExport SEXP KernelPredictionkNNgCV(SEXP DistanceMatrix, SEXP response, SEXP knnlen) {
    
    Rcpp::NumericMatrix Kr(DistanceMatrix);
    Rcpp::NumericVector yr(response);
    int knn_len = Rcpp::as<int>(knnlen);
    int n = Kr.nrow();
    
    arma::mat Kc(Kr.begin(), n, n, false); 
    arma::mat K(n, n);
    arma::colvec y(yr.begin(), yr.size(), false);
    arma::colvec resp(n);
    
    arma::colvec tnzeros = arma::zeros<arma::colvec>(n);
    arma::mat Kones = arma::ones<arma::mat>(n, n);
    
    // setup bandwidth matrix
    arma::mat KNN = sort(Kc);
    KNN.shed_row(0);
    arma::rowvec h;
    
    arma::colvec mse(knn_len);
    arma::u32 index;
    
    for (int i = 0; i < knn_len; i++) {
        h = 0.5 * (KNN.row(i) + KNN.row(i + 1));
        for (int j = 0; j < n; j++) {
            K.col(j) = Kc.col(j) / h(j);
        }
        K = Kones - (K % K);
        K.diag() = tnzeros;
        K = K % (K >= 0);
        K = K % (K <= 1);
        resp = arma::trans(K) * y / arma::trans(arma::sum(K));        
        resp = (resp - y) % (resp - y);
        mse(i) = arma::sum(resp);
    }
    mse.min(index);
    
    h = 0.5 * (KNN.row(index) + KNN.row(index + 1));
    for (int j = 0; j < n; j++) {
        K.col(j) = Kc.col(j) / h(j);
    }
    K = Kones - (K % K);
    K = K % (K >= 0);
    K = K % (K <= 1);
    resp = arma::trans(K) * y / arma::trans(arma::sum(K));
    
    return Rcpp::List::create(Rcpp::Named("kopt") = Rcpp::wrap(index), Rcpp::Named("mse") = Rcpp::wrap(mse), Rcpp::Named("yhat") = Rcpp::wrap(resp));
}


RcppExport SEXP KernelPredictionkNNlCV(SEXP DistanceMatrix, SEXP response, SEXP knnlen) {
    
    Rcpp::NumericMatrix Kr(DistanceMatrix);
    Rcpp::NumericVector yr(response);
    int knn_len = Rcpp::as<int>(knnlen);
    int n = Kr.nrow();
    
    arma::mat Kc(Kr.begin(), n, n, false); 
    arma::colvec y(yr.begin(), yr.size(), false);
    
    arma::colvec tknn_lenones = arma::ones<arma::colvec>(knn_len);
    arma::colvec tnzeros = arma::zeros<arma::colvec>(n);
    
    arma::uvec bandwidth_opt(n);
    arma::colvec yhat(n);
    
    arma::colvec tmp(n);
    
    arma::colvec zz(knn_len + 1);
    arma::colvec criterium(knn_len);
    arma::colvec bandwidth(knn_len);
    arma::colvec z(knn_len);
    arma::colvec ytmp(knn_len);
    arma::colvec ktemp(knn_len);
    arma::colvec yhat_tmp(knn_len);
    
    arma::uvec tmp_ind(n);
    arma::uvec norm_order(knn_len);
    arma::u32 index;    
    
    double denom;
    for (int i = 0; i < n; i++) {
        tmp_ind = sort_index(Kc.col(i));
        tmp = sort(Kc.col(i));
        for (int k = 0; k <= knn_len; k++) {
            zz(k) = tmp(k + 1);
        }
        // get the correct response values
        for (int k = 0; k < knn_len; k++) {
            norm_order(k) = tmp_ind(k + 1);
            z(k) = zz(k);
            ytmp(k) = y(norm_order(k));
            // get bandwidths
            bandwidth(k) = 0.5 * (zz(k) + zz(k + 1));
        }  
        // calculate yhat for each bandwidth
        for (int j = 0; j < knn_len; j++) {
            ktemp = z / bandwidth(j);
            ktemp = tknn_lenones - (ktemp % ktemp);
            ktemp = ktemp % (ktemp >= 0);
            ktemp = ktemp % (ktemp <= 1);
            yhat_tmp(j) = arma::accu(ktemp % ytmp);
            denom = arma::accu(ktemp);
            yhat_tmp(j) = yhat_tmp(j) / denom;
        }
        // choose best one
        criterium = arma::abs(y(i) * tknn_lenones - yhat_tmp);
        criterium.min(index);
        yhat(i) = yhat_tmp(index);
        bandwidth_opt(i) = index;
    }
    // calc mse
    double mse;
    yhat = (yhat - y) % (yhat - y);
    mse = arma::sum(yhat);
    
    return Rcpp::List::create(Rcpp::Named("kopt") = Rcpp::wrap(bandwidth_opt), Rcpp::Named("yhat") = Rcpp::wrap(yhat), Rcpp::Named("mse") = Rcpp::wrap(mse));
}


/*
 Portet from the R-code written by Ferraty and Vieu. More methods are available
 on their website http://www.math.univ-toulouse.fr/staph/npfda/ 
 
 
 kNN-classification with local cross-validation
 
 
 INPUT:
 DistanceMatrix     :    n x n - distance matrix
 classes            :    Vector of classes with length n
 knnlen             :    Number of neighbours for bandwidth selection
 
 OUTPUT:
 kopt               :   Vector of length n with the optimal k
 classes.estimated  :   Estimated classes of the learning set
 Prob.estimated     :   n x numberofclasses probability matrix
 */
RcppExport SEXP KernelClassificationkNNlCV(SEXP DistanceMatrix, SEXP classes, SEXP knnlen) {
    
    Rcpp::NumericMatrix Kr(DistanceMatrix);
    Rcpp::IntegerVector yr(classes);
    int knn_len = Rcpp::as<int>(knnlen);
    int n = Kr.nrow();
    
    arma::mat Kc(Kr.begin(), n, n, false); 
    arma::icolvec y(yr.begin(), yr.size(), false);
    int nbclass = arma::max(y);
    
    arma::colvec tknn_lenones = arma::ones<arma::colvec>(knn_len);
    arma::colvec tnzeros = arma::zeros<arma::colvec>(n);
    
    arma::uvec bandwidth_opt(n);
    arma::uvec classes_est(n);
    
    arma::colvec tmp(n);
    arma::rowvec tmp2(nbclass);
    arma::colvec zz(knn_len + 1);
    arma::colvec criterium(knn_len);
    arma::colvec bandwidth(knn_len);
    arma::colvec z(knn_len);
    arma::colvec ytmp(knn_len);
    arma::colvec ktemp(knn_len);
    
    arma::uvec tmp_ind(n);
    arma::uvec norm_order(knn_len);
    arma::u32 index;        
    
    arma::umat Binary(n, nbclass);
    arma::mat yhat_tmp(knn_len, nbclass);
    arma::mat Probhat(n, nbclass);
    for (int g = 1; g <= nbclass; g++) {
        Binary.col(g - 1) = (y == g);
    }
    double denom;
    for (int i = 0; i < n; i++) {
        tmp_ind = sort_index(Kc.col(i));
        tmp = sort(Kc.col(i));
        for (int k = 0; k <= knn_len; k++) {
            zz(k) = tmp(k + 1);
        }
        for (int k = 0; k < knn_len; k++) {
            norm_order(k) = tmp_ind(k + 1);
            z(k) = zz(k);
            bandwidth(k) = 0.5 * (zz(k) + zz(k + 1));
        }
        for (int j = 0; j < knn_len; j++) {
            ktemp = z / bandwidth(j);
            ktemp = tknn_lenones - (ktemp % ktemp);
            ktemp = ktemp % (ktemp >= 0);
            ktemp = ktemp % (ktemp <= 1);
            denom = arma::accu(ktemp);
            for (int g = 0; g < nbclass; g++) {
                for (int k = 0; k < knn_len; k++) {
                    ytmp(k) = Binary(norm_order(k), g);
                }
                yhat_tmp(j, g) = arma::accu(ktemp % ytmp) / denom;
            }
            criterium[j] = arma::accu((yhat_tmp.row(j) - Binary.row(i)) % (yhat_tmp.row(j) - Binary.row(i)));
        }
        criterium.min(index);
        bandwidth_opt(i) = index;
        tmp2 = yhat_tmp.row(index);
        Probhat.row(i) = tmp2;
        tmp2.max(index);
        classes_est(i) = index + 1;
    }
    return Rcpp::List::create(Rcpp::Named("kopt") = Rcpp::wrap(bandwidth_opt), Rcpp::Named("classes.estimated") = Rcpp::wrap(classes_est), Rcpp::Named("Prob.estimated") = Rcpp::wrap(Probhat));
}

/*
 Portet from the R-code written by Ferraty and Vieu. More methods are available
 on their website http://www.math.univ-toulouse.fr/staph/npfda/ 
 
 
 kNN-classification (prediction)
 
 
 INPUT:
 DistanceMatrix     :    n x m - distance matrix
 classes            :    Vector of classes with length n
 knnlen             :    Number of neighbours for bandwidth selection
 
 OUTPUT:
 classes.predicted  :   Vector of predicted classes with length m
 Prob.predicted     :   m x numberofclasses probability matrix
 */
RcppExport SEXP KernelClassificationkNN(SEXP DistanceMatrix, SEXP classes, SEXP neighbours) {
    
    Rcpp::NumericMatrix Kr(DistanceMatrix);
    Rcpp::IntegerVector yr(classes);
    Rcpp::IntegerVector knnr(neighbours);
    int n = Kr.nrow();
    int m = Kr.ncol();
    
    arma::mat Kc(Kr.begin(), n, m, false); 
    arma::mat K = Kc;
    arma::icolvec y(yr.begin(), yr.size(), false);
    int nbclass = arma::max(y);
    arma::ivec knn(knnr.begin(), knnr.size(), false);
    arma::colvec resp;
    
    arma::rowvec h(m);
    arma::rowvec tmp2(nbclass);
    arma::colvec tmp(n);
    arma::u32 index;
    arma::uvec classes_pred(m);
    
    arma::mat KNN = sort(K);
    
    arma::umat Binary(n, nbclass);
    arma::mat Probhat(m, nbclass);
    for (int g = 1; g <= nbclass; g++) {
        Binary.col(g - 1) = (y == g);
    }
    int k;
    for (int i = 0; i < m; i++) {
        tmp = K.col(i);
        tmp.min(index);
        k = knn(index);
        tmp = KNN.col(i);
        h(i) = 0.5 * (tmp(k) + tmp(k + 1));
        K.col(i) = K.col(i) / h(i);
    }    
    K = arma::ones<arma::mat>(n, m) - (K % K);
    K = K % (K >= 0);
    K = K % (K <= 1);
    for (int g = 0; g < nbclass; g++) {
        Probhat.col(g) = arma::trans(K) * Binary.col(g) / arma::trans(arma::sum(K));
    }
    
    for (int i = 0; i < m; i++) {
        tmp2 = Probhat.row(i);
        tmp2.max(index);
        classes_pred(i) = index + 1;
    }
    
    return Rcpp::List::create(Rcpp::Named("classes.predicted") = Rcpp::wrap(classes_pred), Rcpp::Named("Prob.predicted") = Rcpp::wrap(Probhat));
}

