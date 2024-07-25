// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
double est_entry_FFT(arma::cx_mat V_star_mat, arma::cx_mat V_mat, arma::vec c_vec, int K, int N){
    
    cx_vec Fc = arma::fft(c_vec);
    cx_mat FV_zero = arma::fft(V_mat, c_vec.n_elem);
    cx_mat CV = arma::ifft(FV_zero.each_col() % Fc);
    cx_mat RV = CV.rows(1,N);
    
    double x = arma::norm2est(V_star_mat*RV, 1e-6, 100);

    double output = (1.0/K)*x;

    return output;
}

