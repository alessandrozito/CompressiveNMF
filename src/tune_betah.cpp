#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export]]
double MedianCosineSimilarity(double beta, arma::vec Sig) {
  arma::mat R(1000, Sig.n_elem);
  for(int i = 0; i< Sig.n_elem; i++){
    R.col(i) = arma::randg(1000, arma::distr_param(beta * Sig(i), 1.0));
  }
  arma::mat Y = R * Sig; // matrix product
  arma::mat tot_R = arma::sqrt(arma::diagvec(R * arma::trans(R)));
  arma::mat tot_Sig = arma::sqrt(arma::trans(Sig) * Sig);
  arma::mat res = arma::median(Y/(tot_R * tot_Sig(0,0)));
  return res(0,0);
}

double tune_betah(arma::vec Sig, double target_sim = 0.975, int n_points = 1000) {
  arma::vec betas = arma::logspace(10, 5000, n_points);
  double best_beta = betas(0);
  double distance = std::abs(MedianCosineSimilarity(best_beta, Sig) - target_sim);
  for(int i = 1; i < n_points; i++){
    double new_distance = std::abs(MedianCosineSimilarity(betas(i), Sig) - target_sim);
    if(new_distance < distance){
      distance = new_distance;
      best_beta = betas(i);
    }
  }
  return(best_beta);
}


