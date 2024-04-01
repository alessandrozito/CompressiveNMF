#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]


//// [[Rcpp::export]]
// void sample_weights_cpp(arma::mat& shape_mat, arma::mat& rate_mat, arma::mat& Theta) {
//   int K = shape_mat.n_rows;
//   int J = shape_mat.n_cols;
//   for(int k = 0; k<K; k++){
//     for(int j = 0; j<J; j++){
//       Theta(k, j) = arma::randg(1,arma::distr_param(shape_mat(k, j), 1/rate_mat(k, j)))(0);
//     }
//   }
// }

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

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
#timesTwo(42)
Sig <- cosmic_data$SBS1
tune_betah(Sig)

*/
