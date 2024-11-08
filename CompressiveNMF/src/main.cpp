#include <RcppArmadillo.h>
using namespace Rcpp;


// [[Rcpp::depends(RcppArmadillo)]]

// // [[Rcpp::export]]
// double randbeta(double par1, double par2) {
//   double x1 = arma::randg(1,arma::distr_param(par1,1.0))(0);
//   double x2 = arma::randg(1,arma::distr_param(par2,1.0))(0);
//   return x1/(x1+x2);
// }

// [[Rcpp::export]]
void sample_weights_cpp(arma::mat& shape_mat, arma::mat& rate_mat, arma::mat& Theta) {
  int K = shape_mat.n_rows;
  int J = shape_mat.n_cols;
  for(int k = 0; k<K; k++){
    for(int j = 0; j<J; j++){
      Theta(k, j) = arma::randg(1,arma::distr_param(shape_mat(k, j), 1/rate_mat(k, j)))(0);
    }
  }
}

// [[Rcpp::export]]
void sample_signatures_cpp(arma::mat& Alpha, arma::mat& R) {
  int I = Alpha.n_rows;
  int K = Alpha.n_cols;
  for(int i = 0; i<I; i++){
    for(int k = 0; k<K; k++){
      R(i, k) = arma::randg(1,arma::distr_param(Alpha(i, k), 1.0))(0);
    }
  }
  // Normalize
  arma::rowvec S = sum(R, 0);
  for(int k = 0; k<K; k++){
    R.col(k) = R.col(k) / S(k);
  }
}


// [[Rcpp::export]]
arma::vec randmult(double n, arma::vec& prob){
  int S = n;
  double rho = 1.0;
  prob = prob/sum(prob); // normalize
  arma::vec out(prob.n_elem);
  for(int i = 0; i < prob.n_elem - 1; i++){
    if(rho != 0){
      out(i) = sum(arma::randu(S) <= prob[i]/rho);
    } else{
      out(i) = 0.0;
    }
    S = S - out(i);
    rho = rho - prob[i];
  }
  out(prob.n_elem - 1) = S;
  return(out);
}


// [[Rcpp::export]]
void sample_mu_ARD_cpp (double a, double a0, double b0, arma::mat& Theta, arma::vec& mu){
  int K = Theta.n_rows;
  int J = Theta.n_cols;
  for(int k = 0; k < K; k++){
    mu(k) = 1/arma::randg(1,arma::distr_param(a0 + J * a, 1/(b0 + a * sum(Theta.row(k)))))(0);
  }

}
// [[Rcpp::export]]
void sample_Y_cpp(arma::mat& X, arma::mat& nonzero_ids, arma::mat& R,
                  arma::mat& Theta, arma::cube& Y) {
  int I = Y.n_rows;
  int J = Y.n_cols;
  int K = Y.n_slices;
  arma::vec q_vec(K);
  arma::vec M(K);
  for(int id = 0; id < nonzero_ids.n_rows; id++){
    int i = nonzero_ids(id, 0)-1;
    int j = nonzero_ids(id, 1)-1;
    q_vec = R.row(i).t() % Theta.col(j);
    // Sample the multinomial
    M = randmult(X(i, j),  q_vec);
    // Allocate the output in Y
    for(int k = 0; k < K; k++){
      Y(i, j, k) = M(k);
    }
  }
}

// [[Rcpp::export]]
arma::cube sample_Y_cpp2(arma::mat& X, arma::mat& nonzero_ids, arma::mat& R,
                  arma::mat& Theta, arma::cube& Yt) {
  int I = Yt.n_rows;
  int J = Yt.n_cols;
  int K = Yt.n_slices;
  arma::cube Y = Yt;
  arma::vec q_vec(K);
  arma::vec M(K);
  for(int id = 0; id < nonzero_ids.n_rows; id++){
    int i = nonzero_ids(id, 0)-1;
    int j = nonzero_ids(id, 1)-1;
    q_vec = R.row(i).t() % Theta.col(j);
    // Sample the multinomial
    M = randmult(X(i, j),  q_vec);
    // Allocate the output in Y
    for(int k = 0; k < K; k++){
      Y(i, j, k) = M(k);
    }
  }
  return(Y);
}
// [[Rcpp::export]]
List CompressiveNMF_cpp(arma::mat& X,
                        arma::mat nonzero_ids,
                        arma::mat R,
                        arma::mat Theta,
                        arma::cube Y,
                        arma::vec mu,
                        arma::mat SignaturePrior,
                        int nsamples,
                        int burnin,
                        double a,
                        double a0,
                        double b0){
  int I = X.n_rows;
  int J = X.n_cols;
  int K = Theta.n_rows;
  // Store the output
  arma::mat shape_mat(K, J);
  arma::mat rate_mat(K, J);
  arma::mat Alpha(I, K);
  arma::cube SIGN(I, K, nsamples);
  arma::cube THETA(K, J, nsamples);
  arma::mat MU(nsamples, K);
  arma::mat sumY1;
  int R_show = (nsamples + burnin)/20;
  // Begin sampling now!
  for(int iter = 0; iter < nsamples + burnin; iter++) {
    // Monitor the number of processed points
    if((iter+1)%std::max(1, R_show)==0) {
      int state = round(100 * (iter+1)/(nsamples + burnin));
      Rprintf("Iteration: %i [%i%%] \n", iter+1,  state);
    }

    // ------------------------------ 1. Sample the latent variables from multinomial
    sample_Y_cpp(X, nonzero_ids, R, Theta, Y);
    // ------------------------------ 2. Sample the weights
    shape_mat = a + sum(Y, 0);
    shape_mat = shape_mat.t();
    for(int j = 0; j < J; j++){
      rate_mat.col(j) =  1 +  a / mu;
    }

    sample_weights_cpp(shape_mat, rate_mat, Theta);
    //------------------------------ 3. Sample the signatures
    sumY1 = sum(Y, 1);
    Alpha = SignaturePrior + sumY1;
    sample_signatures_cpp(Alpha, R);
    //------------------------------ 4. Sample the relavance weights
    sample_mu_ARD_cpp(a, a0, b0, Theta, mu);

    // SAVE THE OUTPUT NOW
    if(iter >= burnin){
      SIGN.slice(iter - burnin) = R;
      THETA.slice(iter - burnin) = Theta;
      MU.row(iter - burnin) = mu.t();
    }
  }
  return(List::create(_["Signatures"] = SIGN,
                      _["Weights"] = THETA,
                      _["Mu"] = MU));

}

/***R
#Y <- sample_Y(X, R = init_pars$R, Theta = init_pars$Theta, nonzero_ids = model_pars$nonzero_ids)
#temp <- CompressiveNMF_cpp(X = X,
#                           nonzero_ids = model_pars$nonzero_ids,
#                           R = init_pars$R,
##                           Theta = init_pars$Theta,
#                           mu = init_pars$mu,
#                           SignaturePrior = model_pars$SignaturePrior,
##                           nsamples = 1000,
##                           burnin = 2000,
#                           a = model_pars$a,
#                           a0 = model_pars$a0,
#                           b0 = model_pars$b0,
#                           Y = Y)


#i <- 4; plot(temp$Mu[, i], type = "l")
#round(colMeans(temp$Mu), 5)

#SigMat <- apply(temp$Signatures, c(1,2), mean)[,round(colMeans(temp$Mu), 5) > 1.5 * 0.001]
#plot_matrix_signature_v2(SigMat, add_cosine_cosmic = FALSE)

#Theta <- init_pars$Theta
#R <- init_pars$R
#mu <- init_pars$mu

#sample_Y_cpp2(X, nonzero_ids = model_pars$nonzero_ids,
#             R = init_pars$R,
#             Theta = init_pars$Theta,
#             Yt = Y)

#shape_mat = apply(Y, c(3, 2), sum) + 1
#rate_mat <- 1 + a/matrix(mu)[, rep(1, 21)]

#sample_weights_cpp(shape_mat = shape_mat, rate_mat = rate_mat, Theta = Theta);

#Alpha <- model_pars$SignaturePrior + apply(Y, c(1, 3), sum)
#sample_signatures_cpp(Alpha = Alpha, R = R)

#sample_mu_ARD_cpp(a, a0, b0, Theta, mu)

#mu
*/




