#include <RcppArmadillo.h>
#include <math.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

//[[Rcpp::export]]
double eval_logPosterior(arma::mat &X,
                         arma::mat &R,
                         arma::mat &Theta,
                         arma::vec &Mu,
                         double a,
                         double a0, double b0, arma::mat &SigPrior){
  int I = X.n_rows;
  int J = X.n_cols;
  int K = R.n_cols;
  arma::mat Lambda = R * Theta;
  double logpost = 0.0;
  // Poisson loglikelihood
  logpost += arma::accu(X % log(Lambda) - Lambda);
  // Dirichlet prior
  logpost += arma::accu((SigPrior - 1) % log(R));
  // Gamma prior
  arma::vec a_over_Mu = a * log(a/Mu);
  for(int j = 0; j < J; j++){
    for(int k = 0; k < K; k++){
      logpost += a_over_Mu(k) + (a-1) * log(Theta(k, j)) -
        a * Theta(k, j) / Mu(k);
    }
  }
  // Inverse gamma hyperprior
  logpost += arma::accu(- (a0 + 1) * log(Mu) - b0/Mu);
  return logpost;
}

// [[Rcpp::export]]
List compute_CompressiveNMF_MAP(arma::mat X,
                                arma::mat R_start,
                                arma::mat Theta_start,
                                arma::vec Mu_start,
                                arma::mat &SigPrior,
                                double a0,
                                double b0,
                                double a,
                                int maxiter,
                                double tol,
                                bool use_logpost_for_convergence = false) {

  // Model parameters
  int K = R_start.n_cols;
  int I = X.n_rows;
  int J = X.n_cols;
  // Initialize quantities
  arma::mat R = R_start;
  arma::mat Theta = Theta_start;
  arma::vec Mu = Mu_start;
  arma::vec Mu_new = Mu_start;
  arma::mat Mu_all(K, J);
  // Updated quantities for multiplicative rules
  arma::mat R_upd(I, K);
  arma::mat Theta_upd(K, J);

  double maxdiff = 10;
  int R_show = 1000;
  int it = 0;
  for(int iter = 0; iter < maxiter; iter++){
    if((iter+1)%R_show==0) {
      Rprintf("Iteration %i - diff %.10f \n", iter+1, maxdiff);
    }
    // Update R
    R_upd = (X / (R * Theta)) * Theta.t();
    R = arma::normalise(SigPrior - 1 + R % R_upd, 1, 0);
    R.elem(arma::find(R < arma::datum::eps)).fill(arma::datum::eps);
    R = arma::normalise(R, 1, 0);
    // Update Theta
    Theta_upd = R.t() * (X / (R * Theta));
    for(int j = 0; j < J; j++){
      Mu_all.col(j) = Mu_new;
    }
    Theta = (Mu_all / (a + Mu_all)) % (a - 1 + Theta % Theta_upd);
    Theta.elem(arma::find(Theta < arma::datum::eps)).fill(arma::datum::eps);
    // Update Mu
    Mu_new = (a * sum(Theta, 1) + b0)/ (a * J + a0 + 1);
    // Evaluate the difference on average with relevance weights
    maxdiff = arma::abs(Mu_new/Mu - 1).max();
    Mu = Mu_new;
    if(iter >= maxiter || maxdiff < tol) {
      break;
    }
    it += 1;
  }
  return List::create(_["Theta"] = Theta,
                      _["R"] = R,
                      _["Mu"] = Mu,
                      _["iter"] = it,
                      _["maxdiff"]= maxdiff);
}



