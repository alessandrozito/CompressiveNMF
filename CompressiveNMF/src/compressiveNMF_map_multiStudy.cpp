#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
double eval_logPosterior_MultiStudy(arma::mat &X,
                                    arma::mat &R,
                                    arma::mat &Theta,
                                    arma::mat &Mu,
                                    arma::vec &Mu_low,
                                    arma::mat &SigPrior,
                                    arma::uvec &cohorts_num,
                                    arma::vec &Js,
                                    double a,
                                    double a0,
                                    double b0,
                                    bool hierarchy,
                                    bool compressive,
                                    double lambda = 0.0,
                                    double delta  = 1.0) {
  int J = X.n_cols;
  int K = R.n_cols;
  int S = Mu.n_cols;

  arma::mat Lambda = R * Theta;
  double logpost = 0.0;

  logpost += arma::accu(X % log(Lambda) - Lambda);
  logpost += arma::accu((SigPrior - 1) % log(R));

  for (int j = 0; j < J; j++) {
    int s = cohorts_num(j);
    for (int k = 0; k < K; k++) {
      logpost += a * std::log(a / Mu(k, s))
      + (a - 1) * std::log(Theta(k, j))
      - a * Theta(k, j) / Mu(k, s);
    }
  }

  if (hierarchy) {
    for (int s = 0; s < S; s++) {
      for (int k = 0; k < K; k++) {
        logpost += (a * Js(s) + 1) * std::log(Mu_low(k) * a * Js(s))
        - (a * Js(s) + 2.0) * std::log(Mu(k, s))
        - a * Js(s) * Mu_low(k) / Mu(k, s);
      }
    }
    logpost += arma::accu((a0 - 1) * log(Mu_low) - b0 * Mu_low);
  } else {
    if (compressive) {
      for (int s = 0; s < S; s++) {
        a0 = a * Js[s] + 1;
        b0 = 0.01 * a * Js[s];
        logpost += arma::accu(-(a0 + 1) * log(Mu.col(s)) - b0 / Mu.col(s));
      }
    } else {
      logpost += arma::accu(-(a0 + 1) * log(Mu) - b0 / Mu);
    }
  }

  if (lambda > 0.0) {
    double ld, sign_val;
    arma::log_det(ld, sign_val, R.t() * R + delta * arma::eye(K, K));
    logpost -= lambda * ld;
  }

  return logpost;
}


// [[Rcpp::export]]
List compute_CompressiveNMF_MAP_MultiStudy(arma::mat X,
                                           arma::mat R_start,
                                           arma::mat Theta_start,
                                           arma::mat Mu_start,
                                           arma::mat &SigPrior,
                                           arma::uvec &cohorts_num,
                                           double a0,
                                           double b0,
                                           double a,
                                           int maxiter,
                                           double tol,
                                           bool hierarchy,
                                           bool compressive) {

  int K = R_start.n_cols;
  int I = X.n_rows;
  int J = X.n_cols;
  int S = Mu_start.n_cols;

  arma::mat R      = R_start;
  arma::mat Theta  = Theta_start;
  arma::mat Mu     = Mu_start;
  arma::mat Mu_new = Mu_start;
  arma::mat Mu_all(K, J);
  arma::vec Js(S, arma::fill::zeros);
  arma::vec Mu_low(K);
  arma::vec trace;

  for (int j = 0; j < J; j++) Js(cohorts_num(j)) += 1.0;
  for (int k = 0; k < K; k++)
    Mu_low(k) = (a0 + a * J - 1) / (b0 + a * arma::sum(Js / Mu_new.row(k).t()));

  double lp = eval_logPosterior_MultiStudy(X, R, Theta, Mu, Mu_low,
                                           SigPrior, cohorts_num, Js,
                                           a, a0, b0, hierarchy, compressive);
  double maxdiff = std::numeric_limits<double>::infinity();
  int it = 0;

  for (int iter = 0; iter < maxiter; iter++) {

    if ((iter + 1) % 1000 == 0)
      Rprintf("Iteration %i - logpost rel diff %.10f \n", iter + 1, maxdiff);

    // Update R
    arma::mat R_upd = (X / (R * Theta)) * Theta.t();
    R = SigPrior - 1 + R % R_upd;
    R.elem(arma::find(R < arma::datum::eps)).fill(arma::datum::eps);
    R = arma::normalise(R, 1, 0);

    // Update Theta
    arma::mat Theta_upd = R.t() * (X / (R * Theta));
    for (int j = 0; j < J; j++) Mu_all.col(j) = Mu_new.col(cohorts_num(j));
    Theta = (Mu_all / (a + Mu_all)) % (a - 1 + Theta % Theta_upd);
    Theta.elem(arma::find(Theta < arma::datum::eps)).fill(arma::datum::eps);

    // Update Mu
    arma::mat theta_sums(K, S, arma::fill::zeros);
    for (int j = 0; j < J; j++) theta_sums.col(cohorts_num(j)) += Theta.col(j);

    if (hierarchy) {
      for (int s = 0; s < S; s++)
        Mu_new.col(s) = (a * theta_sums.col(s) + Mu_low * a * Js[s]) / (2 * a * Js[s] + 2.0);
      for (int k = 0; k < K; k++)
        Mu_low(k) = (a0 + a * J - 1) / (b0 + a * arma::sum(Js / Mu_new.row(k).t()));
    } else {
      if (compressive) {
        for (int s = 0; s < S; s++) {
          double a0_s = a * Js[s] + 1;
          double b0_s = 0.01 * a * Js[s];
          Mu_new.col(s) = (a * theta_sums.col(s) + b0_s) / (a * Js[s] + a0_s + 1.0);
        }
      } else {
        for (int s = 0; s < S; s++)
          Mu_new.col(s) = (a * theta_sums.col(s) + b0) / (a * Js[s] + a0 + 1.0);
      }
    }

    double lp_new = eval_logPosterior_MultiStudy(X, R, Theta, Mu_new, Mu_low,
                                                 SigPrior, cohorts_num, Js,
                                                 a, a0, b0, hierarchy, compressive);
    maxdiff = std::abs((lp_new - lp) / lp);
    lp = lp_new;
    Mu = Mu_new;
    it++;
    if (maxdiff < tol) break;
    if ((iter + 1) % 100 == 0)
      trace = arma::join_vert(trace, arma::vec({lp_new}));
  }

  return List::create(_["Theta"]   = Theta,
                      _["R"]       = R,
                      _["Mu"]      = Mu,
                      _["Mu_low"]  = Mu_low,
                      _["iter"]    = it,
                      _["maxdiff"] = maxdiff,
                      _["logpost"] = lp,
                      _["logpost_trace"] = arma::join_vert(trace, arma::vec({lp})));
}


// ─── Min-Vol helpers ──────────────────────────────────────────────────────────

static double eval_minvol_obj(const arma::mat &X,
                              const arma::mat &R,
                              const arma::mat &Theta,
                              double lambda,
                              double delta) {
  arma::mat Lambda = R * Theta;
  double kl = -arma::accu(X % arma::log(Lambda)) + arma::accu(Lambda);
  double ld, sign_val;
  arma::log_det(ld, sign_val, R.t() * R + delta * arma::eye(R.n_cols, R.n_cols));
  return kl + lambda * ld;
}

static void normalize_cols(arma::mat &R, arma::mat &Theta) {
  arma::rowvec col_norms = arma::clamp(arma::sum(R, 0), arma::datum::eps, arma::datum::inf);
  R     = R     * arma::diagmat(1.0 / col_norms);
  Theta = arma::diagmat(col_norms) * Theta;
}


// [[Rcpp::export]]
List compute_MinVol_MAP_MultiStudy(arma::mat X,
                                   arma::mat R_start,
                                   arma::mat Theta_start,
                                   arma::mat Mu_start,
                                   arma::mat &SigPrior,
                                   arma::uvec &cohorts_num,
                                   double a0,
                                   double b0,
                                   double a,
                                   double lambda,
                                   double delta,
                                   int maxiter,
                                   double tol,
                                   bool hierarchy,
                                   bool compressive) {

  int K = R_start.n_cols;
  int I = X.n_rows;
  int J = X.n_cols;
  int S = Mu_start.n_cols;

  arma::mat R      = R_start;
  arma::mat Theta  = Theta_start;
  arma::mat Mu     = Mu_start;
  arma::mat Mu_new = Mu_start;
  arma::mat Mu_all(K, J);
  arma::vec Js(S, arma::fill::zeros);
  arma::vec Mu_low(K);
  arma::vec trace;

  for (int j = 0; j < J; j++) Js(cohorts_num(j)) += 1.0;
  for (int k = 0; k < K; k++)
    Mu_low(k) = (a0 + a * J - 1) / (b0 + a * arma::sum(Js / Mu_new.row(k).t()));

  double gamma   = 1.0;
  double lp      = eval_logPosterior_MultiStudy(X, R, Theta, Mu, Mu_low,
                                                SigPrior, cohorts_num, Js,
                                                a, a0, b0, hierarchy, compressive,
                                                lambda, delta);
  double maxdiff = std::numeric_limits<double>::infinity();
  int    it      = 0;

  for (int iter = 0; iter < maxiter; iter++) {

    if ((iter + 1) % 1000 == 0)
      Rprintf("Iteration %i - logpost rel diff %.10f \n", iter + 1, maxdiff);

    // ── Update R (Algorithm 1, Table II) ─────────────────────────────────────
      {
        arma::mat Y    = arma::inv(R.t() * R + delta * arma::eye(K, K));  // LU, no dpotrf_
        arma::mat Ypos = arma::clamp( Y, 0.0, arma::datum::inf);
        arma::mat Yneg = arma::clamp(-Y, 0.0, arma::datum::inf);

        arma::mat Lambda = R * Theta;
        arma::mat JHt    = arma::repmat(arma::sum(Theta, 1).t(), I, 1);
        arma::mat P      = JHt - 4.0 * lambda * (R * Yneg);
        arma::mat Q      = (X / Lambda) * Theta.t();
        arma::mat Sv     = R * (Ypos + Yneg);
        Sv.elem(arma::find(Sv < arma::datum::eps)).fill(arma::datum::eps);

        arma::mat R_new = R % (arma::sqrt(arma::square(P) + 8.0 * lambda * Sv % Q) - P)
          / (4.0 * lambda * Sv);
        R_new.elem(arma::find(R_new < arma::datum::eps)).fill(arma::datum::eps);

        double F_curr         = eval_minvol_obj(X, R, Theta, lambda, delta);
        arma::mat R_gamma     = R_new;
        arma::mat Theta_gamma = Theta;
        normalize_cols(R_gamma, Theta_gamma);

        int ls_it = 0;
        while (eval_minvol_obj(X, R_gamma, Theta_gamma, lambda, delta) > F_curr
                 && ls_it < 100) {
          gamma       *= 0.8;
          R_gamma      = (1.0 - gamma) * R + gamma * R_new;
          Theta_gamma  = Theta;
          R_gamma.elem(arma::find(R_gamma < arma::datum::eps)).fill(arma::datum::eps);
          normalize_cols(R_gamma, Theta_gamma);
          ls_it++;
        }
        R     = R_gamma;
        Theta = Theta_gamma;
        gamma = std::min(1.0, gamma * 1.2);
      }

      // ── Update Theta (unchanged) ─────────────────────────────────────────────
      arma::mat Theta_upd = R.t() * (X / (R * Theta));
      for (int j = 0; j < J; j++) Mu_all.col(j) = Mu_new.col(cohorts_num(j));
      Theta = (Mu_all / (a + Mu_all)) % (a - 1 + Theta % Theta_upd);
      Theta.elem(arma::find(Theta < arma::datum::eps)).fill(arma::datum::eps);

      // ── Update Mu (unchanged) ────────────────────────────────────────────────
      arma::mat theta_sums(K, S, arma::fill::zeros);
      for (int j = 0; j < J; j++) theta_sums.col(cohorts_num(j)) += Theta.col(j);

      if (hierarchy) {
        for (int s = 0; s < S; s++)
          Mu_new.col(s) = (a * theta_sums.col(s) + Mu_low * a * Js[s]) / (2 * a * Js[s] + 2.0);
        for (int k = 0; k < K; k++)
          Mu_low(k) = (a0 + a * J - 1) / (b0 + a * arma::sum(Js / Mu_new.row(k).t()));
      } else {
        if (compressive) {
          for (int s = 0; s < S; s++) {
            double a0_s = a * Js[s] + 1;
            double b0_s = 0.01 * a * Js[s];
            Mu_new.col(s) = (a * theta_sums.col(s) + b0_s) / (a * Js[s] + a0_s + 1.0);
          }
        } else {
          for (int s = 0; s < S; s++)
            Mu_new.col(s) = (a * theta_sums.col(s) + b0) / (a * Js[s] + a0 + 1.0);
        }
      }

      // ── Convergence ──────────────────────────────────────────────────────────
      double lp_new = eval_logPosterior_MultiStudy(X, R, Theta, Mu_new, Mu_low,
                                                   SigPrior, cohorts_num, Js,
                                                   a, a0, b0, hierarchy, compressive,
                                                   lambda, delta);
      maxdiff = std::abs((lp_new - lp) / lp);
      lp  = lp_new;
      Mu  = Mu_new;
      it++;
      if (maxdiff < tol) break;
      if ((iter + 1) % 100 == 0)
        trace = arma::join_vert(trace, arma::vec({lp_new}));
  }

  return List::create(_["Theta"]   = Theta,
                      _["R"]       = R,
                      _["Mu"]      = Mu,
                      _["Mu_low"]  = Mu_low,
                      _["iter"]    = it,
                      _["maxdiff"] = maxdiff,
                      _["logpost"] = lp,
                      _["logpost_trace"] = arma::join_vert(trace, arma::vec({lp})));
}

// #include <RcppArmadillo.h>
// using namespace Rcpp;
//
// // [[Rcpp::depends(RcppArmadillo)]]
//
// // Log-posterior for the hierarchical multi-study compressive NMF.
// //
// // hierarchy = true:
// //   Theta_{kj} ~ Gamma(a, a/Mu_{ks})
// //   Mu_{ks}    ~ IngGamma(a * J_s + 1, Mu_low_k * a * J_s)
// //   Mu_low_k   ~ Gamma(a0, b0)
// //
// // hierarchy = false:
// //   Theta_{kj} ~ Gamma(a, a/Mu_{ks})
// //   Mu_{ks}    ~ InvGamma(a0, b0)
// //
// // Mu     : K x S matrix
// // Mu_low : K-vector (only used when hierarchy = true)
// // cohorts_num : J-vector mapping each sample to its study (0-indexed)
// // [[Rcpp::export]]
// double eval_logPosterior_MultiStudy(arma::mat &X,
//                                     arma::mat &R,
//                                     arma::mat &Theta,
//                                     arma::mat &Mu,
//                                     arma::vec &Mu_low,
//                                     arma::mat &SigPrior,
//                                     arma::uvec &cohorts_num,
//                                     arma::vec &Js,
//                                     double a,
//                                     double a0,
//                                     double b0,
//                                     bool hierarchy,
//                                     bool compressive) {
//   int J = X.n_cols;
//   int K = R.n_cols;
//   int S = Mu.n_cols;
//
//   arma::mat Lambda = R * Theta;
//   double logpost = 0.0;
//
//   // Poisson log-likelihood
//   logpost += arma::accu(X % log(Lambda) - Lambda);
//
//   // Dirichlet prior on R
//   logpost += arma::accu((SigPrior - 1) % log(R));
//
//   // Gamma prior on Theta_{kj} ~ Gamma(a, a/Mu_{ks})
//   for (int j = 0; j < J; j++) {
//     int s = cohorts_num(j);
//     for (int k = 0; k < K; k++) {
//       logpost += a * std::log(a / Mu(k, s))
//       + (a - 1) * std::log(Theta(k, j))
//       - a * Theta(k, j) / Mu(k, s);
//     }
//   }
//
//   if (hierarchy) {
//     // Inverse gamma prior on Mu_{ks} ~ InvGamma(a * J_s + 1, Mu_low_k * a * J_s)
//     for (int s = 0; s < S; s++) {
//       for (int k = 0; k < K; k++) {
//         logpost += (a * Js(s) + 1) * std::log(Mu_low(k) * a * Js(s))
//         - (a * Js(s) + 2.0) * std::log(Mu(k, s))
//         - a * Js(s) * Mu_low(k) / Mu(k, s);
//       }
//     }
//     // InvGamma prior on Mu_low_k ~ Gamma(a0, b0)
//     logpost += arma::accu((a0 - 1) * log(Mu_low) - b0 * Mu_low);
//   } else {
//     // InvGamma prior on Mu_{ks} ~ InvGamma(a0, b0)
//     if (compressive) {
//       for(int s = 0; s < S; s++){
//         a0 = a * Js[s] + 1;
//         b0 = 0.01 * a * Js[s];
//         logpost += arma::accu(-(a0 + 1) * log(Mu.col(s)) - b0 / Mu.col(s));
//       }
//     } else {
//       logpost += arma::accu(-(a0 + 1) * log(Mu) - b0 / Mu);
//     }
//   }
//
//   return logpost;
// }
//
//
// // [[Rcpp::export]]
// List compute_CompressiveNMF_MAP_MultiStudy(arma::mat X,
//                                            arma::mat R_start,
//                                            arma::mat Theta_start,
//                                            arma::mat Mu_start,
//                                            arma::mat &SigPrior,
//                                            arma::uvec &cohorts_num,
//                                            double a0,
//                                            double b0,
//                                            double a,
//                                            int maxiter,
//                                            double tol,
//                                            bool hierarchy,
//                                            bool compressive) {
//
//   int K = R_start.n_cols;
//   int I = X.n_rows;
//   int J = X.n_cols;
//   int S = Mu_start.n_cols;
//
//   arma::mat R      = R_start;
//   arma::mat Theta  = Theta_start;
//   arma::mat Mu     = Mu_start;
//   arma::mat Mu_new = Mu_start;
//   arma::mat Mu_all(K, J);
//   arma::vec Js(S, arma::fill::zeros);
//   arma::vec Mu_low(K);
//   arma::vec trace;
//
//   for (int j = 0; j < J; j++) {
//     Js(cohorts_num(j)) += 1.0;
//   }
//   for (int k = 0; k < K; k++) {
//     Mu_low(k) = (a0 + a * J - 1) / (b0 + a * arma::sum(Js / Mu_new.row(k).t()));
//   }
//
//   arma::mat R_upd(I, K);
//   arma::mat Theta_upd(K, J);
//
//   double lp = eval_logPosterior_MultiStudy(X, R, Theta, Mu, Mu_low,
//                                            SigPrior, cohorts_num,Js,
//                                            a, a0, b0, hierarchy, compressive);
//   double maxdiff = std::numeric_limits<double>::infinity();
//   int R_show = 1000;
//   int it = 0;
//
//   for (int iter = 0; iter < maxiter; iter++) {
//
//     if ((iter + 1) % R_show == 0) {
//       Rprintf("Iteration %i - logpost rel diff %.10f \n", iter + 1, maxdiff);
//     }
//
//     // Update R
//     R_upd = (X / (R * Theta)) * Theta.t();
//     R = SigPrior - 1 + R % R_upd;
//     R.elem(arma::find(R < arma::datum::eps)).fill(arma::datum::eps);
//     R = arma::normalise(R, 1, 0);
//
//     // Update Theta
//     Theta_upd = R.t() * (X / (R * Theta));
//     for (int j = 0; j < J; j++) {
//       Mu_all.col(j) = Mu_new.col(cohorts_num(j));
//     }
//     Theta = (Mu_all / (a + Mu_all)) % (a - 1 + Theta % Theta_upd);
//     Theta.elem(arma::find(Theta < arma::datum::eps)).fill(arma::datum::eps);
//
//     // Update Mu
//     arma::mat theta_sums(K, S, arma::fill::zeros);
//     for (int j = 0; j < J; j++) {
//       int s = cohorts_num(j);
//       theta_sums.col(s) += Theta.col(j);
//     }
//
//     if (hierarchy) {
//       for (int s = 0; s < S; s++) {
//         Mu_new.col(s) = (a * theta_sums.col(s) + Mu_low * a * Js[s]) / (2 * a * Js[s] + 2.0);
//       }
//       for (int k = 0; k < K; k++) {
//         Mu_low(k) = (a0 + a * J - 1) / (b0 + a * arma::sum(Js / Mu_new.row(k).t()));
//       }
//     } else {
//       if (compressive) {
//         for (int s = 0; s < S; s++) {
//           a0 = a * Js[s] + 1;
//           b0 = 0.01 * a * Js[s];
//           Mu_new.col(s) = (a * theta_sums.col(s) + b0) / (a * Js[s] + a0 + 1.0);
//         }
//       } else {
//         for (int s = 0; s < S; s++) {
//           Mu_new.col(s) = (a * theta_sums.col(s) + b0) / (a * Js[s] + a0 + 1.0);
//         }
//       }
//     }
//
//     // Convergence: relative change in log-posterior
//     double lp_new = eval_logPosterior_MultiStudy(X, R, Theta, Mu_new, Mu_low,
//                                                  SigPrior, cohorts_num, Js,
//                                                  a, a0, b0, hierarchy, compressive);
//     maxdiff = std::abs((lp_new - lp) / lp);
//     lp = lp_new;
//     Mu = Mu_new;
//
//     it += 1;
//     if (maxdiff < tol) break;
//
//     if ((iter + 1) % 100 == 0) {
//       trace = arma::join_vert(trace, arma::vec({lp_new}));
//     }
//   }
//
//   return List::create(_["Theta"]   = Theta,
//                       _["R"]       = R,
//                       _["Mu"]      = Mu,
//                       _["Mu_low"]  = Mu_low,
//                       _["iter"]    = it,
//                       _["maxdiff"] = maxdiff,
//                       _["logpost"] = lp,
//                       _["logpost_trace"] = arma::join_vert(trace, arma::vec({lp})));
// }
