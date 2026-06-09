#' Maximum-a-posteriori for the CompressiveNMF algorithm
#'
#' @param X Matrix of mutations (rows are mutations, columns are patients)
#' @param K Upper bound to the number of signatures
#' @param S Matrix of informative prior for the signatures
#' @param alpha Parameter for the Dirichlet prior in the signature
#' @param a Parameter for the gamma prior over the loadings
#' @param epsilon Parameter controlling the mean of the Compressive hyperprior
#' @param a0 Shape parameter for the inverse gamma compressive hyperprior. Default is \code{a * ncol(X) + 1}
#' @param b0 Rate parameter for the inverse gamma compressive hyperprior. Default is \code{a * epsilon * ncol(X)}
#' @param cutoff Cutoff for excluded signatures
#' @param tol Tolerance for convergence. Set to \code{1e-6}
#' @param maxiter Maximum number of iterations for the algorithm. Set to  \code{1e6}
#' @importFrom Rcpp sourceCpp
#' @export
#' @useDynLib CompressiveNMF
CompressiveNMF_MultiStudy <- function(X,
                                      cohorts_num,
                                       K = 20,
                                       S = NULL,
                                       alpha = 1,
                                       a = 1,
                                       epsilon = 0.001,
                                       cutoff = 5 * epsilon,
                                       hierarchy = FALSE,
                                       compressive = TRUE,
                                       tol = 1e-7,
                                       a0 = epsilon * 1,
                                       b0 = 1,
                                       maxiter = 1e6){

  # Check if the matrix X is made of integers
  if(any(round(X) != X)){
    stop("The matrix X has to be made of integers! Consider using round(X).")
  }
  # Matrix dimensions
  I <- nrow(X); J <- ncol(X)

  # Check if the matrix S is correctly specified
  if(!is.null(S)){
    if(!is.matrix(S)){
      stop("The hyperparameter S must be a matrix!")
    }
    if(any(S<1)){
      warning("Maximum-a-posteriori not guaranteed when any(S<1) is TRUE")
    }
  }
  if(alpha <= 0){
    stop("Must have alpha > 0")
  }
  if(a <= 0){
    stop("Must have a > 0")
  }
  if(K < 0 | round(K) != K){
    stop("Must have K >= 0 and integer.")
  }
  if(alpha < 1 | a < 1) {
    warning("Maximum-a-posteriori not guaranteed when alpha < 1 or a < 1")
  }

  # Fix the prior matrix for the signatures
  SignaturePrior <- cbind(S, matrix(alpha, nrow = I, ncol = K))
  if(K > 0){
    colnames(SignaturePrior) <- c(colnames(S), paste0(paste0("Sig", sprintf("%02d", 1:K))))
  } else {
    colnames(SignaturePrior) <- colnames(S)
  }

  rownames(SignaturePrior) <- rownames(X)
  SignaturePrior <- pmax(SignaturePrior, 1e-10) # correction for Dirichlet

  # Random initialization from the prior
  Ktot <- ncol(SignaturePrior)
  Stot <- length(unique(cohorts_num))
  Mu_start <- 1/matrix(rgamma(n = Stot * Ktot, 1, 1), ncol = Stot)
  R <- sample_signatures(SignaturePrior)
  shape_mat <- matrix(a, nrow = Ktot, ncol = J)
  rate_mat <- as.matrix(a/Mu_start[, 1])[, rep(1, J)]
  Theta <- sample_weights(shape_mat, rate_mat)
  # Run the code in Rcpp
  res <- compute_CompressiveNMF_MAP_MultiStudy(X = X, R_start = R,
                                              Theta_start = Theta,
                                              Mu_start = Mu_start,
                                              hierarchy = hierarchy,
                                              compressive = compressive,
                                              SigPrior = SignaturePrior,
                                              cohorts_num = cohorts_num,
                                              a0 = a0, b0 = b0, a = a,
                                              maxiter = maxiter,
                                              tol = tol)
  # Postprocess results
  #ids <- which(c(res$Mu) > cutoff)
  R <- as.matrix(res$R)
  colnames(R) <- colnames(SignaturePrior); rownames(R) <- rownames(X)
  Theta <- res$Theta
  if(!is.matrix(Theta)){
    Theta <- as.matrix(t(Theta))
  }
  rownames(Theta) <- colnames(SignaturePrior); colnames(Theta) <- colnames(X)
  Mu <- res$Mu
  rownames(Mu) <- colnames(SignaturePrior)
  colnames(Mu) <- unique(cohorts_num)
  Mu_low <- c(res$Mu_low)
  names(Mu_low) <- colnames(SignaturePrior)
  # Calculate the lkogPosterior
  #logpost = calculate_logPosterior(X, R, Theta, Mu, a, a0, b0, SignaturePrior)
  return(list(Signatures = R,
              Theta = Theta,
              Mu = Mu,
              Mu_low = Mu_low,
              mapOutput = res,
              a0 = a0,
              b0 = b0,
              a = a,
              compressive = compressive,
              hierarchy = hierarchy,
              #logPosterior = logpost,
              SignaturePrior = SignaturePrior,
              MutMatrix = X))
}


