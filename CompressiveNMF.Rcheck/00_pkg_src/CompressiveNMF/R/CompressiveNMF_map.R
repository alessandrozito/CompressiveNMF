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
#'
#' @export
#' @useDynLib CompressiveNMF
CompressiveNMF_map <- function(X,
                               K = 30,
                               S = NULL,
                               alpha = 1,
                               a = 1,
                               epsilon = 0.001,
                               a0 = a * ncol(X) + 1,
                               b0 = a * epsilon * ncol(X),
                               cutoff = 5 * epsilon,
                               tol = 1e-6,
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
    colnames(SignaturePrior) <- c(colnames(S), paste0(paste0("Sig_new", c(1:K))))
  } else {
    colnames(SignaturePrior) <- colnames(S)
  }

  rownames(SignaturePrior) <- rownames(X)
  SignaturePrior <- pmax(SignaturePrior, 1e-10) # correction for Dirichlet

  # Random initialization from the prior
  Ktot <- ncol(SignaturePrior)
  R <- sample_signatures(SignaturePrior)
  Mu <- 1/rgamma(Ktot, a0, b0)
  shape_mat <- matrix(a, nrow = Ktot, ncol = J)
  rate_mat <- as.matrix(a/Mu)[, rep(1, J)]
  Theta <- sample_weights(shape_mat, rate_mat)
  # Run the code in Rcpp
  res <- compute_CompressiveNMF_MAP(X, R, Theta, Mu, SignaturePrior,
                                    a0, b0, a,
                                    maxiter, tol)
  # Postprocess results
  ids <- which(c(res$Mu) > cutoff)
  R <- res$R[, ids]
  colnames(R) <- colnames(SignaturePrior[, ids]); rownames(R) <- rownames(X)
  Theta <- res$Theta[ids, ]
  rownames(Theta) <- colnames(SignaturePrior[, ids]); colnames(Theta) <- colnames(X)
  Mu <- c(res$Mu)[ids]
  names(Mu) <- colnames(SignaturePrior[, ids])

  # Calculate the lkogPosterior
  logpost = calculate_logPosterior(X, R, Theta, Mu, a, a0, b0, SignaturePrior)
  return(list(Signatures = R,
              Theta = Theta,
              Mu = Mu,
              mapOutput = res,
              a0 = a0,
              b0 = b0,
              a = a,
              logPosterior = logpost,
              SignaturePrior = SignaturePrior,
              MutMatrix = X))
}





