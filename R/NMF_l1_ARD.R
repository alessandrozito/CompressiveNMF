# This file replicates Algorithm 2 in Tan and Fevotte (2013) for the l1 automatic
# relevance determination of the number of latent dimension in a non-negative 
# factorization. They assume a Poisson likelihood where the entries are all half-normal
# distributed with inverse-gamma prior. The method uses the MAP estimate.
# The function uses their notation for simplicity. 
# Note that beta = 1 is the equivalent of the KL divergence, with is the Poisson 
# likelihood in non-negative matrix factorization.

gamma_beta_ARD <- function(beta){
  if(beta < 1){
    return(1/(2-beta))
  } else if (beta >= 1 & beta <= 2){
    return(1)
  } else {
    return(1/(beta-1))
  }
} 

NMF_l1_ARD <- function(V, K, a = 5, beta = 1, phi = 1, tau = 1e-7,
                       maxit = 1e6, verbose = TRUE, normalize = TRUE){
  # Set the hyperparameters
  dimF <- nrow(V); dimN <- ncol(V)
  c <- dimF + dimN + a + 1
  gamma <- gamma_beta_ARD(beta)
  mu_v <- mean(V)
  b <- sqrt((a-1) * (a-2) * mu_v/K)
  tol <- Inf
  # Perform the updates as described by their algorithm
  H <- matrix(runif(K * dimN), nrow = K)
  W <- matrix(runif(dimF * K), ncol = K)
  WH <-  (W %*% H)
  lambda <- 1/rgamma(K, a, b)
  iter <- 1
  while((tol > tau) & iter < maxit){
    # Update H
    numH <- crossprod(W, WH^(beta-2) * V)
    denH <- crossprod(W, WH^(beta-1)) + phi/matrix(lambda)[, rep(1, dimN)]
    H <- H * (numH / denH)^gamma
    WH <-  (W %*% H)
    # Check for instability (Inf * 0 -> NaN)
    WH[WH < 1e-308] <- 1e-308
    # Update W
    numW <- tcrossprod(WH^(beta-2) * V, H)
    denW <- tcrossprod(WH^(beta-1), H) + phi/t(lambda)[rep(1, dimF),]
    W <- W * (numW / denW)^gamma
    WH <-  (W %*% H)
    WH[WH < 1e-308] <- 1e-308
    # Update Lambda
    lambda_new <- (colSums(W) + rowSums(H) + b)/c
    tol <- max(abs(lambda_new - lambda)/lambda)
    lambda <- lambda_new
    # Update iterations
    iter <- iter + 1
    if(verbose == TRUE){
      # Every 1000 iterations, print the level of tolerance
      if(iter %% 1000 == 0){
        print(paste0("Iteration: ", iter, " - Tolerance: ", round(tol, 7)))
      }
    }
    
  }
  if(iter >= maxit){
    print("The number of iterations excedeed the maximum.")
  }
  # Exclude the signatures that fall below the threshold
  sign_to_include <- (lambda - b/c)/(b/c) > tau
  W <- W[, sign_to_include]
  H <- H[sign_to_include, ]
  lambda <-lambda[sign_to_include]
  # Normalize the signatures
  if(normalize) {
    normalizing_const <- colSums(W)
    W <- W / t(normalizing_const)[rep(1, dimF),]
    H <- H * matrix(normalizing_const)[, rep(1, dimN)]
  } else {
    normalizing_const <- rep(1, K)
  }
  
  return(list(W = W, H = H,
              lambda = lambda,
              normalizing_const = normalizing_const))
}


