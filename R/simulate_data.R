# This function simulates the data from 
simulate_data <- function(J = 100, cosmic_sig = c("SBS1", "SBS2", "SBS5", "SBS13"), 
                          K_new = 0, alpha = 0.25, theta = 100, overdispersion = 0){
  
  # Generate the signatures
  cosmic_data <- CompressiveNMF::COSMIC_SBS96_GRCh37_v3.4
  Rmat_random <- t(LaplacesDemon::rdirichlet(n = K_new, alpha = rep(alpha, 96)))
  if(K_new > 0) colnames(Rmat_random) <- paste0("SBSnew", 1:K_new)
  Rmat_cos <- as.matrix(cosmic_data[, cosmic_sig])
  Rmat <- as.matrix(cbind(Rmat_cos, Rmat_random))
  rownames(Rmat) <- cosmic_data$Channel
  
  # Generate the weights
  K <- ncol(Rmat)
  exposures <- rgamma(K, theta, 1)
  Theta <- matrix(rgamma(K * J, 0.5, 0.5), ncol = J, nrow = K)
  Theta <- apply(Theta, 2, function(x) x * exposures)
  rownames(Theta) <- colnames(Rmat)
  
  # Generate the counts
  Lambda <- Rmat %*% Theta
  X <- matrix(rnbinom(length(Lambda), size = 1/overdispersion, mu = c(Lambda)), nrow = 96, ncol = J)
  rownames(X) <- rownames(Rmat) 
  return(list(X = X, Rmat = Rmat, Theta = Theta))
}


simulate_data_sparse <- function(J = 100, 
                                 cosmic_sig = c("SBS1", "SBS2", "SBS5", "SBS13"), 
                                 K_new = 0, 
                                 alpha = 0.25, 
                                 theta = 100, 
                                 overdispersion = 0, 
                                 pi0 = 0.1){
  
  # Generate the signatures
  cosmic_data <- CompressiveNMF::COSMIC_SBS96_GRCh37_v3.4
  Rmat_random <- t(LaplacesDemon::rdirichlet(n = K_new, alpha = rep(alpha, 96)))
  if(K_new > 0) colnames(Rmat_random) <- paste0("SBSnew", 1:K_new)
  Rmat_cos <- as.matrix(cosmic_data[, cosmic_sig])
  Rmat <- as.matrix(cbind(Rmat_cos, Rmat_random))
  rownames(Rmat) <- cosmic_data$Channel
  
  # Generate the weights
  K <- ncol(Rmat)
  exposures <- rgamma(K, theta, 1)
  Theta <- matrix(rgamma(K * J, 0.5, 0.5), ncol = J, nrow = K)
  Theta <- apply(Theta, 2, function(x) x * exposures)
  rownames(Theta) <- colnames(Rmat)
  
  # Zero-out the weights at random with probability pi0
  is_Theta_zero <- matrix(rbinom(K * J, 1, pi0), ncol = J, nrow = K)
  if(sum(is_Theta_zero) > 0) {
    array_zeros <- which(is_Theta_zero == 1, arr.ind = TRUE)
    for(r in 1:nrow(array_zeros)){
      Theta[array_zeros[r, 1], array_zeros[r, 2]] <- rgamma(1, 0.1, 0.5)
    }
  }
  
  # Generate the counts
  Lambda <- Rmat %*% Theta
  X <- matrix(rnbinom(length(Lambda), size = 1/overdispersion, mu = c(Lambda)), nrow = 96, ncol = J)
  rownames(X) <- rownames(Rmat) 
  
  # Check if all patients have at least one mutation. If not, 
  # substitute with a 1 in the entry has the highest lambda value
  mean_X <- colMeans(X)
  if(any(mean_X == 0)){
    id_X_zero <- which(mean_X == 0)
    for(j in id_X_zero){
      id_max_lambda <- which.max(Lambda[, j])
      X[id_max_lambda, j] <- 1
    }
  }
  
  return(list(X = X, Rmat = Rmat, Theta = Theta))
}


simulate_data_Indels <- function(J = 100, 
                                 cosmic_sig = c("ID1", "ID2", "ID8", "ID9"), 
                                 K_new = 0, alpha = 0.05, theta = 50, overdispersion = 0){
  
  # Generate the signatures
  cosmic_data <- CompressiveNMF::COSMIC_v3.4_ID83_GRCh37
  Rmat_random <- t(LaplacesDemon::rdirichlet(n = K_new, alpha = rep(alpha, nrow(cosmic_data))))
  if(K_new > 0) colnames(Rmat_random) <- paste0("IDnew", 1:K_new)
  Rmat_cos <- as.matrix(cosmic_data[, cosmic_sig])
  Rmat <- as.matrix(cbind(Rmat_cos, Rmat_random))
  rownames(Rmat) <- rownames(cosmic_data)
  
  # Generate the weights
  K <- ncol(Rmat)
  exposures <- rgamma(K, theta, 1)
  Theta <- matrix(rgamma(K * J, 0.5, 0.5), ncol = J, nrow = K)
  Theta <- apply(Theta, 2, function(x) x * exposures)
  rownames(Theta) <- colnames(Rmat)
  
  # Generate the counts
  Lambda <- Rmat %*% Theta
  X <- matrix(rnbinom(length(Lambda), size = 1/overdispersion, mu = c(Lambda)), 
              nrow = nrow(cosmic_data), ncol = J)
  rownames(X) <- rownames(Rmat) 
  return(list(X = X, Rmat = Rmat, Theta = Theta))
}

