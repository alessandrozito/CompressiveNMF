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

