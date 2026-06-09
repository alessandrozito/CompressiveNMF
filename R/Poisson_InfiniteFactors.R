# This function implements the Poisson factorization using the CUSP algorithm
# in Legramanti et al. (2020). 

# Sample the mutational signature matrix
sample_signatures <- function(Alpha) {
  # Alpha is a matrix of I x K
  I <- nrow(Alpha)
  Ktot <- ncol(Alpha)
  R <- matrix(rgamma(I * Ktot, Alpha), nrow = I) + 1e-10 # Small nugget to avoid degeneracies of the gamma prior
  R <- apply(R, 2, function(x) x/sum(x))
  colnames(R) <- colnames(Alpha)
  return(R)
}

# Sample the weight matrix
sample_weights <- function(shape_mat, rate_mat) {
  # Alpha is a matrix of I x K
  K <- nrow(shape_mat)
  J <- ncol(shape_mat)
  Theta <- matrix(rgamma(K * J, shape_mat, rate_mat), nrow = K)
  rownames(Theta) <- rownames(shape_mat)
  return(Theta)
}

# Sample the augmented variables
sample_Y <- function(X, R, Theta, nonzero_ids) {
  I <- nrow(R)
  K <- ncol(R)
  J <- ncol(Theta)
  Y <- array(0, dim = c(I, J, K))
  for(id in 1:nrow(nonzero_ids)){
    i <- nonzero_ids[id, 1]
    j <- nonzero_ids[id, 2]
    q_vec <- R[i, ] * Theta[, j]
    Y[i, j, ] <- rmultinom(n = 1, size = X[i, j], prob = q_vec)
  }
  return(Y)
}


sample_nu <- function(Z, alpha_sp) {
  K <- length(Z)
  a_beta <- 1 + sapply(1:(K-1), function(l) sum(Z == l))
  b_beta <- alpha_sp + sapply(1:(K-1), function(l) sum(Z > l))
  nu <- c(rbeta(K-1, a_beta, b_beta), 1)
  return(nu)
}

sample_Z <- function(nu, sum_Y, sum_Eps, a, mu_inf, a0, b0, KowalCanale = FALSE){
  K <- length(nu)
  # Compute log w
  log_w <- log(nu) + c(0, cumsum(log(1 - nu[-K])))
  # Compute densities
  Z <- rep(0, K)
  for(k in 1:K){
    # Compute the probabilities
    log_prob_z <- rep(0, K)
    if(k < K){
      # Spike components
      if(KowalCanale) {
        log_prob_z[1:k] <- log_w[1:k] +  a0 * log(b0/mu_inf) - lgamma(a0)  +  lgamma(a0 + sum_Y[k]) -
          (a0 + sum_Y[k]) * log(b0/mu_inf + sum_Eps[k])
      } else {
        log_prob_z[1:k] <- log_w[1:k] + sum_Y[k] * log(mu_inf) - mu_inf * sum_Eps[k]
      }
      # Slab components
      log_prob_z[(k+1):K] <- log_w[(k+1):K] +  a0 * log(b0) - lgamma(a0)  +  lgamma(a0 + sum_Y[k]) -
        (a0 + sum_Y[k]) * log(b0 + sum_Eps[k])
    } else {
      if(KowalCanale) {
        log_prob_z[1:k] <- log_w[1:k] +  a0 * log(b0/mu_inf) - lgamma(a0)  +  lgamma(a0 + sum_Y[k]) -
          (a0 + sum_Y[k]) * log(b0/mu_inf + sum_Eps[k])
      } else {
        log_prob_z[1:k] <- log_w[1:k] + sum_Y[k] * log(mu_inf) - mu_inf * sum_Eps[k]
      }
    }
    prob_z <- exp(log_prob_z - max(log_prob_z))
    Z[k] <- sample(1:K, 1, prob = prob_z)
  }
  return(Z)
}

sample_Z_ig <- function(nu, J, sum_Theta, a, mu_inf, a0, b0, KowalCanale = FALSE){
  K <- length(nu)
  # Compute log w
  log_w <- log(nu) + c(0, cumsum(log(1 - nu[-K])))
  # Compute densities
  Z <- rep(0, K)
  for(k in 1:K){
    # Compute the probabilities
    log_prob_z <- rep(0, K)
    if(k < K){
      # Spike components
      if(KowalCanale) {
        log_prob_z[1:k] <- log_w[1:k] +  a0 * log(b0 * mu_inf) - lgamma(a0)  +  
          lgamma(a0 + a * J) - (a0 + a * J) * log(b0 * mu_inf + a * sum_Theta[k])
      } else {
        log_prob_z[1:k] <- log_w[1:k] - J * a * log(mu_inf) - a * sum_Theta[k] / mu_inf
      }
      # Slab components
      log_prob_z[(k+1):K] <- log_w[(k+1):K] +  a0 * log(b0) - lgamma(a0)  +  
        lgamma(a0 + a * J) - (a0 + a * J) * log(b0 + a * sum_Theta[k])
    } else {
      if(KowalCanale) {
        log_prob_z[1:k] <- log_w[1:k] +  a0 * log(b0 * mu_inf) - lgamma(a0)  +  
          lgamma(a0 + a * J) - (a0 + a * J) * log(b0 * mu_inf + a * sum_Theta[k])
      } else {
        log_prob_z[1:k] <- log_w[1:k] - J * a * log(mu_inf) - a * sum_Theta[k] / mu_inf
      }
    }
    prob_z <- exp(log_prob_z - max(log_prob_z))
    Z[k] <- sample(1:K, 1, prob = prob_z)
  }
  return(Z)
}

sample_mu_sp <- function(Z, a0, b0, sum_Y, sum_Eps, mu_inf, KowalCanale = FALSE){
  K <- length(Z)
  mu <- rep(mu_inf, K)
  id_slab <- which(Z > 1:K)
  mu[id_slab] <- rgamma(length(id_slab), a0 + sum_Y[id_slab], b0 + sum_Eps[id_slab])
  if (KowalCanale) {
    mu[-id_slab] <- rgamma(K - length(id_slab), a0 + sum_Y[-id_slab], 
                           b0/mu_inf + sum_Eps[-id_slab])
  }
  return(mu)
}

sample_mu_sp_ig <- function(Z, a, J, a0, b0, sum_Theta, mu_inf, KowalCanale = FALSE){
  K <- length(Z)
  mu <- rep(mu_inf, K)
  id_slab <- which(Z > 1:K)
  mu[id_slab] <- 1/rgamma(length(id_slab), a0 + a * J, b0 + a * sum_Theta[id_slab])
  if (KowalCanale) {
    mu[-id_slab] <- 1/rgamma(K - length(id_slab), a0 + a * J, mu_inf * b0 + a * sum_Theta[-id_slab])
  }
  return(mu)
}

sample_alpha_sp <- function(nu, c0, d0){
  K <- length(nu)
  alpha_sp <- rgamma(1, c0 + K - 1, d0 - sum(log(1 - nu[1:(K-1)])))
  return(alpha_sp)
}

# Main function to sample from the Poisson CUSP
PoissonCUSP <- function(X, K, nsamples = 2000, burnin = 1000,
                        alpha = 0.5, a = 1, a0 = 1, b0 = 1,
                        mu_inf = 0.01,
                        alpha_sp = 5, random_alpha_sp = FALSE, c0 = 2, d0 = 1,
                        alpha0 = -1, alpha1 = -5e-04,
                        adapt_cutoff = 500, KowalCanale = FALSE) {
  t_start <- Sys.time()
  I <- nrow(X)
  J <- ncol(X)
  
  # Store the output
  SIGN <- vector(mode = "list", length = nsamples) 
  THETA <- vector(mode = "list", length = nsamples)
  MU <- vector(mode = "list", length = nsamples)
  ALPHASP <- rep(NA, length = nsamples)
  Zall <- vector(mode = "list", length = nsamples)
  
  # Threshold for sparsity
  nu <- c(rbeta(K-1, 1, alpha_sp), 1)
  
  # Initialization of the sampler
  Alpha <- matrix(alpha, nrow = I, ncol = K)
  # Sample signatures
  R <- sample_signatures(Alpha)
  # Sample weights
  shape_mat <- matrix(a, nrow = K, ncol = J)
  rate_mat <- matrix(a, nrow = K, ncol = J)
  Eps <- sample_weights(shape_mat, rate_mat)
  # Sample the augmented variables
  nonzero_ids <- which(X != 0, arr.ind = TRUE)
  # Sample the global means
  mu <- rgamma(K, a0, b0)
  Theta <- Eps * mu
  verbose_step <- round((nsamples + burnin)/10)
  adapt <- 0
  for(iter in 1:(nsamples + burnin)) {
    if(iter %% verbose_step == 0){
      print(paste0("Iteration: ", iter, " [", round(iter/(nsamples + burnin) * 100), "%]"))
    } 
    
    #------------------------------ 0. Adaptive truncation move
    if(iter > adapt_cutoff){
      if(log(runif(1)) < alpha0 + alpha1 * iter){
        adapt <- 1
        # Check what columns need to be dropped
        K <- length(Z)#length(mu)
        id_to_keep <- which(Z > 1:K)#which(mu != mu_inf)
        K_star <- length(id_to_keep)
        if(K_star < K - 1) {
          # Drop the spike components 
          Theta <- Theta[id_to_keep, ]
          R <- R[, id_to_keep]
          mu <- mu[id_to_keep]
          nu <- nu[id_to_keep]
          # Add a final component sampled from the prior
          Theta <- rbind(Theta, rgamma(J, a, a) * mu_inf)
          r_new <- rgamma(I, alpha) 
          R <- cbind(R, r_new/sum(r_new))
          mu <- c(mu, mu_inf)
          nu <- c(nu, 1)
        } else {
          Theta <- rbind(Theta, rgamma(J, a, a) * mu_inf)
          r_new <- rgamma(I, alpha) 
          R <- cbind(R, r_new/sum(r_new))
          mu <- c(mu, mu_inf)
          nu[length(nu)] <- rbeta(1, 1, alpha)
          nu <- c(nu, 1)
        }
      }
    }
    
    #------------------------------ 1. Sample the latent variables from multinomial
    Y <- sample_Y(X = X, R = R, Theta = Theta, nonzero_ids = nonzero_ids)
    #------------------------------ 2. Sample the weights
    shape_mat <- a + apply(Y, c(3, 2), sum)
    rate_mat <-  a + matrix(mu) [, rep(1, J)]
    Eps <- sample_weights(shape_mat, rate_mat)
    #------------------------------ 3. Sample the signatures
    Alpha <- alpha + apply(Y, c(1, 3), sum)
    R <- sample_signatures(Alpha)
    #------------------------------ 4. Sample the global column mean via spike and slab
    sum_Y <- apply(Y, 3, sum)
    sum_Eps <- rowSums(Eps)
    # Step 1 - Sample Z
    Z <- sample_Z(nu, sum_Y, sum_Eps, a, mu_inf, a0, b0, KowalCanale) 
    # Step 2 - Sample nu
    nu <- sample_nu(Z, alpha_sp)
    # Step 3 - Sample mu
    mu <- sample_mu_sp(Z, a0, b0, sum_Y, sum_Eps, mu_inf, KowalCanale)
    # Update Theta
    Theta <- Eps * mu
    # Optional: update alpha_sp hyperparameter
    if(random_alpha_sp){
      alpha_sp <- sample_alpha_sp(nu, c0, d0)
    }
    #------------------------------ 5. Store the output
    if(iter > burnin) {
      SIGN[[iter - burnin]] <- R
      THETA[[iter - burnin]] <- Theta
      MU[[iter - burnin]] <- mu
      ALPHASP[[iter - burnin]] <- alpha_sp
      Zall[[iter - burnin]] <- Z
    }
  }
  time <- difftime(Sys.time(), t_start, units = "secs")[[1]]
  return(list(Signatures = SIGN, 
              Weights = THETA,
              Mu = MU, 
              alpha_sp = ALPHASP,
              Z = Zall,
              spike = mu_inf,
              time = time))
}


# Main function to sample from the Poisson CUSP
PoissonCUSP_ig <- function(X, K, nsamples = 2000, burnin = 1000,
                           alpha = 0.5, a = 1, a0 = 2.5, b0 = 1.5,
                           mu_inf = 0.01,
                           alpha_sp = 5, random_alpha_sp = FALSE, c0 = 2, d0 = 1,
                           alpha0 = -1, alpha1 = -5e-04,
                           adapt_cutoff = 500, KowalCanale = FALSE) {
  t_start <- Sys.time()
  I <- nrow(X)
  J <- ncol(X)
  
  # Store the output
  SIGN <- vector(mode = "list", length = nsamples) 
  THETA <- vector(mode = "list", length = nsamples)
  MU <- vector(mode = "list", length = nsamples)
  ALPHASP <- rep(NA, length = nsamples)
  Zall <- vector(mode = "list", length = nsamples)
  
  # Threshold for sparsity
  nu <- c(rbeta(K-1, 1, alpha_sp), 1)
  
  # Initialization of the sampler
  Alpha <- matrix(alpha, nrow = I, ncol = K)
  # Sample signatures
  R <- sample_signatures(Alpha)
  # Sample the augmented variables
  nonzero_ids <- which(X != 0, arr.ind = TRUE)
  # Sample the global means
  mu <- rgamma(K, a0, b0)
  # Sample weights
  shape_mat <- matrix(a, nrow = K, ncol = J)
  rate_mat <- matrix(a, nrow = K, ncol = J) / mu
  Theta <- sample_weights(shape_mat, rate_mat)
  # Run the sampler
  verbose_step <- round((nsamples + burnin)/10)
  adapt <- 0
  for(iter in 1:(nsamples + burnin)) {
    if(iter %% verbose_step == 0){
      print(paste0("Iteration: ", iter, " [", round(iter/(nsamples + burnin) * 100), "%]"))
    } 
    
    #------------------------------ 0. Adaptive truncation move
    if(iter > adapt_cutoff){
      if(log(runif(1)) < alpha0 + alpha1 * iter){
        adapt <- 1
        # Check what columns need to be dropped
        K <- length(Z)#length(mu)
        id_to_keep <- which(Z > 1:K)#which(mu != mu_inf)
        K_star <- length(id_to_keep)
        if(K_star < K - 1) {
          # Drop the spike components 
          Theta <- Theta[id_to_keep, ]
          R <- R[, id_to_keep]
          mu <- mu[id_to_keep]
          nu <- nu[id_to_keep]
          # Add a final component sampled from the prior
          Theta <- rbind(Theta, rgamma(J, a, a) * mu_inf)
          r_new <- rgamma(I, alpha) 
          R <- cbind(R, r_new/sum(r_new))
          mu <- c(mu, mu_inf)
          nu <- c(nu, 1)
        } else {
          Theta <- rbind(Theta, rgamma(J, a, a) * mu_inf)
          r_new <- rgamma(I, alpha) 
          R <- cbind(R, r_new/sum(r_new))
          mu <- c(mu, mu_inf)
          nu[length(nu)] <- rbeta(1, 1, alpha_sp)
          nu <- c(nu, 1)
        }
      }
    }
    
    #------------------------------ 1. Sample the latent variables from multinomial
    Y <- sample_Y(X = X, R = R, Theta = Theta, nonzero_ids = nonzero_ids)
    #------------------------------ 2. Sample the weights
    shape_mat <- a + apply(Y, c(3, 2), sum)
    rate_mat <-  1 + a / matrix(mu) [, rep(1, J)]
    Theta <- sample_weights(shape_mat, rate_mat)
    #------------------------------ 3. Sample the signatures
    Alpha <- alpha + apply(Y, c(1, 3), sum)
    R <- sample_signatures(Alpha)
    #------------------------------ 4. Sample the global column mean via spike and slab
    sum_Y <- apply(Y, 3, sum)
    sum_Theta <- rowSums(Theta)
    # Step 1 - Sample Z
    Z <- sample_Z_ig(nu = nu, J = J, sum_Theta = sum_Theta, a =  a, 
                     mu_inf = mu_inf, a0 = a0, b0 = b0, KowalCanale = KowalCanale) 
    # Step 2 - Sample nu
    nu <- sample_nu(Z, alpha_sp)
    # Step 3 - Sample mu
    mu <- sample_mu_sp_ig(Z = Z, a = a, J = J, a0 = a0, b0 = b0, 
                          sum_Theta = sum_Theta, mu_inf = mu_inf, 
                          KowalCanale = KowalCanale)
    # Optional: update alpha_sp hyperparameter
    if(random_alpha_sp){
      alpha_sp <- sample_alpha_sp(nu, c0, d0)
    }
    #------------------------------ 5. Store the output
    if(iter > burnin) {
      SIGN[[iter - burnin]] <- R
      THETA[[iter - burnin]] <- Theta
      MU[[iter - burnin]] <- mu
      ALPHASP[[iter - burnin]] <- alpha_sp
      Zall[[iter - burnin]] <- Z
    }
  }
  time <- difftime(Sys.time(), t_start, units = "secs")[[1]]
  return(list(Signatures = SIGN, 
              Weights = THETA,
              Mu = MU, 
              alpha_sp = ALPHASP,
              Z = Zall,
              spike = mu_inf,
              time = time))
}

# Main function to sample from the Poisson Multiplicative Gamma Process (MGP)
PoissonMGP <- function(X, K, nsamples = 2000, burnin = 1000,
                       alpha = 0.5, a = 1, c0 = 2.1, d0 = 3.1) {
  t_start <- Sys.time()
  I <- nrow(X)
  J <- ncol(X)
  
  # Store the output
  # SIGN <- vector(mode = "list", length = nsamples) 
  # THETA <- vector(mode = "list", length = nsamples)
  # MU <- vector(mode = "list", length = nsamples)
  # NU <- vector(mode = "list", length = nsamples)
  SIGN <- array(NA, dim = c(nsamples, nrow(X), K))
  THETA <- array(NA, dim = c(nsamples, K, ncol(X)))
  MU <- matrix(NA, nrow = nsamples, ncol = K)
  NU <- matrix(NA, nrow = nsamples, ncol = K)

  # Initialization of the sampler
  Alpha <- matrix(alpha, nrow = I, ncol = K)
  # Sample signatures
  R <- sample_signatures(Alpha)
  # Sample weights
  shape_mat <- matrix(a, nrow = K, ncol = J)
  rate_mat <- matrix(a, nrow = K, ncol = J)
  Eps <- sample_weights(shape_mat, rate_mat)
  # Sample the augmented variables
  nonzero_ids <- which(X != 0, arr.ind = TRUE)
  # Sample the global means from the inverse gamma product
  nu <- 1/c(rgamma(1, c0, 1), rgamma(K - 1,  d0, 1))
  mu <- cumprod(nu)
  Theta <- Eps
  verbose_step <- round((nsamples + burnin)/10)
  adapt <- 0
  for(iter in 1:(nsamples + burnin)) {
    if(iter %% verbose_step == 0){
      print(paste0("Iteration: ", iter, " [", round(iter/(nsamples + burnin) * 100), "%]"))
    } 
    
    #------------------------------ 1. Sample the latent variables from multinomial
    Y <- sample_Y(X = X, R = R, Theta = Theta, nonzero_ids = nonzero_ids)
    #------------------------------ 2. Sample the weights
    shape_mat <- a + apply(Y, c(3, 2), sum)
    rate_mat <-  1 + a/matrix(mu) [, rep(1, J)]
    Theta <- sample_weights(shape_mat, rate_mat)
    #------------------------------ 3. Sample the signatures
    Alpha <- alpha + apply(Y, c(1, 3), sum)
    R <- sample_signatures(Alpha)
    #------------------------------ 4. Sample the global means
    ThetaSums <- rowSums(Theta)
    # Update the first entry
    nu[1] <- 1/rgamma(1, c0 + J * K, 1 + sum(a *  nu[1] * ThetaSums/cumprod(nu)))
    # Update all the others
    for (l in 2:K){
      pars2 <- sum(a * nu[l] * ThetaSums[l:K] / cumprod(nu)[l:K])
      nu[l] <- 1/rgamma(1, d0 + J * (K - l + 1), 1 + pars2)
    }
    mu <- cumprod(nu)
    
    #------------------------------ 5. Store the output
    if(iter > burnin) {
      SIGN[iter - burnin, ,] <- R
      THETA[iter - burnin, ,] <- Theta
      NU[iter - burnin, ] <- nu
      MU[iter - burnin, ] <- mu
    }
  }
  time <- difftime(Sys.time(), t_start, units = "secs")[[1]]
  return(list(Signatures = SIGN, 
              Weights = THETA,
              Mu = MU, 
              Nu = NU, 
              time = time))
}





# get_posterior_CUSP_v2 <- function(resCUSP) {
#   J <- ncol(resCUSP$Weights[[1]])
#   Kmax <- max(unlist(lapply(resCUSP$Mu, length)))
#   nspike <- rep(0, Kmax)
#   nsims <- length(resCUSP$Signatures)
#   
#   k0 <- length(resCUSP$Mu[[1]])
#   R_hat <- cbind(resCUSP$Signatures[[1]] , matrix(1 / 96, nrow = 96, ncol = Kmax - k0))
#   Theta_hat <- rbind(resCUSP$Weights[[1]], matrix(0, ncol = J, nrow =  Kmax - k0))
#   Mu_hat <- c(resCUSP$Mu[[1]], rep(resCUSP$spike, Kmax - k0))
#   mu_seq <- t(c(resCUSP$Mu[[1]], rep(resCUSP$spike, Kmax - k0)))
#   
#   for (i in 2:nsims) {
#     mu_temp <- resCUSP$Mu[[i]]; mu_old <- resCUSP$Mu[[i-1]]; Zold <- resCUSP$Z[[i-1]]
#     # Check if a column was dropped
#     k_new <- length(mu_temp); k_old <- length(mu_old)
#     if(k_new < k_old){
#       # Check which column was dropped
#       id_dropped <- which(Zold <= 1:k_old)#which(mu_old == resCUSP$spike) # These are the columns that were dropped. We need to take the resulting mean and move it to the end
#       R_hat <- cbind(R_hat[, -id_dropped], R_hat[, id_dropped])
#       Theta_hat <- rbind(Theta_hat[-id_dropped, ], Theta_hat[id_dropped, ])
#       Mu_hat <- c(Mu_hat[-id_dropped], Mu_hat[id_dropped])
#       if(i == 2){
#         mu_seq <- cbind(t(mu_seq[, -id_dropped]), t(mu_seq[, id_dropped]))  
#       } else {
#         mu_seq <- cbind(mu_seq[, -id_dropped], mu_seq[, id_dropped])
#       }
#       nspike <- c(nspike[-id_dropped], nspike[id_dropped])
#     }
#     # Check global mean
#     mu <- c(mu_temp, rep(resCUSP$spike, Kmax - length(mu_temp)))
#     mu_seq <- rbind(mu_seq, mu)
#     Mu_hat <- Mu_hat + mu
#     nspike <- nspike + 1 * (mu == resCUSP$spike)
#     # Check signatures
#     mat_temp <- resCUSP$Signatures[[i]]
#     R_hat <- R_hat + cbind(mat_temp, matrix(1 / 96, nrow = 96, ncol = Kmax - length(mu_temp)))
#     # Check weights
#     Theta_temp <- resCUSP$Weights[[i]]
#     Theta_hat <- Theta_hat + rbind(Theta_temp, matrix(0, ncol = ncol(Theta_temp),
#                                                       nrow =  Kmax - length(mu_temp)))
#     #mu_seq <- rbind(mu_seq, mu)
#   }
#   R_hat <- R_hat / nsims
#   nspike <- nspike / nsims
#   Theta_hat <- Theta_hat / nsims
#   Mu_hat <- Mu_hat / nsims
#   return(list("Theta_hat" = Theta_hat, "R_hat" = R_hat,
#               "Mu_hat" = Mu_hat, "nspike" = nspike, "mu_seq" = mu_seq))
# }
# 
# 
# Postprocess_PoissonCUSP_v2 <- function(resCUSP, data) {
#   # Step 1 - find the number of signatures
#   Kchain <- unlist(lapply(resCUSP$Z, function(x) sum(x > 1:length(x))))
#   K <- mean(Kchain)
#   # Step 2 - Calculate RMSE wrt to the count matrix and the rate matrix
#   Lambda <- get_Lambda_CUSP(resCUSP)
#   Lambda_true <- data$Rmat %*% data$Theta
#   rmse_Lambda <- sqrt(mean((Lambda_true - Lambda)^2))
#   rmse_Counts <- sqrt(mean((data$X - Lambda)^2))
#   # Step 3 - Infer the mutational signatures
#   R_true <- apply(data$Rmat, 2, function(x) x / sum(x))
#   post_CUSP <- get_posterior_CUSP_v2(resCUSP)
#   nspike <- post_CUSP$nspike
#   R_hat <- post_CUSP$R_hat[, nspike < 0.05]
#   matchedSign <- match_MutSign(R_true = R_true, R_hat = R_hat)
#   cos_sim <- mean(get_cosine_similarity(matchedSign))
#   # Step 4 - Calculate RMSE for signatures and Weights
#   Theta_hat <- post_CUSP$Theta_hat[nspike < 0.05, ]
#   rmse_R <- compute_RMSE_Signature(R_hat = matchedSign$R_hat, R_true = matchedSign$R_true)  
#   rmse_Theta <- compute_RMSE_Theta(Theta_true = data$Theta, Theta_hat = Theta_hat, matchedSign$match)
#   # Step 5 - calculate the sensitivity and precision
#   sens_prec  <- Compute_sensitivity_precision(R_hat = R_hat, data$Rmat)
#   # Step 6 - add Effective sample sizes
#   effsize <- get_ESS_PoissonCUSP(resCUSP)
#   return(list("Lambda" = Lambda,
#               "R_hat" = R_hat, #post_CUSP$R_hat,
#               "Theta_hat" = Theta_hat, #post_CUSP$Theta_hat,
#               "Mu_hat" = post_CUSP$Mu_hat,
#               "Kchain" = Kchain,
#               "Mu_chain"= post_CUSP$mu_seq,
#               "nspike" = nspike,
#               "signatures" = matchedSign,
#               "results" = c("K" = K, 
#                             "rmse_Lambda" = rmse_Lambda, 
#                             "rmse_Counts" = rmse_Counts, 
#                             rmse_R, 
#                             rmse_Theta, 
#                             sens_prec,
#                             "cos_sim" = cos_sim, 
#                             "time" = resCUSP$time, 
#                             effsize)
#   ))
# }
