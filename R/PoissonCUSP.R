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


sample_Z <- function(nu, sum_Y, sum_Eps, a, mu_inf, a0, b0){
  K <- length(nu)
  # Compute log w
  log_w <- log(nu) + c(0, cumsum(log(1 - nu[-K])))
  # Compute densities
  Z <- rep(0, K)
  for(k in 1:K){
    # Compute the probabilities
    log_prob_z <- rep(0, K)
    if(k < K){
      log_prob_z[1:k] <- log_w[1:k] + sum_Y[k] * log(mu_inf) - mu_inf * sum_Eps[k]
      log_prob_z[(k+1):K] <- log_w[(k+1):K] +  a0 * log(b0) - lgamma(a0)  +  lgamma(a0 + sum_Y[k]) -
        (a0 + sum_Y[k]) * log(b0 + sum_Eps[k])
    } else {
      log_prob_z[1:k] <- log_w[1:k] + sum_Y[k] * log(mu_inf) - mu_inf * sum_Eps[k] 
    }
    prob_z <- exp(log_prob_z - max(log_prob_z))
    Z[k] <- sample(1:K, 1, prob = prob_z)
  }
  return(Z)
}

sample_mu_sp <- function(Z, a0, b0, sum_Y, sum_Eps, mu_inf){
  K <- length(Z)
  mu <- rep(mu_inf, K)
  id_slab <- which(Z > 1:K)
  mu[id_slab] <- rgamma(length(id_slab), a0 + sum_Y[id_slab], b0 + sum_Eps[id_slab])
  return(mu)
}



# Main function to sample from the Poisson CUSP
PoissonCUSP <- function(X, K, nsamples = 2000, burnin = 1000,
                        alpha = 0.5, a = 1, a0 = 1, b0 = 1,
                        mu_inf = 0.01,
                        alpha_sp = 5, alpha0 = -1, alpha1 = -5e-04,
                        adapt_cutoff = 500) {
  t_start <- Sys.time()
  I <- nrow(X)
  J <- ncol(X)
  
  # Store the output
  SIGN <- vector(mode = "list", length = nsamples) 
  THETA <- vector(mode = "list", length = nsamples)
  MU <- vector(mode = "list", length = nsamples)
  ID_Dropped <- vector(mode = "list", length = nsamples)
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
        K <- length(mu)
        id_to_keep <- which(mu != mu_inf)
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
    Z <- sample_Z(nu, sum_Y, sum_Eps, a, mu_inf, a0, b0) 
    # Step 2 - Sample nu
    nu <- sample_nu(Z, alpha_sp)
    # Step 3 - Sample mu
    mu <- sample_mu_sp(Z, a0, b0, sum_Y, sum_Eps, mu_inf)
    # Update Theta
    Theta <- Eps * mu
    #------------------------------ 5. Store the output
    if(iter > burnin) {
      SIGN[[iter - burnin]] <- R
      THETA[[iter - burnin]] <- Theta
      MU[[iter - burnin]] <- mu
    }
  }
  time <- difftime(Sys.time(), t_start, units = "secs")[[1]]
  return(list(Signatures = SIGN, 
              Weights = THETA,
              Mu = MU, 
              spike = mu_inf,
              time = time))
}













