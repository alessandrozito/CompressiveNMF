
generate_data <- function(I = 96, J, K0) {
  # I = 96 mutational signatures
  # J = number of patients (or units)
  # K0 = number of latent signatures
  #----------------------------------- Signature matrix
  Rmat <- matrix(rgamma(n = I * K0, 0.5, 0.5), nrow = I, ncol = K0)
  #Rmat <- apply(Rmat, 2, function(x) x /sum(x))
  #----------------------------------- Weights
  Theta <- matrix(rgamma(n = K0 * J, 1, 0.5), nrow = K0, ncol = J)
  #----------------------------------- Generate the counts
  Lambda <- Rmat %*% Theta
  X <- matrix(rpois(n = length(Lambda), c(Lambda)), nrow = I, ncol = J)
  return(list(X = X, Rmat = Rmat, Theta = Theta))
}

generate_data2 <- function(I = 96, J, K0) {
  # I = 96 mutational signatures
  # J = number of patients (or units)
  # K0 = number of latent signatures
  #----------------------------------- Signature matrix
  Rmat <- matrix(rgamma(n = I * K0, 10/I, 1), nrow = I, ncol = K0)
  Rmat <- apply(Rmat, 2, function(x) x /sum(x))
  #----------------------------------- Weights
  Theta <- matrix(rgamma(n = K0 * J, 10, 10/200), nrow = K0, ncol = J)
  #----------------------------------- Generate the counts
  Lambda <- Rmat %*% Theta
  X <- matrix(rpois(n = length(Lambda), c(Lambda)), nrow = I, ncol = J)
  return(list(X = X, Rmat = Rmat, Theta = Theta))
}



# Sample the mutational signature matrix
sample_signatures <- function(Alpha) {
  # Alpha is a matrix of I x K
  I <- nrow(Alpha)
  K <- ncol(Alpha)
  R <- matrix(rgamma(I * K, Alpha), nrow = I)
  return(apply(R, 2, function(x) x/sum(x)))
}

# Sample the weight matrix
sample_weights <- function(shape_mat, rate_mat) {
  # Alpha is a matrix of I x K
  K <- nrow(shape_mat)
  J <- ncol(shape_mat)
  Theta <- matrix(rgamma(K * J, shape_mat, rate_mat), nrow = K)
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


compute_logw <- function(nu){
  log_w <- log(nu) + c(0, cumsum(log(1 - nu[-K])))
  return(log_w)
}

sample_nu <- function(Z, alpha_sp) {
  K <- length(Z)
  a_beta <- 1 + sapply(1:(K-1), function(l) sum(Z == l))
  b_beta <- alpha_sp + sapply(1:(K-1), function(l) sum(Z > l))
  nu <- c(rbeta(K-1, a_beta, b_beta), 1)
  return(nu)
}


sample_Z <- function(nu, Theta, a, mu_inf, a0, b0){
  K <- nrow(Theta)
  # Compute log w
  log_w <- log(nu) + c(0, cumsum(log(1 - nu[-K])))
  # Compute densities
  Z <- rep(0, K)
  for(k in 1:K){
    # Compute the probabilities
    log_prob_z <- rep(0, K)
    if(k < K){
      log_prob_z[1:k] <- log_w[1:k] + sum(dgamma(Theta[k,], a, a/mu_inf, log = TRUE))  
      log_prob_z[(k+1):K] <- log_w[(k+1):K] + J * a * log(a) + lgamma(a0 + J * a) - J * lgamma(a) - lgamma(a0)  + 
        a0 * log(b0) - (a0 + J * a) * log(a * sum(Theta[k, ]) + b0) + (a - 1) * sum(log(Theta[k,])) 
    } else {
      log_prob_z[1:k] <- log_w[1:k] + sum(dgamma(Theta[k,], a, a/mu_inf, log = TRUE))  
    }
    prob_z <- exp(log_prob_z - max(log_prob_z))
    Z[k] <- sample(1:K, 1, prob = prob_z)
  }
  return(Z)
}

sample_mu_sp <- function(Z, a0, b0, Theta, mu_inf){
  K <- length(Z)
  mu <- rep(mu_inf, length(Z))
  id_slab <- which(Z > 1:K)
  mu[id_slab] <- 1/rgamma(length(id_slab), rep(a0 + J * a, K), b0 + a * rowSums(Theta)[id_slab])
  return(mu)
}

# This function runs the Gibbs sampler for Bayesian NMF 
BayesNMF <- function(X, K, nsamples = 2000, burnin = 1000,
                     alpha = 1, a = 1, b = 1, shrink = "none",
                     a0 = 2, b0 = 1, mu_inf = 0.01, alpha_sp = 5) {
  
  if(!shrink %in% c("none", "ARD", "CUSP")){
    stop(paste0("Method shrink = ", shrink, " not supported. Options are 'none', 'ARD' and 'CUSP'"))
  }
  
  I <- nrow(X)
  J <- ncol(X)
  
  # Store the output
  SIGN <- array(dim = c(nsamples, I, K))
  THETA <- array(dim = c(nsamples, K, J))
  Y_means <- matrix(nrow = nsamples, ncol = K)
  if(shrink != "none") {
    MU <- matrix(nrow = nsamples, ncol = K)
  } else {
    MU <- NULL
  }
  
  if(shrink == "CUSP") {
    nu <- c(rbeta(K-1, 1, alpha_sp), 1)
  }
  
  # Initialization of the sampler
  Alpha <- matrix(alpha, nrow = I, ncol = K)
  # Sample signatures
  R <- sample_signatures(Alpha)
  # Sample weights
  shape_mat <- matrix(a, nrow = K, ncol = J)
  rate_mat <- matrix(a/mu_inf, nrow = K, ncol = J)
  Theta <- sample_weights(shape_mat, rate_mat)
  # Sample the augmented variables
  nonzero_ids <- which(X != 0, arr.ind = TRUE)
  
  if(shrink != "none") {
    mu <- rep(0.01, K)#1/rgamma(K, a0, b0)
  }
  
  verbose_step <- round((nsamples + burnin)/10)
  for(iter in 1:(nsamples + burnin)) {
    if(iter %% verbose_step == 0){
      print(paste0("Iteration: ", iter, " [", round(iter/(nsamples + burnin) * 100), "%]"))
    } 
    
    #------------------------------ 1. Sample the latent variables from multinomial
    Y <- sample_Y(X = X, R = R, Theta = Theta, nonzero_ids = nonzero_ids)
    #------------------------------ 2. Sample the weights
    shape_mat <- a + apply(Y, c(3, 2), sum)
    if(shrink != "none"){
      rate_mat <- 1 + matrix(a / mu) [, rep(1, J)]
    } else {
      rate_mat <- matrix(b + 1, nrow = K,  ncol = J)
    }
    Theta <- sample_weights(shape_mat, rate_mat)
    #------------------------------ 3. Sample the signatures
    Alpha <- alpha + apply(Y, c(1, 3), sum)
    R <- sample_signatures(Alpha)
    #------------------------------ 4. Optional step: sample the global column mean
    if(shrink == "ARD") {
      mu <- 1/rgamma(K, rep(a0 + J * a, K), b0 + a * rowSums(Theta))
    } else if (shrink == "CUSP") {
      # Step 1 - Sample Z
      Z <- sample_Z(nu, Theta, a, mu_inf, a0, b0)
      # Step 2 - Sample nu
      nu <- sample_nu(Z, alpha_sp)
      # Step 3 - Sample mu
      mu <- sample_mu_sp(Z, a0, b0, Theta, mu_inf, K)
    }
    #------------------------------ 5. Store the output
    if(iter > burnin) {
      SIGN[iter - burnin, ,] <- R
      THETA[iter - burnin, ,] <- Theta
      Y_means[iter - burnin, ] <- colMeans(apply(Y, c(2,3), sum))
      if(shrink != "none"){
        MU[iter - burnin, ] <- mu
      }
    }
  }
  return(list(Signatures = SIGN, 
              Weights = THETA,
              Mu = MU, 
              Y_means = Y_means))
}


# This function runs the Gibbs sampler for Bayesian NMF 
BayesNMF_adapt <- function(X, K, nsamples = 2000, burnin = 1000,
                     alpha = 1, a = 1, a0 = 1, b0 = 1, mu_inf = 0.01,
                     alpha_sp = 5, alpha0 = 1, alpha1 = -5e-4, adapt_cutoff = 500) {
  I <- nrow(X)
  J <- ncol(X)
  
  # Store the output
  SIGN <- vector(mode = "list", length = nsamples) #array(dim = c(nsamples, I, K))
  THETA <- vector(mode = "list", length = nsamples)#array(dim = c(nsamples, K, J))
  MU <- vector(mode = "list", length = nsamples)
  
  nu <- c(rbeta(K-1, 1, alpha_sp), 1)
  # Initialization of the sampler
  Alpha <- matrix(alpha, nrow = I, ncol = K)
  # Sample signatures
  R <- sample_signatures(Alpha)
  # Sample weights
  shape_mat <- matrix(a, nrow = K, ncol = J)
  rate_mat <- matrix(a/mu_inf, nrow = K, ncol = J)
  Theta <- sample_weights(shape_mat, rate_mat)
  # Sample the augmented variables
  nonzero_ids <- which(X != 0, arr.ind = TRUE)
  mu <- rep(0.01, K)
  verbose_step <- round((nsamples + burnin)/10)
  
  for(iter in 1:(nsamples + burnin)) {
    if(iter %% verbose_step == 0){
      print(paste0("Iteration: ", iter, " [", round(iter/(nsamples + burnin) * 100), "%]"))
    } 
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
          Theta <- rbind(Theta, rgamma(J, a, a/mu_inf))
          r_new <- rgamma(I, alpha) 
          R <- cbind(R, r_new/sum(r_new))
          mu <- c(mu, mu_inf)
          nu <- c(nu, 1)
        } else {
          Theta <- rbind(Theta, rgamma(J, a, a/mu_inf))
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
    rate_mat <- 1 + matrix(a / mu) [, rep(1, J)]
    Theta <- sample_weights(shape_mat, rate_mat)
    #------------------------------ 3. Sample the signatures
    Alpha <- alpha + apply(Y, c(1, 3), sum)
    R <- sample_signatures(Alpha)
    #------------------------------ 4. Optional step: sample the global column mean
    # Step 1 - Sample Z
    Z <- sample_Z(nu, Theta, a, mu_inf, a0, b0)
    # Step 2 - Sample nu
    nu <- sample_nu(Z, alpha_sp)
    # Step 3 - Sample mu
    mu <- sample_mu_sp(Z, a0, b0, Theta, mu_inf)
    #------------------------------ 5. Store the output
    if(iter > burnin) {
      SIGN[[iter - burnin]] <- R
      THETA[[iter - burnin]] <- Theta
      MU[[iter - burnin]] <- mu
    }
  }
  return(list(Signatures = SIGN, 
              Weights = THETA,
              Mu = MU))
}

### Other useful functions

plot_Lambda <- function(data, out, title = NULL){
  Theta <- apply(out$Weights, c(2,3), mean)
  R <- apply(out$Signatures, c(2,3), mean)
  
  Lambda <- R %*% Theta
  
  Lambda_true <- data$Rmat %*% data$Theta
  
  plot(Lambda, Lambda_true, main = title, asp = 1)
  abline(b = 1, col = "red", a = 0)
}


#
I <- 96
J <- 100
K0 <- 5
K <- 15
a <- 1
a0 <- a * J + 1
b0 <- 0.01 * (a0 - 1)
alpha_sp <- 5
alpha <- 1/I
nsamples <- 2000
burnin <- 1000
set.seed(10)
data <- generate_data(I = I, J = J, K0 = K0)

out <- BayesNMF(X = data$X, K = K, nsamples = nsamples,
                burnin = burnin, alpha = alpha, a = a,
                a0 = a0, b0 = b0, shrink = "ARD")

round(colMeans(out$Mu), 4)
plot_Lambda(data, out)



out2 <- BayesNMF(X = data$X, K = K, nsamples = nsamples,
                 burnin = burnin, alpha = 1, a = 1,
                 a0 = a0, b0 = b0, shrink = "ARD")
sort(round(colMeans(out2$Mu), 4))
plot_Lambda(data, out2)


plot(colMeans(out$Signatures[,,1]), data$Rmat[,5])
plot(colMeans(out$Signatures[,,1]))
plot(colMeans(out$Signatures[,,2]))
plot(colMeans(out$Signatures[,,5]))


# Match signatures via greedy algorithm
R_hat <- apply(out$Signatures, c(2,3), mean)
R_hat <- R_hat[,colMeans(out$Mu) > 0.04]
R_true <- apply(data$Rmat, 2, function(x) x /sum(x))

round(R_hat[,1], 5)

match_MutSign <- function(R_true, R_hat){
  # Need to make sure that R_hat and R_true are the same
  k_true <- ncol(R_true)
  k_hat <- ncol(R_hat)
  if(k_hat > k_true){
    R_true <- cbind(R_true, 1/nrow(R_true))  
  } else if (k_hat < k_true){
    R_hat <- cbind(R_hat, 1/nrow(R_hat))  
  }
  k <- ncol(R_hat)
  R_new <- matrix(NA, ncol = ncol(R_hat), nrow = nrow(R_hat))
  cols <- 1:ncol(R_new)
  for(k in 1:ncol(R_new)){
    r <- R_true[, k]
    dist <- apply(R_hat, 2, function(x) sum(abs(x - r)))
    i <- which.min(dist[cols])
    R_new[, k] <- R_hat[, cols][, i]
    cols <- cols[cols!=i]
  }
  return(R_new)
}

R_new <- match_MutSign(R_true, R_hat)
j <-3
plot(R_new[, j], R_true[, j])



# 
# out2 <- BayesNMF(X = data$X, K = K, nsamples = nsamples,
#                 burnin = burnin, alpha = 1, a = a,
#                 a0 = a0, b0 = b0, shrink = "ARD")
# plot_Lambda(data, out2)
# round(sort(colMeans(out2$Mu)), 4)
# 
# out_old_adp <- BayesNMF_adapt(X = data$X, K = K, nsamples = nsamples,
#                               burnin = burnin, alpha = 0.5, a = 1,
#                               a0 = a0, b0 = b0, mu_inf = 0.01, alpha_sp = alpha_sp)
# 
# table(unlist(lapply(out_old_adp$Mu, function(x) sum(x != 0.01))))
# 
# 
# plot(unlist(lapply(out_old_adp$Mu, function(mu) mu[2])), type = "l")
# 
# out_oldnospike <- BayesNMF_adapt(X = data$X, K = K, nsamples = nsamples,
#                                  burnin = burnin, alpha = 0.5, a = 1,
#                                  a0 = 1, b0 = 1, mu_inf = 0.01, alpha_sp = alpha_sp)
# 
# table(unlist(lapply(out_oldnospike$Mu, function(x) sum(x != 0.01))))
# 
# 
# 


fit <- nmf(data2$X, rank = 20)
dim(fit@fit@H)
dim(fit@fit@W)
fit@fit@H
rowMeans(fit@fit@H)








