# This function runs the gibbs sampler for the compressive non-negative matrix
# factorization algorithm. We consider the case when also known signatures are 
# included in the analysis, with the possibility to update them.

source("R/rInvKummer.R")

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


# Function to run our method
CompressiveNMF <- function(X, K, nsamples = 2000, burnin = 1000,
                           alpha = 0.5, a = 1, a0 = NULL, b0 = NULL, S = NULL,
                           betah = 100, useInvKummer = FALSE) {
  # Find the data
  I <- nrow(X)
  J <- ncol(X)
  if(!is.null(S)){
    H <- ncol(S)
    K <- K + H
  }
  # Store the output
  SIGN <- array(dim = c(nsamples, I, K))
  THETA <- array(dim = c(nsamples, K, J))
  MU <- matrix(nrow = nsamples, ncol = K)
  
  # Set the hyperparameters for the Compressive hyperprior
  if(is.null(a0)){
    a0 <- a * J + 1
    b0 <- a # implies that the mean for mu is 1/J, i.e. decreases with the sample size
  }
  
  # Initialization of the sampler
  Alpha <- matrix(alpha, nrow = I, ncol = K)
  if(!is.null(S)){
    Alpha[, 1:H] <- betah * S
  }
  # Sample signatures
  R <- sample_signatures(Alpha)
  # Sample weights
  shape_mat <- matrix(a, nrow = K, ncol = J)
  rate_mat <- matrix(a/0.01, nrow = K, ncol = J)
  Theta <- sample_weights(shape_mat, rate_mat)
  # Sample the augmented variables
  nonzero_ids <- which(X != 0, arr.ind = TRUE)
  mu <- rep(0.01, K)
  verbose_step <- round((nsamples + burnin)/10)
  for(iter in 1:(nsamples + burnin)) {
    if(iter %% verbose_step == 0){
      print(paste0("Iteration: ", iter, " [", round(iter/(nsamples + burnin) * 100), "%]"))
    } 
    #------------------------------ 1. Sample the latent variables from multinomial
    Y <- sample_Y(X = X, R = R, Theta = Theta, nonzero_ids = nonzero_ids)
    #------------------------------ 2. Sample the weights
    shape_mat <- a + apply(Y, c(3, 2), sum)
    rate_mat <- 1 + matrix(a / mu) [, rep(1, J)]
    Theta <- sample_weights(shape_mat, rate_mat)
    #------------------------------ 3. Sample the signatures
    Alpha <- alpha + apply(Y, c(1, 3), sum)
    if(!is.null(S)){
      Alpha[, 1:H] <- betah * S +  Alpha[, 1:H] - alpha
    }
    R <- sample_signatures(Alpha)
    if(!is.null(S)){
      colnames(R[,1:ncol(S)]) <- colnames(S)
    }
    #------------------------------ 4. Sample the global column mean
    if(useInvKummer){
      gamma_kum <- apply(Y, 3, sum) +  a * J
      mu <- sapply(1:K, function(i) rInvKummer(nsamples = 1,
                                               alpha = a0 + a * J,
                                               beta = b0,
                                               gamma = gamma_kum[i],
                                               delta = a))
    } else {
      mu <- 1/rgamma(K, rep(a0 + J * a, K), b0 + a * rowSums(Theta))  
    }
    
    #------------------------------ 5. Store the output
    if(iter > burnin) {
      SIGN[iter - burnin, ,] <- R
      THETA[iter - burnin, ,] <- Theta
      MU[iter - burnin, ] <- mu
    }
  }
  return(list(Signatures = SIGN, 
              Weights = THETA,
              Mu = MU))
}


# 
# ## Try the Cpp implementation
# library(RcppArmadillo)
# library(Rcpp)
# sourceCpp("src/BayesNFM_cpp.cpp")
# 
# CompressiveNMF_fast <- function(X, K, nsamples = 2000, burnin = 1000,
#                            alpha = 0.5, a = 1, a0 = NULL, b0 = NULL) {
#   # Find the data
#   I <- nrow(X)
#   J <- ncol(X)
#   
#   # Store the output
#   SIGN <- array(dim = c(nsamples, I, K))
#   THETA <- array(dim = c(nsamples, K, J))
#   MU <- matrix(nrow = nsamples, ncol = K)
#   
#   # Set the hyperparameters for the Compressive hyperprior
#   if(is.null(a0)){
#     a0 <- a * J + 1
#     b0 <- a # implies that the mean for mu is 1/J, i.e. decreases with the sample size
#   }
#   
#   # Initialization of the sampler
#   Alpha <- matrix(alpha, nrow = I, ncol = K)
#   # Sample signatures
#   R <- sample_signatures(Alpha)
#   # Sample weights
#   shape_mat <- matrix(a, nrow = K, ncol = J)
#   rate_mat <- matrix(a/0.01, nrow = K, ncol = J)
#   Theta <- sample_weights(shape_mat, rate_mat)
#   # Sample the augmented variables
#   nonzero_ids <- which(X != 0, arr.ind = TRUE)
#   mu <- rep(0.01, K)
#   b <- 1
#   # Call the cpp function
#   cpp_out <- CompressiveNMF_cpp(X, nonzero_ids, R, Theta, Y, mu, shape_mat, rate_mat, Alpha,
#                                 nsamples, burnin, alpha, a, b, a0, b0)
#  
#   return(cpp_out)
# }
# 
# 
# 
# # nsamples <- 2000
# # burnin <- 1000
# # 
# # out <- CompressiveNMF(X, K = 20, a0 = J + 1, b0 = 1, nsamples = nsamples, burnin = burnin)
# # sort(round(colMeans(out$Mu), 5))
# # 
# # outKu <- CompressiveNMF(X, K = 20, a0 = J + 1, b0 = 1, nsamples = nsamples, burnin = burnin, useInvKummer = TRUE)
# # 
# # sort(round(colMeans(out$Mu), 5))
# # sort(round(colMeans(outKu$Mu), 5))
# # 
# # # Greedily order and sort
# # postprocess_relevance_weights <- function(Mu, tau = 0.04){
# #   mu_hat <- colMeans(Mu)
# #   Mu <- Mu[, mu_hat >= tau]
# #   return(Mu[, order(mu_hat[mu_hat >= tau])])
# # }
# # 
# # Mu <- postprocess_relevance_weights(out$Mu)
# # MuKu <- postprocess_relevance_weights(outKu$Mu)
# # 
# # 
# # plot_autocorr <- function(x1, x2){
# #   plot(x1, type = "l")
# #   lines(x2, col = "red")
# # }
# # 
# # 
# # id <- 3
# # plot_autocorr(Mu[, id], MuKu[, id])
# # effectiveSize(Mu)
# # effectiveSize(MuKu)
# # 
# # 
# # 
