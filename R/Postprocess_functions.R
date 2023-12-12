# This file contains all the post-processing functions for the methods available
# Postprocess functions for the all the methods 



get_Lambda_Comp <- function(resComp, fast = FALSE) {
  Lambda <- 0
  if (fast) {
    nsamples <- dim(resComp$Signatures)[3]
    for (i in 1:nsamples) {
      Lambda <- Lambda + resComp$Signatures[, , i] %*% resComp$Weights[, , i]
    }
  } else {
    nsamples <- dim(resComp$Signatures)[1]
    for (i in 1:nsamples) {
      Lambda <- Lambda + resComp$Signatures[i, , ] %*% resComp$Weights[i, , ]
    }
  }
  return(Lambda / nsamples)
}


get_Lambda_CUSP <- function(resCUSP) {
  Lambda <- 0
  nsamples <- length(resCUSP$Signatures)
  for (i in 1:nsamples) {
    Lambda <- Lambda + resCUSP$Signatures[[i]] %*% resCUSP$Weights[[i]]
  }
  return(Lambda / nsamples)
}

get_cosine_similarity <- function(matchedSign){
  sims <- sapply(1:ncol(matchedSign$R_hat), function(i) cosine(matchedSign$R_hat[, i], matchedSign$R_true[, i]))
  return(sims) 
}


match_MutSign <- function(R_true, R_hat) {
  # Need to make sure that R_hat and R_true are the same
  k_true <- ncol(R_true)
  k_hat <- ncol(R_hat)
  k_tot <- max(c(k_hat, k_true))
  I <- nrow(R_true)
  mat0 <- matrix(1 / I, nrow = I, ncol = abs(k_hat - k_true))
  if (k_hat > k_true) {
    R_true <- cbind(R_true, mat0)
  } else if (k_hat < k_true) {
    R_hat <- cbind(R_hat, mat0)
  }

  # Match mutational signatures using the Hungarian algorithm
  CosMat <- matrix(1, k_tot, k_tot)
  for (i in 1:k_tot) {
    for (j in 1:k_tot) {
      CosMat[i, j] <- 1 - cosine(R_true[, i], R_hat[, j])
    }
  }
  match <- RcppHungarian::HungarianSolver(CosMat)$pairs[, 2]
  return(list("R_hat" = R_hat[, match], "R_true" = R_true, "match" = match))
}


Postprocess_Compressive <- function(resComp, data, cutoff = 0.03, fast = FALSE) {
  id_drop <- colMeans(resComp$Mu) > cutoff
  # Step 1 - calculate the number of inferred signatures
  K <- sum(id_drop)
  # Step 2 - Calculate RMSE wrt to the count matrix and the rate matrix
  Lambda <- get_Lambda_Comp(resComp, fast = fast)
  Lambda_true <- data$Rmat %*% data$Theta
  rmse_Lambda <- sqrt(mean((Lambda_true - Lambda)^2))
  rmse_Counts <- sqrt(mean((data$X - Lambda)^2))
  # Step 3 - calculate the cosine similarity between the true and the inferred signatures
  R_true <- apply(data$Rmat, 2, function(x) x / sum(x))
  if (fast) {
    R_hat <- apply(resComp$Signatures, c(1, 2), mean)[, id_drop]
  } else {
    R_hat <- apply(resComp$Signatures, c(2, 3), mean)[, id_drop]
  }
  matchedSign <- match_MutSign(R_true = R_true, R_hat = R_hat)
  cos_sim <- mean(get_cosine_similarity(matchedSign))
  return(list(
    Lambda = Lambda,
    R_hat = R_hat,
    Theta_hat = apply(resComp$Weights, c(2,3), mean),
    Mu_hat = colMeans(resComp$Mu),
    signatures = matchedSign,
    results = c("K" = K, "rmse_Lambda" = rmse_Lambda, "rmse_Counts" = rmse_Counts, "cos_sim" = cos_sim)
  ))
}


Postprocess_ARD <- function(resARD, data) {
  # Step 1 - find the number of signatures
  K <- ncol(resARD$W)
  # Step 2 - Calculate RMSE wrt to the count matrix and the rate matrix
  Lambda <- resARD$W %*% resARD$H
  Lambda_true <- data$Rmat %*% data$Theta
  rmse_Lambda <- sqrt(mean((Lambda_true - Lambda)^2))
  rmse_Counts <- sqrt(mean((data$X - Lambda)^2))
  # Step 3 - Calculate cosine similarity
  matchedSign <- match_MutSign(R_true = data$Rmat, R_hat = resARD$W)
  cos_sim <- mean(get_cosine_similarity(matchedSign))
  return(list(
    "Lambda" = Lambda,
    "R_hat" = resARD$W,
    "Theta_hat" = resARD$H,
    "Mu_hat" = resARD$lambda,
    "signatures" = matchedSign,
    "results" = c("K" = K, "rmse_Lambda" = rmse_Lambda, "rmse_Counts" = rmse_Counts, "cos_sim" = cos_sim)
  ))
}



Postprocess_signeR <- function(resSigneR, data) {
  # Step 1 - find the number of signatures
  K <- resSigneR$Nsign
  # Step 2 - Calculate RMSE wrt to the count matrix and the rate matrix
  Lambda <- resSigneR$Phat %*% resSigneR$Ehat
  Lambda_true <- data$Rmat %*% data$Theta
  rmse_Lambda <- sqrt(mean((Lambda_true - Lambda)^2))
  rmse_Counts <- sqrt(mean((data$X - Lambda)^2))
  # Step 3 - Infer the mutational signatures
  matchedSign <- match_MutSign(R_true = data$Rmat, R_hat = resSigneR$Phat)
  cos_sim <- mean(get_cosine_similarity(matchedSign))
  return(list(
    "Lambda" = Lambda,
    "R_hat" = resSigneR$Phat,
    "Theta_hat" = resSigneR$Ehat,
    "signatures" = matchedSign,
    "results" = c("K" = K, "rmse_Lambda" = rmse_Lambda, "rmse_Counts" = rmse_Counts, "cos_sim" = cos_sim)
  ))
}


# resCUSP <- PoissonCUSP(X = data$X,
#                        K = K, nsamples = nsamples,
#                        burnin = burnin, a0 = 0.5, b0 = 0.5, alpha_sp = ncol(data$Rmat))

get_posterior_CUSP <- function(resCUSP) {
  Kmax <- max(unlist(lapply(resCUSP$Mu, length)))
  nspike <- rep(0, Kmax)
  nsims <- length(resCUSP$Signatures)
  R_hat <- 0
  Theta_hat <- 0
  Mu_hat <- 0
  for (i in 1:nsims) {
    # Check global mean
    mu_temp <- resCUSP$Mu[[i]]
    mu <- c(mu_temp, rep(0.01, Kmax - length(mu_temp)))
    Mu_hat <- Mu_hat + mu
    nspike <- nspike + 1 * (mu == 0.01)
    # Check signatures
    mat_temp <- resCUSP$Signatures[[i]]
    R_hat <- R_hat + cbind(mat_temp, matrix(1 / 96, nrow = 96,
                                            ncol = Kmax - length(mu_temp)))
    # Check weights
    Theta_temp <- resCUSP$Weights[[i]]
    Theta_hat <- Theta_hat + rbind(Theta_temp, matrix(0, ncol = ncol(Theta_temp),
                                                      nrow =  Kmax - length(mu_temp)))
  }
  R_hat <- R_hat / nsims
  nspike <- nspike / nsims
  Theta_hat <- Theta_hat / nsims
  Mu_hat <- Mu_hat / nsims
  return(list("Theta_hat" = Theta_hat, "R_hat" = R_hat,
              "Mu_hat" = Mu_hat, "nspike" = nspike))
}


Postprocess_PoissonCUSP <- function(resCUSP, data) {
  # Step 1 - find the number of signatures
  Kchain <- unlist(lapply(resCUSP$Mu, function(x) sum(x != 0.01)))
  K <- mean(Kchain)
  # Step 2 - Calculate RMSE wrt to the count matrix and the rate matrix
  Lambda <- get_Lambda_CUSP(resCUSP)
  Lambda_true <- data$Rmat %*% data$Theta
  rmse_Lambda <- sqrt(mean((Lambda_true - Lambda)^2))
  rmse_Counts <- sqrt(mean((data$X - Lambda)^2))
  # Step 3 - Infer the mutational signatures
  R_true <- apply(data$Rmat, 2, function(x) x / sum(x))
  post_CUSP <- get_posterior_CUSP(resCUSP)
  nspike <- post_CUSP$nspike
  R_hat <- post_CUSP$R_hat[, nspike < 0.05]
  matchedSign <- match_MutSign(R_true = R_true, R_hat = R_hat)
  cos_sim <- mean(get_cosine_similarity(matchedSign))
  return(list("Lambda" = Lambda,
         "R_hat" = post_CUSP$R_hat,
         "Theta_hat" = post_CUSP$Theta_hat,
         "Mu_hat" = post_CUSP$Mu_hat,
         "Kchain" = Kchain,
         "nspike" = nspike,
         "signatures" = matchedSign,
         "results" = c("K" = K, "rmse_Lambda" = rmse_Lambda, "rmse_Counts" = rmse_Counts, "cos_sim" = cos_sim)))
}

