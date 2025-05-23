# This file contains all the post-processing functions for the methods available
# Postprocess functions for the all the methods 

Compute_sensitivity_precision <- function(R_hat, R_true, cos_cutoff = 0.9){
  # Sensitivity: proportion of ground truth signatures that were estimated with 
  #              sufficiently high similarity
  sig_sens <- sapply(1:ncol(R_true), function(i) 
    max(sapply(1:ncol(R_hat), function(j) cosine(R_true[, i], R_hat[, j]))))
  # Precision: proportion of estimated signatures that were sufficiently similar to 
  #            ground truth signatures
  sig_prec <- sapply(1:ncol(R_hat), function(i) 
    max(sapply(1:ncol(R_true), function(j) cosine(R_true[, j], R_hat[, i]))))
  return(c("Sensitivity" = mean(sig_sens > cos_cutoff), "Precision" = mean(sig_prec > cos_cutoff)))
}


# Function to find the true value of Lambda
get_Lambda_Comp <- function(resComp,  samples = NULL) {
  Lambda <- 0
  id <- resComp$selected_chain
  if(is.null(samples)){
      nsamples <- dim(resComp$mcmc_out[[id]]$Signatures)[1]
      samples <- 1:nsamples
  }
  for (i in samples) {
    Lambda <- Lambda + resComp$mcmc_out[[id]]$Signatures[i, , ] %*% resComp$mcmc_out[[id]]$Weights[i, , ]
  }
  return(Lambda / length(samples))
}


get_Lambda_CUSP <- function(resCUSP) {
  Lambda <- 0
  nsamples <- length(resCUSP$Signatures)
  for (i in 1:nsamples) {
    Lambda <- Lambda + resCUSP$Signatures[[i]] %*% resCUSP$Weights[[i]]
  }
  return(Lambda / nsamples)
}

# Calculate the cosine similarity between two sets of signatures
get_cosine_similarity <- function(matchedSign){
  sims <- sapply(1:ncol(matchedSign$R_hat), function(i) cosine(matchedSign$R_hat[, i], matchedSign$R_true[, i]))
  return(sims) 
}

# Match the mutational signatures with the true ones using the Hungarian algorithm 
# to calculate the average cosine similarity 
match_MutSign <- function(R_true, R_hat) {
  # Need to make sure that R_hat and R_true are the same
  k_true <- ncol(R_true)
  k_hat <- ncol(R_hat)
  k_tot <- max(c(k_hat, k_true))
  I <- nrow(R_true)
  mat0 <- matrix(100, nrow = I, ncol = abs(k_hat - k_true)) # 
  if (k_hat > k_true) {
    colnames(mat0) <- paste0("new_extra", 1:ncol(mat0))
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
  R_hat_matched <- R_hat[, match]
  
  # Change 100 with 0s
  R_hat_matched[R_hat_matched == 100] <- 0
  R_true[R_true == 100] <- 0
  
  return(list("R_hat" = R_hat_matched, "R_true" = R_true, "match" = match))
}

# Calculate the MSE between the true Theta and the estimated one
match_Theta <- function(Theta_true, Theta_hat, match){
  # Need to make sure that Theta_hat and Theta_true are the same
  k_true <- nrow(Theta_true)
  k_hat <- nrow(Theta_hat)
  k_tot <- max(c(k_hat, k_true))
  J <- ncol(Theta_true)
  mat0 <- matrix(0, nrow = abs(k_hat - k_true), ncol = J)
  if (k_hat > k_true) {
    rownames(mat0) <- paste0("new_extra", 1:nrow(mat0))
    Theta_true <- rbind(Theta_true, mat0)
  } else if (k_hat < k_true) {
    Theta_hat <- rbind(Theta_hat, mat0)
  }
  Theta_hat <- Theta_hat[match, ]
  return(Theta_hat)
}
# Calculate the MSE between the true Theta and the estimated one
compute_RMSE_Theta <- function(Theta_true, Theta_hat, match){
  # Need to make sure that Theta_hat and Theta_true are the same
  k_true <- nrow(Theta_true)
  k_hat <- nrow(Theta_hat)
  k_tot <- max(c(k_hat, k_true))
  J <- ncol(Theta_true)
  mat0 <- matrix(0, nrow = abs(k_hat - k_true), ncol = J)
  if (k_hat > k_true) {
    rownames(mat0) <- paste0("new_extra", 1:nrow(mat0))
    Theta_true <- rbind(Theta_true, mat0)
  } else if (k_hat < k_true) {
    Theta_hat <- rbind(Theta_hat, mat0)
  }
  Theta_hat <- Theta_hat[match, ]
  # Now, match the matrices
  is_novel <- grepl("new", rownames(Theta_true))
  # RMSE between true cosmic 
  if(any(!is_novel)){
    if(sum(!is_novel) == 1){
      Theta_temp <- t(Theta_hat[!is_novel, ])
      Theta_temp_true <- t(Theta_true[!is_novel, ])
    } else {
      Theta_temp <- as.matrix(Theta_hat[!is_novel, ])
      Theta_temp_true <- as.matrix(Theta_true[!is_novel, ])
    }
    rmse_cosmic <- mean(sapply(1:sum(!is_novel), function(i) 
      sqrt(mean((Theta_temp[i, ] - Theta_temp_true[i, ])^2))))
  } else {
    rmse_cosmic  <- -1 
  }
  # RMSE between novel
  if(any(is_novel)){
    if(sum(is_novel) == 1){
      Theta_temp <- t(Theta_hat[is_novel, ])
      Theta_temp_true <- t(Theta_true[is_novel, ])
    } else {
      Theta_temp <- as.matrix(Theta_hat[is_novel, ])
      Theta_temp_true <- as.matrix(Theta_true[is_novel, ])
    }
    rmse_novel <- mean(sapply(1:sum(is_novel), function(i) 
      sqrt(mean(Theta_temp[i, ] - Theta_temp_true[i, ])^2)))
  } else {
    rmse_novel  <- -1 
  }
  rmse_all <- c(rmse_novel, rmse_cosmic)
  
  return(c("rmse_novel_Theta" = rmse_novel, "rmse_cosmic_Theta" = rmse_cosmic, 
           "rmse_Weights" = mean(rmse_all[rmse_all!=-1])))
  return(rmse)
}

compute_RMSE_Signature <- function(R_hat, R_true){
  is_novel <- grepl("new", colnames(R_true))
  # RMSE between true cosmic 
  if(any(!is_novel)){
    rmse_cosmic <- mean(sapply(1:sum(!is_novel), function(i) 
      sqrt(mean((as.matrix(R_hat[, !is_novel])[, i] - as.matrix(R_true[, !is_novel])[, i])^2))))
  } else {
    rmse_cosmic  <- -1 
  }
  # RMSE between novel
  if(any(is_novel)){
    rmse_novel <- mean(sapply(1:sum(is_novel), function(i) 
      sqrt(mean((as.matrix(R_hat[, is_novel])[, i] - as.matrix(R_true[, is_novel])[, i])^2))))
  } else {
    rmse_novel  <- -1 
  }
  rmse_all <- c(rmse_novel, rmse_cosmic)
  
  return(c("rmse_novel_Sig" = rmse_novel, "rmse_cosmic_Sig" = rmse_cosmic, 
           "rmse_Signatures" = mean(rmse_all[rmse_all!=-1])))
}


get_ESS_Comp <- function(resComp){
  
  R_all <- resComp$mcmc_out[[resComp$selected_chain]]$Signatures
  Theta_all <- resComp$mcmc_out[[resComp$selected_chain]]$Weights
  Mu_all <- resComp$mcmc_out[[resComp$selected_chain]]$Mu
  select <- colMeans(Mu_all) > 0.005
  
  EffectiveSigs <- apply(R_all[, , select], c(2,3), function(x) coda::effectiveSize(x))
  EffectiveTheta <- apply(Theta_all[, select, ], c(2,3), function(x) coda::effectiveSize(x))
  EffectiveRelW <- coda::effectiveSize(Mu_all[, select])
  effsize <- c("ESS_Sig_mean" = mean(colMeans(EffectiveSigs)), 
               "ESS_Sig_sd" = sd(colMeans(EffectiveSigs)),
               "ESS_Theta_mean" = mean(rowMeans(EffectiveTheta)), 
               "ESS_Theta_sd" = sd(rowMeans(EffectiveTheta)), 
               "ESS_relweight_mean" = mean(EffectiveRelW), 
               "ESS_relweight_sd" = sd(EffectiveRelW))
  
  return(effsize)
}

get_ESS_signeR <- function(resSigneR, which_samples = 1001:2000){
  
  R_all <- resSigneR$SignExposures@Sign[, , which_samples]
  Theta_all <- resSigneR$SignExposures@Exp[, , which_samples]
  EffectiveSigs <- apply(R_all, c(1,2), function(x) coda::effectiveSize(x))
  EffectiveTheta <- apply(Theta_all, c(1, 2), function(x) coda::effectiveSize(x))
  effsize <- c("ESS_Sig_mean" = mean(colMeans(EffectiveSigs)), 
               "ESS_Sig_sd" = sd(colMeans(EffectiveSigs)),
               "ESS_Theta_mean" = mean(rowMeans(EffectiveTheta)), 
               "ESS_Theta_sd" = sd(rowMeans(EffectiveTheta)), 
               "ESS_relweight_mean" = NA, 
               "ESS_relweight_sd" = NA)
  
  return(effsize)
}

get_ESS_PoissonCUSP <- function(out_CUSP){
  pos_cusp <- get_posterior_CUSP(out_CUSP)
  Mu_all <- pos_cusp$mu_seq[, pos_cusp$nspike < 0.05]
  sigMat <- pos_cusp$R_hat[, pos_cusp$nspike < 0.05]
  ThetaMat <- pos_cusp$Theta_hat[pos_cusp$nspike < 0.05, ]
  nsamples <- nrow(pos_cusp$mu_seq)
  R_all <- array(NA, dim = c(nsamples, nrow(sigMat), ncol(sigMat)))
  Theta_all <- array(NA, dim = c(nsamples, ncol(sigMat), ncol(ThetaMat)))
  for(i in 1:nsamples){
    # match the signature with closest similarity
    cosMat <- sigminer::cosine(sigMat, out_CUSP$Signatures[[i]])
    match <- RcppHungarian::HungarianSolver(cosMat)$pairs[, 1]
    R_all[i, ,] <- out_CUSP$Signatures[[i]][, match]
    Theta_all[i, , ] <- out_CUSP$Weights[[i]][match, ]
  }
  # Calculate effective sample sizes
  EffectiveSigs <- apply(R_all, c(2,3), function(x) coda::effectiveSize(x))
  EffectiveTheta <- apply(Theta_all, c(2,3), function(x) coda::effectiveSize(x))
  EffectiveRelW <- coda::effectiveSize(Mu_all)
  effsize <- c("ESS_Sig_mean" = mean(colMeans(EffectiveSigs)), 
               "ESS_Sig_sd" = sd(colMeans(EffectiveSigs)),
               "ESS_Theta_mean" = mean(rowMeans(EffectiveTheta)), 
               "ESS_Theta_sd" = sd(rowMeans(EffectiveTheta)), 
               "ESS_relweight_mean" = mean(EffectiveRelW), 
               "ESS_relweight_sd" = sd(EffectiveRelW))
  return(effsize)
}

#---------------------------------------------------------------------- CompressiveNMF
Postprocess_Compressive <- function(resComp, data) {
  # Step 1 - calculate the number of inferred signatures
  K <- length(resComp$RelWeights)
  # Step 2 - Calculate RMSE wrt to the count matrix and the rate matrix
  Lambda <- get_Lambda_Comp(resComp)
  Lambda_true <- data$Rmat %*% data$Theta
  rmse_Lambda <- sqrt(mean((Lambda_true - Lambda)^2))
  rmse_Counts <- sqrt(mean((data$X - Lambda)^2))
  # Step 3 - calculate the cosine similarity between the true and the inferred signatures
  R_true <- apply(data$Rmat, 2, function(x) x / sum(x))
  R_hat <- resComp$Signatures
  matchedSign <- match_MutSign(R_true = R_true, R_hat = R_hat)
  cos_sim <- mean(get_cosine_similarity(matchedSign))
  # Step 4 - calculate the RMSE between Theta and the rest
  Theta_hat <- resComp$Weights
  rmse_R <- compute_RMSE_Signature(R_hat = matchedSign$R_hat, R_true = matchedSign$R_true)  
  rmse_Theta <- compute_RMSE_Theta(Theta_true = data$Theta, Theta_hat = Theta_hat, matchedSign$match)  
  # Step 5 - calculate the sensitivity and precision
  sens_prec  <- Compute_sensitivity_precision(R_hat = R_hat, R_true)
  # Step 6 - add Effective sample sizes
  effsize <- get_ESS_Comp(resComp)
  return(list(
    Lambda = Lambda,
    R_hat = R_hat,
    Theta_hat = Theta_hat,
    Mu_hat = resComp$RelWeights,
    signatures = matchedSign,
    results = c("K" = K, 
                "rmse_Lambda" = rmse_Lambda, 
                "rmse_Counts" = rmse_Counts, 
                rmse_R, 
                rmse_Theta, 
                sens_prec,
                "cos_sim" = cos_sim, 
                "time" = resComp$time, 
                effsize)
  ))
}

#---------------------------------------------------------------------- SignatureAnalyzer
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
  # Step 4 - calculate the RMSE between Theta and the rest
  rmse_R <- compute_RMSE_Signature(R_hat = matchedSign$R_hat, R_true = matchedSign$R_true)  
  rmse_Theta <- compute_RMSE_Theta(Theta_true = data$Theta, Theta_hat = resARD$H, matchedSign$match)  
  # Step 5 - calculate the sensitivity and precision
  sens_prec  <- Compute_sensitivity_precision(R_hat = resARD$W, data$Rmat)
  return(list(
    "Lambda" = Lambda,
    "R_hat" = resARD$W,
    "Theta_hat" = resARD$H,
    #"Mu_hat" = resARD$lambda,
    "signatures" = matchedSign,
    "results" = c("K" = K, 
                  "rmse_Lambda" = rmse_Lambda, 
                  "rmse_Counts" = rmse_Counts, 
                  rmse_R, 
                  rmse_Theta, 
                  sens_prec,
                  "cos_sim" = cos_sim, 
                  "time" = resARD$time)
  ))
}



#---------------------------------------------------------------------- SigneR
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
  # RMSE
  rmse_R <- compute_RMSE_Signature(R_hat = matchedSign$R_hat, R_true = matchedSign$R_true)  
  rmse_Theta <- compute_RMSE_Theta(Theta_true = data$Theta, Theta_hat = resSigneR$Ehat, matchedSign$match)
  # Step 5 - calculate the sensitivity and precision
  sens_prec  <- Compute_sensitivity_precision(R_hat = resSigneR$Phat, data$Rmat)
  # Step 6 - add Effective sample sizes
  effsize <- get_ESS_signeR(resSigneR)
  return(list(
    "Lambda" = Lambda,
    "R_hat" = resSigneR$Phat,
    "Theta_hat" = resSigneR$Ehat,
    "signatures" = matchedSign,
    "results" = c("K" = K, 
                  "rmse_Lambda" = rmse_Lambda, 
                  "rmse_Counts" = rmse_Counts, 
                  rmse_R, 
                  rmse_Theta, 
                  sens_prec,
                  "cos_sim" = cos_sim, 
                  "time" = resSigneR$time, 
                  effsize)
  ))
}

#---------------------------------------------------------------------- SigProfiler
Postprocess_SigProfiler <- function(resSigPro, data) {
  # Step 1 - find the number of signatures
  K <- resSigPro$K
  # Step 2 - Calculate RMSE wrt to the count matrix and the rate matrix
  Lambda <- resSigPro$Signatures %*% resSigPro$Weights
  Lambda_true <- data$Rmat %*% data$Theta
  rmse_Lambda <- sqrt(mean((Lambda_true - Lambda)^2))
  rmse_Counts <- sqrt(mean((data$X - Lambda)^2))
  # Step 3 - Infer the mutational signatures
  matchedSign <- match_MutSign(R_true = data$Rmat, R_hat = resSigPro$Signatures)
  cos_sim <- mean(get_cosine_similarity(matchedSign))
  # RMSE
  rmse_R <- compute_RMSE_Signature(R_hat = matchedSign$R_hat, R_true = matchedSign$R_true)  
  rmse_Theta <- compute_RMSE_Theta(Theta_true = data$Theta, Theta_hat = resSigPro$Weights, matchedSign$match)
  # Step 5 - calculate the sensitivity and precision
  sens_prec  <- Compute_sensitivity_precision(R_hat = resSigPro$Signatures, data$Rmat)
  return(list(
    "Lambda" = Lambda,
    "R_hat" = resSigPro$Signatures,
    "Theta_hat" = resSigPro$Weights,
    "signatures" = matchedSign,
    "results" = c("K" = K, 
                  "rmse_Lambda" = rmse_Lambda, 
                  "rmse_Counts" = rmse_Counts, 
                  rmse_R, 
                  rmse_Theta, 
                  sens_prec,
                  "cos_sim" = cos_sim, 
                  "time" = resSigPro$time)
  ))
}



get_posterior_CUSP <- function(resCUSP) {
  J <- ncol(resCUSP$Weights[[1]])
  Kmax <- max(unlist(lapply(resCUSP$Mu, length)))
  nspike <- rep(0, Kmax)
  nsims <- length(resCUSP$Signatures)
  
  k0 <- length(resCUSP$Mu[[1]])
  R_hat <- cbind(resCUSP$Signatures[[1]] , matrix(1 / 96, nrow = 96, ncol = Kmax - k0))
  Theta_hat <- rbind(resCUSP$Weights[[1]], matrix(0, ncol = J, nrow =  Kmax - k0))
  Mu_hat <- c(resCUSP$Mu[[1]], rep(resCUSP$spike, Kmax - k0))
  mu_seq <- t(c(resCUSP$Mu[[1]], rep(resCUSP$spike, Kmax - k0)))
  
  for (i in 2:nsims) {
    mu_temp <- resCUSP$Mu[[i]]; mu_old <- resCUSP$Mu[[i-1]]
    # Check if a column was dropped
    k_new <- length(mu_temp); k_old <- length(mu_old)
    if(k_new < k_old){
      # Check which column was dropped
      id_dropped <- which(mu_old == resCUSP$spike) # These are the columns that were dropped. We need to take the resulting mean and move it to the end
      R_hat <- cbind(R_hat[, -id_dropped], R_hat[, id_dropped])
      Theta_hat <- rbind(Theta_hat[-id_dropped, ], Theta_hat[id_dropped, ])
      Mu_hat <- c(Mu_hat[-id_dropped], Mu_hat[id_dropped])
      if(i == 2){
        mu_seq <- cbind(t(mu_seq[, -id_dropped]), t(mu_seq[, id_dropped]))  
      } else {
        mu_seq <- cbind(mu_seq[, -id_dropped], mu_seq[, id_dropped])
      }
      nspike <- c(nspike[-id_dropped], nspike[id_dropped])
    }
    # Check global mean
    mu <- c(mu_temp, rep(resCUSP$spike, Kmax - length(mu_temp)))
    mu_seq <- rbind(mu_seq, mu)
    Mu_hat <- Mu_hat + mu
    nspike <- nspike + 1 * (mu == resCUSP$spike)
    # Check signatures
    mat_temp <- resCUSP$Signatures[[i]]
    R_hat <- R_hat + cbind(mat_temp, matrix(1 / 96, nrow = 96, ncol = Kmax - length(mu_temp)))
    # Check weights
    Theta_temp <- resCUSP$Weights[[i]]
    Theta_hat <- Theta_hat + rbind(Theta_temp, matrix(0, ncol = ncol(Theta_temp),
                                                      nrow =  Kmax - length(mu_temp)))
    #mu_seq <- rbind(mu_seq, mu)
  }
  R_hat <- R_hat / nsims
  nspike <- nspike / nsims
  Theta_hat <- Theta_hat / nsims
  Mu_hat <- Mu_hat / nsims
  return(list("Theta_hat" = Theta_hat, "R_hat" = R_hat,
              "Mu_hat" = Mu_hat, "nspike" = nspike, "mu_seq" = mu_seq))
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
  # Step 4 - Calculate RMSE for signatures and Weights
  Theta_hat <- post_CUSP$Theta_hat[nspike < 0.05, ]
  rmse_R <- compute_RMSE_Signature(R_hat = matchedSign$R_hat, R_true = matchedSign$R_true)  
  rmse_Theta <- compute_RMSE_Theta(Theta_true = data$Theta, Theta_hat = Theta_hat, matchedSign$match)
  # Step 5 - calculate the sensitivity and precision
  sens_prec  <- Compute_sensitivity_precision(R_hat = R_hat, data$Rmat)
  # Step 6 - add Effective sample sizes
  effsize <- get_ESS_PoissonCUSP(resCUSP)
  return(list("Lambda" = Lambda,
              "R_hat" = R_hat, #post_CUSP$R_hat,
              "Theta_hat" = Theta_hat, #post_CUSP$Theta_hat,
              "Mu_hat" = post_CUSP$Mu_hat,
              "Kchain" = Kchain,
              "Mu_chain"= post_CUSP$mu_seq,
              "nspike" = nspike,
              "signatures" = matchedSign,
              "results" = c("K" = K, 
                            "rmse_Lambda" = rmse_Lambda, 
                            "rmse_Counts" = rmse_Counts, 
                            rmse_R, 
                            rmse_Theta, 
                            sens_prec,
                            "cos_sim" = cos_sim, 
                            "time" = resCUSP$time, 
                            effsize)
         ))
}

#---------------------------------------------------------------------- BayesNMF
Postprocess_BayesNMF <- function(resBayesNMF, data) {
  # Step 1 - calculate the number of inferred signatures
  K <- ncol(resBayesNMF$Signatures)
  # Step 2 - Calculate RMSE wrt to the count matrix and the rate matrix
  Lambda <- resBayesNMF$Signatures %*% resBayesNMF$Theta
  Lambda_true <- data$Rmat %*% data$Theta
  rmse_Lambda <- sqrt(mean((Lambda_true - Lambda)^2))
  rmse_Counts <- sqrt(mean((data$X - Lambda)^2))
  # Step 3 - calculate the cosine similarity between the true and the inferred signatures
  R_true <- apply(data$Rmat, 2, function(x) x / sum(x))
  R_hat <- resBayesNMF$Signatures
  matchedSign <- match_MutSign(R_true = R_true, R_hat = R_hat)
  cos_sim <- mean(get_cosine_similarity(matchedSign))
  # Step 4 - calculate the RMSE between Theta and the rest
  Theta_hat <- resBayesNMF$Theta
  rmse_R <- compute_RMSE_Signature(R_hat = matchedSign$R_hat, R_true = matchedSign$R_true)  
  rmse_Theta <- compute_RMSE_Theta(Theta_true = data$Theta, Theta_hat = Theta_hat, matchedSign$match)  
  # Step 5 - calculate the sensitivity and precision
  sens_prec  <- Compute_sensitivity_precision(R_hat = R_hat, R_true)
  # Step 6 - calculate effective sample sizes of R and Theta on average
  effsize <- c("ESS_Sig_mean" = mean(colMeans(resBayesNMF$EffectiveSigs)), 
               "ESS_Sig_sd" = sd(colMeans(resBayesNMF$EffectiveSigs)),
               "ESS_Theta_mean" = mean(colMeans(resBayesNMF$EffectiveTheta)), 
               "ESS_Theta_sd" = sd(colMeans(resBayesNMF$EffectiveTheta)), 
               "ESS_relweight_mean" = mean(resBayesNMF$EffectiveLambda), 
               "ESS_relweight_sd" = sd(resBayesNMF$EffectiveLambda))
  return(list(
    Lambda = Lambda,
    R_hat = R_hat,
    Theta_hat = Theta_hat,
    Mu_hat = resBayesNMF$RelWeights,
    signatures = matchedSign,
    results = c("K" = K, 
                "rmse_Lambda" = rmse_Lambda, 
                "rmse_Counts" = rmse_Counts, 
                rmse_R, 
                rmse_Theta, 
                sens_prec,
                "cos_sim" = cos_sim, 
                "time" = resBayesNMF$time, 
                effsize)
  ))
}


#---------------------------------------------------------------------- Extract results

add_ess <- function(x) {
  names_ess <- c("ESS_Sig_mean", "ESS_Sig_sd", "ESS_Theta_mean",
                 "ESS_Theta_sd", "ESS_relweight_mean", "ESS_relweight_sd")
  add_na <- rep(NA, length(names_ess))
  names(add_na) <- names_ess
  if(!any(names_ess %in% names(x))){
    x <- c(x, add_na)
  }
  x
}  
  
extract_results <- function(out, name){
  df <- data.frame("Method" = name, "Simulation" = 1:length(out))
  df <- cbind(df, as.data.frame(t(sapply(1:length(out), function(i) add_ess(out[[i]]$results), simplify = TRUE))))
  if(name == "3.CUSP"){
    df$K = unlist(lapply(out, function(x) sum(x$nspike < 0.05)))
  }
  return(df)
}


# results functions

aggregate_results <- function(out_CompressiveNMF = NULL,
                              out_CompressiveNMF_cosmic = NULL,
                              out_PoissonCUSP = NULL, 
                              out_ARD = NULL, 
                              out_signeR = NULL, 
                              out_sigPro = NULL, 
                              out_BayesNMF = NULL) {
  df_res <- NULL
  #---------- CompressiveNMF
  if(!is.null(out_CompressiveNMF)){
    df_res <- rbind(df_res, extract_results(out_CompressiveNMF, name = "2.CompNMF")) 
  }
  #---------- CompressiveNMF with cosmic database 
  if(!is.null(out_CompressiveNMF_cosmic)){
    df_res <- rbind(df_res, extract_results(out_CompressiveNMF_cosmic, name = "1.CompNMFcos")) 
  }
  #---------- PoissonCUSP
  if(!is.null(out_PoissonCUSP)){
    df_res <- rbind(df_res, extract_results(out_PoissonCUSP, name = "6.CUSP")) 
  }
  #---------- SignatureAnalyzer
  if(!is.null(out_ARD)){
    df_res <- rbind(df_res, extract_results(out_ARD, name = "5.ARD")) 
  }
  #---------- signeR
  if(!is.null(out_signeR)){
    df_res <- rbind(df_res, extract_results(out_signeR, name = "3.signeR")) 
  }
  #---------- sigProfiler
  if(!is.null(out_sigPro)){
    df_res <- rbind(df_res, extract_results(out_sigPro, name = "4.SigPro")) 
  }
  #---------- BayesianNMF
  if(!is.null(out_BayesNMF)){
    df_res <- rbind(df_res, extract_results(out_BayesNMF, name = "7.BayesNMF")) 
  }
  return(df_res)
}



################################################################################
# Functions to match to cosmic + uncertainty quantification
################################################################################
match_to_cosmic <- function(R){
  load("~/CompressiveNMF/data/Cosmic_data.rdata")
  df_res <- data.frame()
  if(is.null(colnames(R))){
    colnames(R) <- paste0("Sig", 1:ncol(R))
  }
  for(k in 1:ncol(R)){
    signature <- R[, k]
    dist <- apply(cosmic_data[, -c(1,2,3)], 2, function(x) cosine(signature, x))
    df_res <- rbind(df_res, 
                    data.frame("signature" = colnames(R)[k],
                               "best_cosmic" = names(which.max(dist)),
                               "cosine_sim" = round(max(dist), 3)))
  }
  return(df_res)
}


# CompressiveNMF
match_to_cosmic_uncertainty_CompNMF <- function(out_CompNMF) {
  load("~/CompressiveNMF/data/Cosmic_data.rdata")
  # Find signatures
  sigMat <- out_CompNMF$Signatures
  # Get the best chain
  nchains <- length(out_CompNMF$mcmc_out)
  # post <-lapply(out_CompNMF$mcmc_out, function(x) postprocess_mcmc_out(x, 0.05))
  # logposterior <- sapply(1:nchains, function(i) post[[i]]$logpost)
  # id_best <- which.max(logposterior)
  # # Get the posterior samples for all signatures
  chain <- out_CompNMF$mcmc_out[[out_CompNMF$selected_chain]]
  nonzero_sign <- which(colMeans(chain$Mu) > 0.05)
  R_chain <- chain$Signatures[, , nonzero_sign]
  # Match signatures
  df_match <- match_to_cosmic(sigMat)
  df_match$lowCI <- NA
  df_match$highCI <- NA
  for(s in 1:nrow(df_match)){
    range_cosine <- apply(R_chain[, , s], 1, function(x) cosine(cosmic_data[, df_match$best_cosmic[s]], x))
    df_match$lowCI[s] <- c(quantile(range_cosine, c(0.05)))
    df_match$highCI[s] <- c(quantile(range_cosine, c(0.95)))
  }
  return(df_match)
}


match_to_cosmic_uncertainty_PoissonCUSP <- function(out_CUSP) {
  load("~/CompressiveNMF/data/Cosmic_data.rdata")
  # Find signatures
  pos_cusp <- get_posterior_CUSP(out_CUSP)
  sigMat <- pos_cusp$R_hat[, pos_cusp$nspike < 0.05]
  # Get the chain. Instead of re-aligning it, we simply get the columns with highest 
  # cosine similarity in each sample
  df_match <- match_to_cosmic(sigMat)
  df_match$lowCI <- NA
  df_match$highCI <- NA
  for(s in 1:nrow(df_match)){
    range_cosine <- unlist(lapply(out_CUSP$Signatures, function(Xmat) max(apply(Xmat, 2, function(x) cosine(cosmic_data[, df_match$best_cosmic[s]], x)))))
    df_match$lowCI[s] <- c(quantile(range_cosine, c(0.05)))
    df_match$highCI[s] <- c(quantile(range_cosine, c(0.95)))
  }
  return(df_match)
}

match_to_cosmic_uncertainty_signeR <- function(out_signeR) {
  load("~/CompressiveNMF/data/Cosmic_data.rdata")
  # Find signatures
  sigMat <- out_signeR$Phat
  # Get the chain. 
  df_match <- match_to_cosmic(sigMat)
  df_match$lowCI <- NA
  df_match$highCI <- NA
  for(s in 1:nrow(df_match)){
    range_cosine <- apply(out_signeR$SignExposures@Sign[, s, ], 2, function(x) cosine(cosmic_data[, df_match$best_cosmic[s]], x))
    df_match$lowCI[s] <- c(quantile(range_cosine, c(0.05)))
    df_match$highCI[s] <- c(quantile(range_cosine, c(0.95)))
  }
  return(df_match)
}


Compute_F1_range <- function(R_hat, R_true, cutoff_range){
  Sens_prec <- sapply(1:length(cutoff_range), function(i) Compute_sensitivity_precision(R_hat, R_true, cos_cutoff = cutoff_range[i]))
  F1 <- apply(Sens_prec, 2, function(x) 2 * prod(x)/sum(x))
  F1[is.nan(F1)] <- 0
  return(F1)
}

extract_F1_range <- function(res, cutoff_range = seq(0.6, 0.9999, length.out = 150), method = ""){
  F1_output <- lapply(res, function(results) {
    if(!is.null(results)){
      # Extract estimated matrices
      R_hat <- results$R_hat
      R_true <- results$signatures$R_true
      # Eliminate the extra columns
      R_true <- R_true[, !grepl("new_extra", colnames(R_true))]
      # Calculate F1
      Compute_F1_range(R_hat, R_true, cutoff_range)
    } 
  })
  F1_output <- do.call("cbind", F1_output)
  res <- t(apply(F1_output, 1, function(x) c("mean" = mean(x), 
                                             "lowCI" = quantile(x, 0.05), 
                                             "median" = quantile(x, 0.5),
                                             "highCI" = quantile(x, 0.95))))
  data.frame("Method" = method,  "cutoff" = cutoff_range, res)
}


