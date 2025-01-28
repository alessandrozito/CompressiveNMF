# This file does the pre- and post-processing for the BayesNMF ARD in 
# Brouwer, Fressner and Lio (2017), which contains a Gibbs sampler for the 
# BayesianNMF with ARD prior that assumes a Gaussian likelihood. 

# Their code is written in Python 2.7 (which is no longer maintained). We download
# the package from their github at https://github.com/ThomasBrouwer/BNMTF_ARD
# and test their model on our simulated data 

# Load libraries
library(tidyverse)
library(CompressiveNMF)
library(coda)
library(lsa)

source("R/Postprocess_functions.R")

# Useful function
create_directory <- function(dir) {
  if (!dir.exists(dir)) {
    dir.create(dir)
  }
}

reshape_python_array <- function(A, dimension = 96, K = 20){
  A_all <- array(NA, dim = c(nrow(A)/dimension, dimension, K))
  n <- nrow(A)
  id_all <- seq(1, n, by = dimension)
  for(s in 1:(n/dimension - 1)) A_all[s, , ] <- A[id_all[s]:(id_all[s+1] - 1), ]
  A_all[n/dimension, ,] <- A[id_all[s + 1]:n, ] 
  return(A_all)
}


#------------------------------------------------ 
# Step 0 
#------------------------------------------------
# Load the simulated data (previously stored in .rds file) and save them
# in .txt file for better opening

simulation_dir <- "~/CompressiveNMF/output/main_simulation/"

overdispersion_list <- c(0, 0.15)
J_list <- c(50, 100, 200)
K_new_list <- c(2, 6)
theta <- 100
# overdispersion <- 0
# J <- 50
# K_new <- 2
create_directories <- FALSE # <------ Set to TRUE if run for the first time

if(create_directories) {
  for(K_new in K_new_list){
    for(J in J_list) {
      for(overdispersion in overdispersion_list){
        # Set output directory
        out_dir <- paste0(simulation_dir, "Scenario_", J, "_overdisp_", 
                          overdispersion, "_Knew_", K_new, "_theta_", theta)
        data <- readRDS(paste0(out_dir, "/data.rds.gzip"))
        create_directory(paste0(out_dir, "/BayesNMF_brouwer/"))
        print(out_dir)
        # Re-save the data as .txt
        out_dir_data <- paste0(out_dir, "/BayesNMF_brouwer/data/")
        out_dir_output <- paste0(out_dir, "/BayesNMF_brouwer/output/")
        
        create_directory(out_dir_data)
        create_directory(out_dir_output)
        
        for(i in 1:length(data)){
          # Get the data
          X <- data[[i]]$X
          nm <- paste0(out_dir_data, "X_", i, ".txt")
          write.table(X, file = nm)
          # Prepare the output folder
          nm_out <- paste0(out_dir_output, "sim_", i, "/")
          create_directory(nm_out)
        }
      }
    }
  }
}

#------------------------------------------------ 
# Step 1 - Run in python
#------------------------------------------------

# in the console, in the directory run_BayesNMF_python (where the ARD method stored)
# nohup python simulation.py &

#------------------------------------------------ 
# Step 2 - Postprocess output
#------------------------------------------------

ReshapeResults_BayesNMF_brouwer <- function(i, K_new, J, overdispersion){
  simulation_dir <- "~/CompressiveNMF/output/main_simulation/"
  # Directories
  out_dir <- paste0(simulation_dir, "Scenario_", J, "_overdisp_", 
                    overdispersion, "_Knew_", K_new, "_theta_", theta)
  out_dir_data <- paste0(out_dir, "/BayesNMF_brouwer/data/")
  out_dir_output <- paste0(out_dir, "/BayesNMF_brouwer/output/")
  # Data
  nm <- paste0(out_dir_data, "X_", i, ".txt")
  X <- read.table(file = nm)
  # Model
  nm_out <- paste0(out_dir_output, "sim_", i, "/")
  
  # R matrix  
  R <- as.matrix(read_table(paste0(nm_out, "Signatures.txt"), 
                            col_names = FALSE, show_col_types = FALSE))
  R_all <- reshape_python_array(R, dimension = 96)
  Rmean <- apply(R_all, c(2, 3), mean)
  
  # Theta
  Theta <- as.matrix(read_table(paste0(nm_out, "Loadings.txt"), 
                                col_names = FALSE, show_col_types = FALSE))
  Theta_all <- reshape_python_array(Theta, dimension = J)
  Theta_mean <- apply(Theta_all, c(2, 3), mean)
  
  # Relevance weights
  Lambda <- as.matrix(read_table(paste0(nm_out, "lambda.txt"), 
                                 col_names = FALSE, show_col_types = FALSE))
  
  # time for execution
  time_exec <- as.matrix(read_table(paste0(nm_out, "times.txt"), 
                                 col_names = FALSE, show_col_types = FALSE))
  time_exec <- time_exec[nrow(time_exec), 1]
  
  # Calculate the effective sample sizes
  EffectiveSigs <- apply(R_all, c(2,3), function(x) coda::effectiveSize(x))
  EffectiveTheta <- apply(Theta_all, c(2,3), function(x) coda::effectiveSize(x))
  EffectiveLambda <- coda::effectiveSize(Lambda)
  
  # Select the number of signatures by excluding the ones that are plainly flat
  norm_weight <- colSums(Rmean)
  sigs_norm <- sapply(1:20, function(i) Rmean[, i]/norm_weight[i])
  Theta_norm <- t(sapply(1:20, function(i) Theta_mean[, i]*norm_weight[i]))
  select <- apply(sigs_norm, 2, function(x) cosine(x, rep(1, 96))) < 0.975
  
  return(list(
    Signatures = sigs_norm[, select],
    Theta = Theta_norm[select, ], 
    RelWeights = colMeans(Lambda)[select],
    EffectiveSigs = EffectiveSigs[, select],
    EffectiveTheta = EffectiveTheta[, select], 
    EffectiveLambda = EffectiveLambda[select],
    time = unname(time_exec)
    ))
}

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

# Save the output in the same format for the function Postprocess_*model* in the 
# Postprocess_functions.R
for(K_new in K_new_list){
  for(J in J_list) {
    for(overdispersion in overdispersion_list){
      out_dir <- paste0(simulation_dir, "Scenario_", J, "_overdisp_", 
                        overdispersion, "_Knew_", K_new, "_theta_", theta)
      for(i in 1:20){
        print(c(i, K_new, overdispersion, J))
        results <- ReshapeResults_BayesNMF_brouwer(i, K_new, J, overdispersion)
        nm <- paste0(out_dir, "/BayesNMF_brouwer/results_", i, ".rds.gzip")
        saveRDS(results, file = nm, compress = "gzip")
      }
      data_all <- readRDS(paste0(out_dir, "/data.rds.gzip"))
      # Now, postprocess the output and save it
      res <- lapply(1:20, function(i){
        # Read the result
        nm <- paste0(out_dir, "/BayesNMF_brouwer/results_", i, ".rds.gzip")
        resBayesNMF <- readRDS(nm)
        Postprocess_BayesNMF(resBayesNMF, data_all[[i]])
      })
      # Save the merged output
      saveRDS(res, file = paste0(out_dir, "/BayesNMF_brouwer.rds.gzip"))
    }
  }
}
