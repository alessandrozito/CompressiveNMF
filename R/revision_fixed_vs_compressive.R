# Simulation - impact of strength-matching vs fixed choice of hyperprior. 
# This will show the benefit of having the compressive property

library(CompressiveNMF)
library(tidyverse)
library(mcmcse)
library(foreach)
library(lsa)
library(doParallel)

source("R/Postprocess_functions.R")
source("R/simulate_data.R")

# Useful functions
create_directory <- function(dir) {
  if (!dir.exists(dir)) {
    dir.create(dir)
  }
}

open_rds_file <- function(file){
  if(file.exists(file)){
    out <- readRDS(file)
  } else {
    out <- NULL
  }
  return(out)
}

KL.div <- function(X, R, Theta){
  Xest <- R %*% Theta
  sum(X * (log(X + 1e-8) - log(Xest + 1e-8)) - X + Xest)
}


# Hyperparameters of the simulation *(fixed vs compressiv)
epsilon <- 0.01
strength_fixed <- 10
K <- 20
ndatasets <- 20
K_new <- 6


# Range of sample size values 
J_range <- c(20, 50, 100, 200, 300, 400, 500)
overd_list <-c(0, 0.15)

# Simulate the array of data
regenerate_data <- FALSE

main_dir <- "~/CompressiveNMF/output/compressive_vs_fixed_simulation/"
create_directory(main_dir)

if(regenerate_data){
  
  for(overd in overd_list){
    set.seed(10)
    for(j in 1:length(J_range)){
      J <- J_range[j]
      data_all <- lapply(1:ndatasets, function(i) 
        simulate_data(J = J, 
                      K_new = K_new, 
                      overdispersion = 0))
      directory <- paste0(main_dir, "Scenario_", J, "_overd_", overd, "/")
      create_directory(directory)
      # Save data
      saveRDS(data_all, file = paste0(directory, "data.rds.gzip"), compress = "gzip")
    }
  }
}



# Function to run the simulation
run_simulation_fixed_vs_compressive <- function(J, overd, case, 
                                                epsilon = 0.01, 
                                                strength_fixed = 10, 
                                                nsamples = 500, 
                                                nburn = 3000){
  
  subdir <- paste0(main_dir, "Scenario_", J, "_overd_", overd, "/")
  # Load the data
  data_all <- readRDS(paste0(subdir, "data.rds.gzip"))
  # Run all models in parallel
  registerDoParallel(length(data_all))
  output <- foreach(j = 1:length(data_all), .combine = "rbind") %dopar% {
    # Load data
    data <- data_all[[j]]
    
    # Estimate the model
    if(case == "compressive"){
      res <- CompressiveNMF(X = data$X,
                           a0 = ncol(data$X) + 1, 
                           b0 = epsilon * ncol(data$X), 
                           cutoff_excluded = 0,
                           K = K, 
                           epsilon = epsilon,
                           nsamples = nsamples, 
                           burnin = nburn, 
                           verbose = FALSE)
      
    } else if (case == "fixed") {
      res <- CompressiveNMF(X = data$X,
                            a0 = strength_fixed + 1, 
                            b0 = epsilon * strength_fixed, 
                            cutoff_excluded = 0,
                            K = K, 
                            epsilon = epsilon,
                            nsamples = nsamples, 
                            burnin = nburn, 
                            verbose = FALSE)
    }
    
    # Get results
    Lambda <- get_Lambda_Comp(res)
    Lambda_true <- data$Rmat %*% data$Theta
    rmse_Lambda <- sqrt(mean((Lambda_true - Lambda)^2))
    rmse_Counts <- sqrt(mean((data$X - Lambda)^2))
    # Find estimated signatures
    id_selected <- res$RelWeights > 5 * epsilon
    K_selected <- sum(id_selected)
    # Compare true vs estimated
    R_true <- apply(data$Rmat, 2, function(x) x / sum(x))
    R_hat <- res$Signatures[, id_selected]
    matchedSign <- match_MutSign(R_true = R_true, R_hat = R_hat)
    cos_sim <- mean(get_cosine_similarity(matchedSign))
    # Step 4 - calculate the RMSE between Theta and the rest
    Theta_hat <- res$Weights[id_selected, ]
    rmse_R <- compute_RMSE_Signature(R_hat = matchedSign$R_hat, R_true = matchedSign$R_true)  
    rmse_Theta <- compute_RMSE_Theta(Theta_true = data$Theta, Theta_hat = Theta_hat, matchedSign$match)  
    # Step 5 - calculate the sensitivity and precision
    sens_prec  <- Compute_sensitivity_precision(R_hat = R_hat, R_true)
    # Step 6 - look also at Kl divergence
    kl_counts <- KL.div(data$X, R_hat, Theta_hat)
    
    # Results
    results <- c("Kest" = K_selected,
                      "mean" = mean(res$RelWeights[!id_selected]), 
                      "lowCI" = quantile(res$RelWeights[!id_selected], probs = c(0.05)), 
                      "highCI" = quantile(res$RelWeights[!id_selected], probs = c(0.95)), 
                      sens_prec, 
                      rmse_R, 
                      rmse_Theta, 
                      "mean_cos_sim" = cos_sim, 
                      "rmse_Lambda" = rmse_Lambda, 
                      "rmse_Counts" = rmse_Counts,
                      "kl_counts" = kl_counts)
    
    results
  }
  
  # Add information on the simulation
  output<- as.data.frame(output)
  output$overd = epsilon
  output$J = J
  output$overd = overd
  output$case = case
  
  # Save the output
  write_csv(output, file = paste0(subdir, "results_", case, ".csv"))
  
}


# Sampler parameters
nsamples <- 10
nburn <- 30

#-------------------------------------------- Run compressive case
run_compressive <- FALSE
if(run_compressive){
  for(J in J_range){
    for(overd in overd_list) {
      run_simulation_fixed_vs_compressive(J = J, 
                                          overd = overd, 
                                          case = "compressive", 
                                          nsamples = nsamples, 
                                          nburn = nburn)
    }
  }
}

#-------------------------------------------- Run fixed case
run_fix <- FALSE
if(run_fix){
  for(J in J_range){
    for(overd in overd_list) {
      run_simulation_fixed_vs_compressive(J = J, 
                                          overd = overd, 
                                          case = "fixed", 
                                          nsamples = nsamples, 
                                          nburn = nburn)
    }
  }
}








