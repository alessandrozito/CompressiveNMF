# This file runs a sensitivity analysis for the 21 breast cancer in 
# CompressiveNMF and CompressiveNMF + cosmic.

library(CompressiveNMF)
library(doParallel)
library(tidyverse)
library(foreach)
library(lsa)

source("R/Postprocess_functions.R")

# Load the 21 breast cancer
load("data/21breast.rdata")


# Possible range of values for the sensitivity analysis of the 21 breast cancer
a_range <- c(1, 2, 0.5)
alpha_range <- c(0.5, 1, 2)
epsilon <- c(0.01, 0.1, 0.25)
K_range <- c(15, 20, 25, 40)
values_all <- expand_grid(a_range, alpha_range, epsilon, K_range)
  
# Parameters of the simulation
nsamples <- 2000
burnin <- 10000

# Find the default solution (we use one chain)
set.seed(10)
out_default <- CompressiveNMF(X = X, K = 15, 
                              nsamples = nsamples, burnin = burnin, 
                              epsilon = 0.01)

# Run the solution across a grid of values in parallel
registerDoParallel(parallel::detectCores() - 1)
set.seed(10, kind = "L'Ecuyer-CMRG")
df_results <- foreach(i = 1:nrow(values_all), .combine = "rbind") %dopar%{
  
  # Parameters
  K <- values_all$K_range[i]
  a <- values_all$a_range[i]
  alpha <- values_all$a_range[i]
  epsilon <- values_all$epsilon[i]
  
  # Run the model
  out_temp <- CompressiveNMF(X = X, 
                             K = K, 
                             a = a, 
                             alpha = alpha,
                             nsamples = nsamples, 
                             burnin = burnin, 
                             epsilon = epsilon)
  
  # Compare the output with the baseline
  K_est <- ncol(out_temp$Signatures)
  sens_prec <- Compute_sensitivity_precision(R_hat = out_temp$Signatures, 
                                             R_true = out_default$Signatures, 
                                             cos_cutoff = 0.95)
  matched_sig <- match_MutSign(R_true = out_default$Signatures, 
                               R_hat = out_temp$Signatures)
  # RMSE_sig
  rmse_sig <- compute_RMSE_Signature(matched_sig$R_hat, matched_sig$R_true)
  # RMSE_weights
  rmse_loadings <- compute_RMSE_Theta(matched_sig$R_hat, matched_sig$R_true)
  # Return the value
  c(K = K, a = a, alpha = alpha, epsilon = epsilon, "Kest" = K_est, sens_prec, 
    rmse_sig, rmse_loadings)
}

df_results <- as.data.frame(df_results)
write_csv(df_results, file = "output/Application_21brca/Sensitivity_21brca_CompNMF.csv")


