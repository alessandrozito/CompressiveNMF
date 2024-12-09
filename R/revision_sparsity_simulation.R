# This file evaluates the performance of the ARD models under in increasing
# level of sparsity
library(CompressiveNMF)
library(tidyverse)
library(foreach)
library(lsa)
library(doParallel)
library(sigminer)

source("R/Postprocess_functions.R")
source("R/simulate_data.R")
source("R/SignatureAnalyzer.R")

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

# Parameters of the simulation
ndatasets <- 20
pi0 <- c(0, 0.25, 0.5, 0.75, 0.9)
J <- 100
K_new <- 2

# Simulate the array of data
regenerate_data <- FALSE

main_dir <- "~/CompressiveNMF/output/sparsity_simulation/"
create_directory(main_dir)

if(regenerate_data){
  # Generate the data with increasing sparsity level
  for(j in 1:length(pi0)){
    data_all <- lapply(1:ndatasets, function(i) 
      simulate_data_sparse(J = J, 
                    K_new = K_new, 
                    overdispersion = 0, pi0 = pi0[j]))
    directory <- paste0(main_dir, "Sparsity_", pi0[j], "/")
    create_directory(directory)
    # Save data
    saveRDS(data_all, file = paste0(directory, "data.rds.gzip"), compress = "gzip")
  }
}

#######################################################
# Part 1 - Run the simulation
#######################################################
run_compNMF <- FALSE
run_compNMF_cos <- TRUE
run_ARD <- FALSE

# Number of samples
nsamples <- 1000
burnin <- 3000

registerDoParallel(min(c(ndatasets, parallel::detectCores() - 1)))
#---------------------------------- 1 - CompressiveNMF 
if (run_compNMF) {
  for(j in 1:length(pi0)){
    print(pi0[j])
    # Load the data
    out_dir <-  paste0(main_dir, "Sparsity_", pi0[j], "/")
    data_all <- readRDS(paste0(out_dir, "data.rds.gzip"))
    # Run the model in parallell
    out_CompressiveNMF <- foreach(i = 1:length(data_all)) %dopar% {
      out <- CompressiveNMF(X = data_all[[i]]$X, 
                            K = 20, 
                            epsilon = 0.01, 
                            nsamples = nsamples, 
                            burnin = burnin, 
                            nchains = 1, 
                            ncores = 1)
      res <- Postprocess_Compressive(out, data = data_all[[i]])
      #saveRDS(out, file = paste0(out_name, "results_", i, ".rds.gzip"), compress = "gzip")
      res
    }
    saveRDS(out_CompressiveNMF, file = paste0(out_dir, "CompressiveNMF.rds.gzip"), compress = "gzip")
  }
}

#---------------------------------- 2 - CompressiveNMF with COSMIC data and de novo signatures
if (run_compNMF_cos) {
  for(j in 1:length(pi0)){
    print(pi0[j])
    # Load the data
    out_dir <-  paste0(main_dir, "Sparsity_", pi0[j], "/")
    data_all <- readRDS(paste0(out_dir, "data.rds.gzip"))
    # Run the model in parallell
    out_CompressiveNMF_cos <- foreach(i = 1:length(data_all)) %dopar% {
      out <- CompressiveNMF(X = data_all[[i]]$X, 
                            K = 10,
                            use_cosmic = TRUE,
                            epsilon = 0.01, 
                            nsamples = nsamples, 
                            burnin = burnin, 
                            nchains = 1, 
                            ncores = 1)
      res <- Postprocess_Compressive(out, data = data_all[[i]])
      res
    }
    saveRDS(out_CompressiveNMF_cos, file = paste0(out_dir, "CompressiveNMF_cos.rds.gzip"), compress = "gzip")
  }
}

#---------------------------------- 3 - ARD with SignatureAnalyzer
if (run_ARD) {
  for(j in 1:length(pi0)){
    print(pi0[j])
    # Load the data
    out_dir <-  paste0(main_dir, "Sparsity_", pi0[j], "/")
    data_all <- readRDS(paste0(out_dir, "data.rds.gzip"))
    # Run the model in parallell
    out_ARD <- lapply(1:length(data_all), function(i) {
      print(i)
      out <- SignatureAnalyzer(X = data_all[[i]]$X, Kmax = 25)
      res <- Postprocess_ARD(out, data = data_all[[i]])
      res
    })
    saveRDS(out_ARD, file = paste0(out_dir, "ARD.rds.gzip"), compress = "gzip")
  }
}

#---------------------------------- 4 - ARD with BayesNMF (use Python script)
# Consider the python code


#######################################################
# Part 2 - Make plots
#######################################################


# out <- SignatureAnalyzer(X = data_all[[1]]$X, Kmax = 25)
# res <- Postprocess_ARD(out, data = data_all[[10]])
# 
# unlist(lapply(out_CompressiveNMF_cos, function(x) x$results[1]))
# unlist(lapply(out_CompressiveNMF, function(x) x$results[1]))
# 
# plot(out$Signature.norm %*% out$Exposure, data_all[[1]]$X)
# 
# i <- 1
# plot(out_ARD[[i]]$Lambda, data_all[[i]]$X)
# out_ARD[[1]]$results
# sqrt(mean((out_ARD[[i]]$Lambda- data_all[[i]]$X)^2))






