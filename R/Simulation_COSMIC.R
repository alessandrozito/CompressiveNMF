# This file implements a simulation using the signatures inferred in the COSMIC database.
# We use version of October 2023
library(tidyverse)
library(SparseSignatures)
library(signeR) # Note that this package constains v 3.2 of the cosmic dataset. We will be usign v3.4 (october 2023)
library(coda)
library(parallel)
library(lsa) # for cosine similarity
library(rvest) # to scrape the data

# Source the functions
source("R/CompressiveNMF.R")
source("R/NMF_l1_ARD.R")
source("R/PoissonCUSP.R")
source("R/plot_signatures.R")
source("R/Postprocess_functions.R")

################################################################################
# Simulation functions
################################################################################
simulate_mutations_from_cosmic <- function(K0 = 5, J = 100, signal = 100) {
  # Load the cosmic dataset containing the signautures
  load("data/Cosmic_data.rdata")
  I <- nrow(cosmic_data)
  # Randomly sample K0 signature without replacement
  signature_names <- colnames(cosmic_data)[-c(1, 2, 3)]
  sampled_sign <- sample(signature_names, size = K0)
  Rmat <- as.matrix(cosmic_data[, sampled_sign])
  rownames(Rmat) <- cosmic_data$Channel
  # Randomly sample the weight matrix
  Theta <- signal * matrix(rgamma(n = K0 * J, 1, 0.5), nrow = K0, ncol = J)
  # Generate the counts
  Lambda <- Rmat %*% Theta
  X <- matrix(rpois(n = length(Lambda), c(Lambda)), nrow = I, ncol = J)
  return(list(X = X, Rmat = Rmat, Theta = Theta))
}


simulate_data <- function(K_new = NULL, K_cosmic = NULL, J = 100, signal = 100,
                          perturb_cosmic = TRUE, alpha_new = 0.1, sd_perturb = 0.005) {
  I <- 96
  Rmat <- NULL
  # Add a perturbed version of the cosmic signatures
  if (!is.null(K_cosmic)) {
    load("data/Cosmic_data.rdata")
    signature_names <- colnames(cosmic_data)[-c(1, 2, 3)]
    sampled_sign <- sample(signature_names, size = K_cosmic)
    Rcosmic <- as.matrix(cosmic_data[, sampled_sign])
    if (perturb_cosmic) {
      perturbations <- matrix(rnorm(I * K_cosmic,mean = 0, sd = sd_perturb), nrow = I, ncol = K_cosmic)
      Rtemp <- Rcosmic + perturbations
      Rtemp[Rtemp < 0] <- 0
      Rcosmic <- apply(Rtemp, 2, function(x) x / sum(x))
    }
    Rmat <- cbind(Rmat, Rcosmic)
  }
  # Add additional extra signatures
  if (!is.null(K_new)) {
    Rnew <- matrix(rgamma(n = I * K_new, alpha_new, 0.5), nrow = I, ncol = K_new)
    Rnew <- apply(Rnew, 2, function(x) x / sum(x))
    colnames(Rnew) <- paste0("SBSnew", 1:K_new)
    Rmat <- cbind(Rmat, Rnew)
  }
  # Sample the weights
  Theta <- signal * matrix(rgamma(n = ncol(Rmat) * J, 1, 0.5), nrow = ncol(Rmat), ncol = J)
  # Generate the counts
  Lambda <- Rmat %*% Theta
  X <- matrix(rpois(n = length(Lambda), c(Lambda)), nrow = I, ncol = J)
  return(list(X = X, Rmat = Rmat, Theta = Theta))
}


#### Functions to post-process the output
extract_results <- function(out, name){
  df <- data.frame("Method" = name, "Simulation" = 1:length(out))
  df <- cbind(df, as.data.frame(t(sapply(1:length(out), function(i) out[[i]]$results, simplify = TRUE))))
  return(df)
}


# results functions
aggregate_results <- function(out_CompressiveNMF = NULL,
                              out_CompressiveNMF_cosmic = NULL,
                              out_PoissonCUSP = NULL, 
                              out_ARD = NULL, out_signeR = NULL) {
  df_res <- NULL
  #---------- CompressiveNMF
  if(!is.null(out_CompressiveNMF)){
    df_res <- rbind(df_res, extract_results(out_CompressiveNMF, name = "CompressiveNMF")) 
  }
  #---------- CompressiveNMF with cosmic database 
  if(!is.null(out_CompressiveNMF_cosmic)){
    df_res <- rbind(df_res, extract_results(out_CompressiveNMF_cosmic, name = "CompressiveNMF_cosmic")) 
  }
  #---------- PoissonCUSP
  if(!is.null(out_PoissonCUSP)){
    df_res <- rbind(df_res, extract_results(out_PoissonCUSP, name = "PoissonCUSP")) 
  }
  #---------- SignatureAnalyzer
  if(!is.null(out_ARD)){
    df_res <- rbind(df_res, extract_results(out_ARD, name = "SignatureAnalyzer")) 
  }
  #---------- signeR
  if(!is.null(out_signeR)){
    df_res <- rbind(df_res, extract_results(out_signeR, name = "signeR")) 
  }
  return(df_res)
}

# Set the number of cores to use
ncores <- 20
################################################################################
# Simulation 1 
################################################################################
# We let K_cosmic = 5, K_new = 3, J = 50, signal = 50
# Set the number of replicates
nsims <- 40
nsamples <- 3000
burnin <- 2000

# Simulate the data and save
set.seed(10)
data_all <- lapply(1:nsims, function(i) simulate_data(K_new = 3, K_cosmic = 5, J = 50, signal = 50))
saveRDS(data_all, file = "output/Simulation1/data.rds.gzip", compress = "gzip")

#---------------------------------- 1 - CompressiveNMF with no prior information
K <- 20
alpha <- 0.4
set.seed(10, kind = "L'Ecuyer-CMRG")
out_CompressiveNMF <- mclapply(X = c(1:nsims),
                               FUN = function(i) { 
                                 out <- CompressiveNMF(X = data_all[[i]]$X, K = K, 
                                                       nsamples = nsamples, burnin = burnin, alpha = alpha)
                                 Postprocess_Compressive(out, data = data_all[[i]])
                                 },
                                mc.cores = ncores,
                               mc.preschedule = FALSE)
saveRDS(out_CompressiveNMF, file = "output/Simulation1/CompressiveNMF.rds.gzip", compress = "gzip")


#---------------------------------- 2 - CompressiveNMF with whole COSMIC data
load("data/Cosmic_data.rdata")
K <- 10
S <- as.matrix(cosmic_data[, -c(1,2,3)])
alpha <- 0.4
set.seed(10, kind = "L'Ecuyer-CMRG")
out_CompressiveNMF_cosmic <- mclapply(X = c(1:nsims),
                               FUN = function(i) { 
                                 out <- CompressiveNMF(X = data_all[[i]]$X, K = K, S = S, 
                                                       nsamples = nsamples, burnin = burnin, alpha = alpha)
                                 Postprocess_Compressive(out, data = data_all[[i]])
                               },
                                mc.cores = ncores,
                               mc.preschedule = FALSE)
saveRDS(out_CompressiveNMF_cosmic, file = "output/Simulation1/CompressiveNMF_cosmic.rds.gzip", compress = "gzip")

#---------------------------------- 3 - PoissonCUSP
K <- 20
alpha <- 0.4
set.seed(10, kind = "L'Ecuyer-CMRG")
out_PoissonCUSP <- mclapply(X = c(1:nsims),
                                      FUN = function(i) { 
                                        out <- PoissonCUSP(X = data_all[[i]]$X, K = K, a0 = 2, b0 = 2, 
                                                              nsamples = nsamples, burnin = burnin, alpha = alpha)
                                        Postprocess_PoissonCUSP(out, data = data_all[[i]])
                                      },
                                       mc.cores = ncores,
                                      mc.preschedule = FALSE)
saveRDS(out_PoissonCUSP, file = "output/Simulation1/PoissonCUSP.rds.gzip", compress = "gzip")


#---------------------------------- 4 - ARD with SignatureAnalyzer
K <- 20
set.seed(10, kind = "L'Ecuyer-CMRG")
out_ARD <- mclapply(X = c(1:nsims),
                                      FUN = function(i) { 
                                        out <- NMF_l1_ARD(V = data_all[[i]]$X, a = 5, K = K, normalize = TRUE)
                                        Postprocess_ARD(out, data = data_all[[i]])
                                      },
                                       mc.cores = ncores,
                                      mc.preschedule = FALSE)
saveRDS(out_ARD, file = "output/Simulation1/ARD.rds.gzip", compress = "gzip")


#---------------------------------- 5 - signeR
Kmin <- 3
Kmax <- 20
set.seed(10, kind = "L'Ecuyer-CMRG")
out_signeR <- mclapply(X = c(1:nsims),
                       FUN = function(i) { 
                                        out <- signeR(M = data_all[[i]]$X, samples = "col", nlim = c(Kmin, Kmax))
                                        Postprocess_signeR(out, data = data_all[[i]])
                                      },
                                       mc.cores = ncores,
                                      mc.preschedule = FALSE)
saveRDS(out_signeR, file = "output/Simulation1/signeR.rds.gzip", compress = "gzip")


################################################################################
# Simulation 2
################################################################################
# We let K_cosmic = 5, K_new = 5, J = 100, signal = 50
# Set the number of replicates
nsims <- 40
nsamples <- 3000
burnin <- 2000
out_dir <- "output/Simulation2/"

# Simulate the data and save
set.seed(10)
data_all <- lapply(1:nsims, function(i) simulate_data(K_new = 5, K_cosmic = 5, J = 100, signal = 50))
saveRDS(data_all, file = paste0(out_dir, "data.rds.gzip"), compress = "gzip")

#---------------------------------- 1 - CompressiveNMF with no prior information
K <- 20
alpha <- 0.4
set.seed(10, kind = "L'Ecuyer-CMRG")
out_CompressiveNMF <- mclapply(X = c(1:nsims),
                               FUN = function(i) { 
                                 out <- CompressiveNMF(X = data_all[[i]]$X, K = K, 
                                                       nsamples = nsamples, burnin = burnin, alpha = alpha)
                                 Postprocess_Compressive(out, data = data_all[[i]])
                               },
                                mc.cores = ncores,
                               mc.preschedule = FALSE)
saveRDS(out_CompressiveNMF, file = paste0(out_dir, "CompressiveNMF.rds.gzip"), compress = "gzip")


#---------------------------------- 2 - CompressiveNMF with whole COSMIC data
load("data/Cosmic_data.rdata")
K <- 10
S <- as.matrix(cosmic_data[, -c(1,2,3)])
alpha <- 0.4
set.seed(10, kind = "L'Ecuyer-CMRG")
out_CompressiveNMF_cosmic <- mclapply(X = c(1:nsims),
                                      FUN = function(i) { 
                                        out <- CompressiveNMF(X = data_all[[i]]$X, K = K, S = S, 
                                                              nsamples = nsamples, burnin = burnin, alpha = alpha)
                                        Postprocess_Compressive(out, data = data_all[[i]])
                                      },
                                       mc.cores = ncores,
                                      mc.preschedule = FALSE)
saveRDS(out_CompressiveNMF_cosmic, file = paste0(out_dir, "CompressiveNMF_cosmic.rds.gzip"), compress = "gzip")

#---------------------------------- 3 - PoissonCUSP
K <- 20
alpha <- 0.4
set.seed(10, kind = "L'Ecuyer-CMRG")
out_PoissonCUSP <- mclapply(X = c(1:nsims),
                            FUN = function(i) { 
                              out <- PoissonCUSP(X = data_all[[i]]$X, K = K, a0 = 2, b0 = 2, 
                                                 nsamples = nsamples, burnin = burnin, alpha = alpha)
                              Postprocess_PoissonCUSP(out, data = data_all[[i]])
                            },
                             mc.cores = ncores,
                            mc.preschedule = FALSE)
saveRDS(out_PoissonCUSP, file = paste0(out_dir, "PoissonCUSP.rds.gzip"), compress = "gzip")


#---------------------------------- 4 - ARD with SignatureAnalyzer
K <- 20
set.seed(10, kind = "L'Ecuyer-CMRG")
out_ARD <- mclapply(X = c(1:nsims),
                    FUN = function(i) { 
                      out <- NMF_l1_ARD(V = data_all[[i]]$X, a = 5, K = K, normalize = TRUE)
                      Postprocess_ARD(out, data = data_all[[i]])
                    },
                     mc.cores = ncores,
                    mc.preschedule = FALSE)
saveRDS(out_ARD, file = paste0(out_dir, "ARD.rds.gzip"), compress = "gzip")


#---------------------------------- 5 - signeR
Kmin <- 3
Kmax <- 20
set.seed(10, kind = "L'Ecuyer-CMRG")
out_signeR <- mclapply(X = c(1:nsims),
                       FUN = function(i) { 
                         out <- signeR(M = data_all[[i]]$X, samples = "col", nlim = c(Kmin, Kmax))
                         Postprocess_signeR(out, data = data_all[[i]])
                       },
                        mc.cores = ncores,
                       mc.preschedule = FALSE)
saveRDS(out_signeR, file = paste0(out_dir, "signeR.rds.gzip"), compress = "gzip")


################################################################################
# Simulation 3
################################################################################
# We let K_cosmic = 15, K_new = 5, J = 400, signal = 80
# Set the number of replicates
nsims <- 40
nsamples <- 3000
burnin <- 2000
out_dir <- "output/Simulation3/"

# Simulate the data and save
set.seed(10)
data_all <- lapply(1:nsims, function(i) simulate_data(K_new = 5, K_cosmic = 15, J = 400, signal = 80))
saveRDS(data_all, file = paste0(out_dir, "data.rds.gzip"), compress = "gzip")

#---------------------------------- 1 - CompressiveNMF with no prior information
K <- 25
alpha <- 0.4
set.seed(10, kind = "L'Ecuyer-CMRG")
out_CompressiveNMF <- mclapply(X = c(1:nsims),
                               FUN = function(i) { 
                                 out <- CompressiveNMF(X = data_all[[i]]$X, K = K, 
                                                       nsamples = nsamples, burnin = burnin, alpha = alpha)
                                 Postprocess_Compressive(out, data = data_all[[i]])
                               },
                                mc.cores = ncores,
                               mc.preschedule = FALSE)
saveRDS(out_CompressiveNMF, file = paste0(out_dir, "CompressiveNMF.rds.gzip"), compress = "gzip")


#---------------------------------- 2 - CompressiveNMF with whole COSMIC data
load("data/Cosmic_data.rdata")
K <- 10
S <- as.matrix(cosmic_data[, -c(1,2,3)])
alpha <- 0.4
set.seed(10, kind = "L'Ecuyer-CMRG")
out_CompressiveNMF_cosmic <- mclapply(X = c(1:nsims),
                                      FUN = function(i) { 
                                        out <- CompressiveNMF(X = data_all[[i]]$X, K = K, S = S, 
                                                              nsamples = nsamples, burnin = burnin, alpha = alpha)
                                        Postprocess_Compressive(out, data = data_all[[i]])
                                      },
                                       mc.cores = ncores,
                                      mc.preschedule = FALSE)
saveRDS(out_CompressiveNMF_cosmic, file = paste0(out_dir, "CompressiveNMF_cosmic.rds.gzip"), compress = "gzip")

#---------------------------------- 3 - PoissonCUSP
K <- 25
alpha <- 0.4
set.seed(10, kind = "L'Ecuyer-CMRG")
out_PoissonCUSP <- mclapply(X = c(1:nsims),
                            FUN = function(i) { 
                              out <- PoissonCUSP(X = data_all[[i]]$X, K = K, a0 = 2, b0 = 2, 
                                                 nsamples = nsamples, burnin = burnin, alpha = alpha)
                              Postprocess_PoissonCUSP(out, data = data_all[[i]])
                            },
                             mc.cores = ncores,
                            mc.preschedule = FALSE)
saveRDS(out_PoissonCUSP, file = paste0(out_dir, "PoissonCUSP.rds.gzip"), compress = "gzip")


#---------------------------------- 4 - ARD with SignatureAnalyzer
K <- 25
set.seed(10, kind = "L'Ecuyer-CMRG")
out_ARD <- mclapply(X = c(1:nsims),
                    FUN = function(i) { 
                      out <- NMF_l1_ARD(V = data_all[[i]]$X, a = 5, K = K, normalize = TRUE)
                      Postprocess_ARD(out, data = data_all[[i]])
                    },
                    mc.cores = ncores,
                    # mc.cores = ncores,
                    mc.preschedule = FALSE)
saveRDS(out_ARD, file = paste0(out_dir, "ARD.rds.gzip"), compress = "gzip")


#---------------------------------- 5 - signeR
Kmin <- 3
Kmax <- 25
set.seed(10, kind = "L'Ecuyer-CMRG")
out_signeR <- mclapply(X = c(1:nsims),
                       FUN = function(i) { 
                         out <- signeR(M = data_all[[i]]$X, samples = "col", nlim = c(Kmin, Kmax))
                         Postprocess_signeR(out, data = data_all[[i]])
                       },
                        mc.cores = ncores,
                       mc.preschedule = FALSE)
saveRDS(out_signeR, file = paste0(out_dir, "signeR.rds.gzip"), compress = "gzip")


