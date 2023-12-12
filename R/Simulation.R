# This file runs the simulation in the paper
library(tidyverse)
library(SparseSignatures)
library(signeR)
library(coda)
library(parallel)
library(lsa) # for cosine similarity

# Source the functions
source("R/CompressiveNMF.R")
source("R/NMF_l1_ARD.R")
source("R/PoissonCUSP.R")

# Function to simulate the data
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

data <- generate_data(I = 96, J = 100, K0 = 5)

X <- data$X

# We are making a comparison between the following methods
# 1) CompressiveNMF (our method)
# 2) Cumulative Shrinkage prior (CUSP) in Legramanti et. al (2020). We reforumate the model so that it selects the spike
# 3) Tan and Fevotte (2013) ARD method - already used in SignatureAnalyzer
# 4) Lahl et al. (2021) - SparseSignatures, uses the LASSO and derives the number of signatures via CV.
# 5) Rosales et al (2013) - signeR, which is a Bayesian method.

run_models <- function(J, K0, K, nsamples, burnin){
  out <- NULL
  # Generate data
  data <- generate_data(J = J, K0 = K0)
  out$data <- data
  # Model 1 - Compressive hyperprior
  out$CompressiveNMF <- CompressiveNMF(data$X, K = K, nsamples = nsamples, burnin = burnin)
  # Model 2 - SignatureAnalyzer 
  out$NMF_l1 <- NMF_l1_ARD(V = data$X, K = K)
  # Model 3 - CUSP
  out$PoissonCUSP <- PoissonCUSP(X = data$X, K = K, nsamples = nsamples,
                                 burnin = burnin, a0 = 0.5, b0 = 0.5, alpha_sp = K0)
  return(out)
}

# Run the simulation
nsamples <- 3000
burnin <- 1000
nsims <- 20

#------------------------------------------------------------ 100 data points 
J <- 100
K0 <- 5
K <- 20
set.seed(10, kind = "L'Ecuyer-CMRG")
out100 <- mclapply(X = c(1:nsims), 
                FUN = function(i) run_models(J, K0, K, nsamples, burnin), 
                affinity.list = 1:nsims,
                mc.preschedule = FALSE)
# Save the output
saveRDS(out100, file = paste0("output/Simulation_", J, ".rds.gzip"), compress = "gzip")

#------------------------------------------------------------ 200 data points 
J <- 250
K0 <- 5
K <- 20
set.seed(10, kind = "L'Ecuyer-CMRG")
out200 <- mclapply(X = c(1:nsims), 
                   FUN = function(i) run_models(J, K0, K, nsamples, burnin), 
                   affinity.list = 1:nsims,
                   mc.preschedule = FALSE)
# Save the output
saveRDS(out200, file = paste0("output/Simulation_", J, ".rds.gzip"), compress = "gzip")

#------------------------------------------------------------ 400 data points 
J <- 500
K0 <- 10
K <- 20
set.seed(10, kind = "L'Ecuyer-CMRG")
out500 <- mclapply(X = c(1:nsims), 
                   FUN = function(i) run_models(J, K0, K, nsamples, burnin), 
                   affinity.list = 1:nsims,
                   mc.preschedule = FALSE)
# Save the output
saveRDS(out500, file = paste0("output/Simulation_", J, ".rds.gzip"), compress = "gzip")

########## Postprocess and make a table to evaluate the performance
# We will check: RMSE, cosine similarity between the inferred signatures, number of inferred signatures

# Postprocess CompressiveNMF
tab_comp_100 <- do.call("rbind", lapply(out100, function(x) Postprocess_Compressive(x$CompressiveNMF, data = x$data)))
df_comp_100 <- data.frame(tab_comp_100, "method" = "CompressiveNMF", J = 100)
tab_comp_200 <- do.call("rbind", lapply(out200, function(x) Postprocess_Compressive(x$CompressiveNMF, data = x$data)))
df_comp_200 <- data.frame(tab_comp_200, "method" = "CompressiveNMF", J = 200)

# Postprocess SignatureAnalyzer
tab_ard_100 <- do.call("rbind", lapply(out100, function(x) Postprocess_ARD(x$NMF_l1, data = x$data)))
df_ard_100 <- data.frame(tab_ard_100, "method" = "ARD", J = 100)
tab_ard_200 <- do.call("rbind", lapply(out200, function(x) Postprocess_ARD(x$NMF_l1, data = x$data)))
df_ard_200 <- data.frame(tab_ard_200, "method" = "ARD", J = 200)

# Postprocess PoissonCUSP
tab_cusp_100 <- do.call("rbind", lapply(out100, function(x) Postprocess_PoissonCUSP(x$PoissonCUSP, data = x$data)))
df_cusp_100 <- data.frame(tab_cusp_100, "method" = "CUSP", J = 100)
tab_cusp_200 <- do.call("rbind", lapply(out200, function(x) Postprocess_PoissonCUSP(x$PoissonCUSP, data = x$data)))
df_cusp_200 <- data.frame(tab_cusp_200, "method" = "CUSP", J = 200)



results <- rbind(df_comp_100, df_comp_200, df_ard_100, df_ard_200)
results$J <- as.factor(results$J)
ggplot(results) +
  geom_boxplot(aes(x = J, y = cos_sim, col = method))

ggplot(results) +
  geom_boxplot(aes(x = J, y = rmse_Lambda, col = method))


results %>%
  group_by(J, method) %>%
  summarize_all(mean) 



library(SparseSignatures)

library("BSgenome.Hsapiens.1000genomes.hs37d5")
bsg = BSgenome.Hsapiens.1000genomes.hs37d5
data(ssm560_reduced)
head(ssm560_reduced)
data(mutation_categories)
head(mutation_categories)

imported_data = import.trinucleotides.counts(data=ssm560_reduced,reference=bsg)
head(imported_data)

patients.plot(trinucleotides_counts=imported_data,samples="PD10010a")

data(patients)
head(patients)

starting_betas = startingBetaEstimation(x=patients,K = 3:12, background_signature = NULL)
dim(starting_betas)

dim(starting_betas[[5]])

cv = nmfLassoCV(x=patients,K=3:6)








