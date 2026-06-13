# This file runs a simulation that compares the results of CompressiveNMF vs
# the methods that allow for an infinite number of factors. We consider five 
# alternatives
#
# --- 0) CompressiveNMF
# --- 1) PoissonCUSP with gamma spike and slab
# --- 2) PoissonCUSP with gamma spike and slab, as in Kowal and Canale (2023)
# --- 3) PoissonCUSP with inverse gamma spike and slab
# --- 4) PoissonCUSP with inverse gamma spike and slab, as in Kowal and Canale (2023)
# --- 5) Multiplicative Gamma Process, Bhattacharya and Dunson (2011)
#
# We use the same data as in the main simulation. 
library(patchwork)
library(tidyverse)
library(CompressiveNMF)
library(sigminer)
library(foreach)
library(doParallel)
library(LaplacesDemon)
# Source useful functions
source("~/CompressiveNMF/R/Postprocess_functions.R")
source("~/CompressiveNMF/R/Poisson_InfiniteFactors.R")

# Useful functions
#---------------------------------------------------------------------------
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


simulate_data <- function(J = 100, cosmic_sig = c("SBS1", "SBS2", "SBS5", "SBS13"), 
                          K_new = 0, alpha = 0.25, theta = 100, overdispersion = 0){
  
  # Generate the signatures
  load("data/Cosmic_data_no_artifacts.rdata")
  Rmat_random <- t(rdirichlet(n = K_new, alpha = rep(alpha, 96)))
  if(K_new > 0) colnames(Rmat_random) <- paste0("SBSnew", 1:K_new)
  Rmat_cos <- as.matrix(cosmic_data[, cosmic_sig])
  Rmat <- as.matrix(cbind(Rmat_cos, Rmat_random))
  rownames(Rmat) <- cosmic_data$Channel
  
  # Generate the weights
  K <- ncol(Rmat)
  exposures <- rgamma(K, theta, 1)
  Theta <- matrix(rgamma(K * J, 0.5, 0.5), ncol = J, nrow = K)
  Theta <- apply(Theta, 2, function(x) x * exposures)
  rownames(Theta) <- colnames(Rmat)
  
  # Generate the counts
  Lambda <- Rmat %*% Theta
  X <- matrix(rnbinom(length(Lambda), size = 1/overdispersion, mu = c(Lambda)), nrow = 96, ncol = J)
  rownames(X) <- rownames(Rmat) 
  return(list(X = X, Rmat = Rmat, Theta = Theta))
}

#---------------------------------------------------------------------------

run_models_Comp_vs_Infinite <- function(J, overd, K_new, theta,
                                        simulation_dir,
                                        K_start, epsilon, mu_inf, 
                                        nsamples, burnin, 
                                        nsims, ncores, 
                                        cosmic_sig,
                                        generate_data = TRUE,
                                        runCompressive = TRUE, 
                                        runCUSP_gamma = TRUE, 
                                        runCUSP_gammaKC = TRUE,
                                        runCUSP_ig = TRUE, 
                                        runCUSP_ig_KC = TRUE, 
                                        runMGP = TRUE) {
  
  registerDoParallel(ncores)
  
  # Create a directory
  out_dir <- paste0(simulation_dir, "/Scenario_", J, "_overdisp_", overd, "_Knew_", K_new, "_theta_", theta)
  create_directory(out_dir)
  name_run <- paste(J, overd, K_new, theta, sep = " ")
  
  # Simulate the data
  set.seed(10)
  if(generate_data){
    data_all <- lapply(1:nsims, function(i) {
      data <- simulate_data(J = J, cosmic_sig = cosmic_sig, 
                            K_new = K_new, alpha = 0.25, theta = theta, overdispersion = overd)
    })
    saveRDS(data_all, file = paste0(out_dir, "/data.rds.gzip"), compress = "gzip")  
  } else {
    data_all <- readRDS(paste0(out_dir, "/data.rds.gzip"))
  }
  
  
  #-------- Model 0 - CompressiveNMF
  set.seed(10, kind = "L'Ecuyer-CMRG")
  if(runCompressive) {
    out_CompressiveNMF <- foreach(i = 1:nsims) %dopar% {
      out_name <- paste0(out_dir, "/CompressiveNMF")
      create_directory(out_name)
      out <- CompressiveNMF(X = data_all[[i]]$X, K = K_start, 
                            epsilon = epsilon, 
                            nsamples = nsamples, burnin = burnin, 
                            nchains = 1, ncores = 1, verbose = FALSE)
      res <- Postprocess_Compressive(out, data = data_all[[i]])
      saveRDS(out, file = paste0(out_name, "/results_", i, ".rds.gzip"), compress = "gzip")
      res
    }
    saveRDS(out_CompressiveNMF, file = paste0(out_dir, "/CompressiveNMF.rds.gzip"), compress = "gzip")
  }
  
  
  #-------- Model 1 - Poisson CUSP, gamma
  if (runCUSP_gamma) {
    set.seed(10, kind = "L'Ecuyer-CMRG")
    print(paste0(name_run, "- Run CUSP gamma"))
    out_PoissonCUSP <- foreach(i = 1:nsims) %dopar% {
      out_name <- paste0(out_dir, "/PoissonCUSP_gamma")
      create_directory(out_name)
      out <- PoissonCUSP(
        X = data_all[[i]]$X, K = K_start, a0 = 1, b0 = 1,
        KowalCanale = FALSE, random_alpha_sp = FALSE,
        nsamples = nsamples, burnin = burnin, 
        alpha = 0.5, mu_inf = mu_inf, 
        alpha_sp = 5)
      res <- Postprocess_PoissonCUSP_v2(out, data = data_all[[i]])
      saveRDS(out, file = paste0(out_name, "/results_", i, ".rds.gzip"), compress = "gzip")
      res
    }
    saveRDS(out_PoissonCUSP, file = paste0(out_dir, "/PoissonCUSP_gamma.rds.gzip"), compress = "gzip")
  }
  
  #-------- Model 2 - Poisson CUSP, gamma, Kolwal Canale
  if (runCUSP_gammaKC) {
    set.seed(10, kind = "L'Ecuyer-CMRG")
    print(paste0(name_run, "- Run CUSP gamma Kolwal Canale"))
    out_PoissonCUSP_KC <- foreach(i = 1:nsims) %dopar% {
      out_name <- paste0(out_dir, "/PoissonCUSP_gammaKC")
      create_directory(out_name)
      out <- PoissonCUSP(
        X = data_all[[i]]$X, K = K_start, a0 = 1, b0 = 1, 
        KowalCanale = TRUE, random_alpha_sp = TRUE,
        nsamples = nsamples, burnin = burnin, 
        alpha = 0.5, mu_inf = mu_inf, 
        alpha_sp = 5)
      res <- Postprocess_PoissonCUSP_v2(out, data = data_all[[i]])
      saveRDS(out, file = paste0(out_name, "/results_", i, ".rds.gzip"), compress = "gzip")
      res
    }
    saveRDS(out_PoissonCUSP_KC, file = paste0(out_dir, "/PoissonCUSP_gammaKC.rds.gzip"), compress = "gzip")
  }
  
  #-------- Model 3 - Poisson CUSP, inverse gamma
  if (runCUSP_ig) {
    set.seed(10, kind = "L'Ecuyer-CMRG")
    print(paste0(name_run, "- Run CUSP invgamma"))
    out_PoissonCUSP_ig <- foreach(i = 1:nsims) %dopar% {
      out_name <- paste0(out_dir, "/PoissonCUSP_invgamma")
      create_directory(out_name)
      out <- PoissonCUSP_ig(
        X = data_all[[i]]$X, K = K_start, a0 = 2, b0 = 1, 
        KowalCanale = FALSE, random_alpha_sp = FALSE,
        nsamples = nsamples, burnin = burnin, 
        alpha = 0.5, mu_inf = mu_inf, 
        alpha_sp = 5)
      res <- Postprocess_PoissonCUSP_v2(out, data = data_all[[i]])
      saveRDS(out, file = paste0(out_name, "/results_", i, ".rds.gzip"), compress = "gzip")
      res
    }
    saveRDS(out_PoissonCUSP_ig, file = paste0(out_dir, "/PoissonCUSP_invgamma.rds.gzip"), compress = "gzip")
  }
  
  #-------- Model 4 - Poisson CUSP, inverse gamma, Kolwal Canale
  if (runCUSP_ig_KC) {
    set.seed(10, kind = "L'Ecuyer-CMRG")
    print(paste0(name_run, "- Run CUSP invgamma Kolwal Canale"))
    out_PoissonCUSP_igKC <- foreach(i = 1:nsims) %dopar% {
      out_name <- paste0(out_dir, "/PoissonCUSP_invgammaKC")
      create_directory(out_name)
      out <- PoissonCUSP_ig(
        X = data_all[[i]]$X, K = K_start, a0 = 2, b0 = 1, 
        KowalCanale = TRUE, random_alpha_sp = TRUE,
        nsamples = nsamples, burnin = burnin, 
        alpha = 0.5, mu_inf = mu_inf, 
        alpha_sp = 5)
      res <- Postprocess_PoissonCUSP_v2(out, data = data_all[[i]])
      saveRDS(out, file = paste0(out_name, "/results_", i, ".rds.gzip"), compress = "gzip")
      res
    }
    saveRDS(out_PoissonCUSP_igKC, file = paste0(out_dir, "/PoissonCUSP_invgammaKC.rds.gzip"), compress = "gzip")
  }
  
  #-------- Model 5 - MGP
  if (runMGP) {
    set.seed(10, kind = "L'Ecuyer-CMRG")
    print(paste0(name_run, "- Multiplicative Gamma Process"))
    out_MGP <- foreach(i = 1:nsims) %dopar% {
      out_name <- paste0(out_dir, "/MGP")
      create_directory(out_name)
      #out <- PoissonMGP(X = data_all[[i]]$X, K = K_start, 
      #                  c0 = 2.1, d0 = 3.1, nsamples = nsamples, burnin = burnin)
      out <- readRDS(paste0(out_name, "/results_", i, ".rds.gzip"))
      res <- Postprocess_PoissonMGP(out, data = data_all[[i]])
      #saveRDS(out, file = paste0(out_name, "/results_", i, ".rds.gzip"), compress = "gzip")
      res
    }
    saveRDS(out_MGP, file = paste0(out_dir, "/PoissonMGP.rds.gzip"), compress = "gzip")
  }
}

#--- Run the simulation, iterating across all scenarios

# Now, run the simulation
theta_list <- c(100)
overdispersion_list <- c(0, 0.15)
J_list <- c(50, 100, 200)
K_new_list <- c(2, 6)
simulation_dir <- "~/CompressiveNMF/output/Compressive_vs_InfiniteFactors/"
create_directory(simulation_dir)

nsamples <- 1000
burnin <- 4000
nsims <- ncores <- 20
K_start <- 20
rerun <- FALSE # <----- Set to true to re-run. 
if(rerun){
  cat("start \n", file = paste0(simulation_dir, "check.txt"))
  # Run the simulation
  for(theta in theta_list){
    for(K_new in K_new_list){
      for(overd in overdispersion_list){
        for(J in J_list){
          try(run_models_Comp_vs_Infinite(J = J, 
                                          K_new = K_new, 
                                          theta = theta,
                                          overd = overd,
                                          simulation_dir = simulation_dir, 
                                          K_start = K_start, 
                                          cosmic_sig = c("SBS1", "SBS2", "SBS5", "SBS13"),
                                          epsilon = 0.01, 
                                          mu_inf = 0.01, 
                                          nsamples = nsamples, 
                                          burnin = burnin, 
                                          nsims = nsims, 
                                          ncores = ncores, 
                                          generate_data = TRUE,
                                          runCompressive = TRUE, 
                                          runCUSP_gamma = TRUE, 
                                          runCUSP_gammaKC = TRUE,
                                          runCUSP_ig = TRUE, 
                                          runCUSP_ig_KC = TRUE, 
                                          runMGP = TRUE),
              silent = FALSE, outFile = "log.txt")
        }
      }
    }
  }
  # Save now all results.
  df_all <- data.frame()
  for(theta in theta_list){
    for(K_new in K_new_list){
      for(overd in overdispersion_list){
        for(J in J_list){
          print(c(theta, K_new, overd, J))
          out_dir <- paste0(simulation_dir, "/Scenario_", J, "_overdisp_", overd, "_Knew_", K_new, "_theta_", theta)
          # Open results files
          files <- list.files(out_dir)
          files <- files[grepl("rds.gzip", files) & !grepl("data", files)]
          df_tmp <- data.frame()
          for(f in files){
            out <- open_rds_file(paste0(out_dir, "/", f))
            res <- extract_results(out, name = gsub(".rds.gzip", "", f))
            df_tmp <- rbind(df_tmp, res)
          }
          df_tmp$J = J
          df_tmp$theta = theta
          df_tmp$K_new = K_new
          df_tmp$overd = overd 
          df_all <- rbind(df_all, df_tmp)
        }
      }
    }
  }
  write_tsv(df_all, file = "~/CompressiveNMF/output/Compressive_vs_InfiniteFactors/Results_CompNMF_InfiniteFactors.tsv")
  
}

df_all <- read_tsv("~/CompressiveNMF/output/Compressive_vs_InfiniteFactors/Results_CompNMF_InfiniteFactors.tsv")

df_all$Kds <- paste0("K = ", df_all$K_new + 4)
df_all$Kds <- factor(df_all$Kds, levels = paste0("K = ", c(6, 10)))
df_all$Ktrue <- df_all$K_new + 4
df_all$Nds <-  paste0("N = ", df_all$J)
df_all$Nds <-  factor(df_all$Nds, levels = paste0("N = ", c(50, 100, 200)))
df_all$tau <-  paste0("tau = ", df_all$overd)

df_all$Method2 <- case_when(df_all$Method == "CompressiveNMF" ~ "CompNMF", 
                            TRUE ~ gsub("_", "-", gsub("Poisson", "", df_all$Method)))
 
# Estimated K
plotK <- ggplot(df_all, aes(x = Method2, y = K, fill = Method2, color = Method2)) +
  geom_hline(aes(yintercept = Ktrue), linetype = "dashed", color = "gray30") +
  geom_boxplot(alpha = 0.5) +
  facet_grid(tau ~ Kds + Nds)+
  geom_jitter(width = 0.05, height = 0.0) +
  theme_bw()+
  theme(legend.position = "none") +
  ylab("Estimated K") +
  theme(axis.title.x = element_blank(), 
        axis.text.x  = element_text(angle = 50, hjust = 1, size = 8.5))

# F1 score
plotF1 <- ggplot(df_all, aes(x = Method2, y = 2*(Precision * Sensitivity)/(Precision + Sensitivity), 
                   fill = Method2, color = Method2)) +
  geom_boxplot(alpha = 0.5) +
  facet_grid(tau ~ Kds + Nds)+
  geom_jitter(width = 0.05, height = 0.0) +
  theme_bw()+
  theme(legend.position = "none") +
  ylab("F1 score") +
  theme(axis.title.x = element_blank(), 
        axis.text.x  = element_text(angle = 50, hjust = 1, size = 8.5))
  

plotK / plotF1
ggsave(filename = "~/CompressiveNMF/figures/Compressive_vs_InfiniteFactors.pdf", 
       width = 9.72, height = 7.91)













