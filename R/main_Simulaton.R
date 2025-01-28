library(LaplacesDemon)
library(tidyverse)
library(parallel)
library(coda)

# Source the functions
source("R/SignatureAnalyzer.R")
source("R/SigProfilerExtractor.R")
source("R/signeR.R")
source("R/PoissonCUSP.R")
source("R/CompressiveNMF.R")
source("R/Postprocess_functions.R")
#source("R/plot_signatures.R")

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


# This function simulates the data from 
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


# This function runs all models
main <- function(J, K_new, theta, overdispersion, 
                 cosmic_sig = c("SBS1", "SBS2", "SBS5", "SBS13"), 
                 nsims = 20, ncores = 20, 
                 simulation_dir = "~/CompressiveNMF/output/main_simulation/", 
                 generate_data = TRUE,
                 run_compNMF =  TRUE,
                 run_compNMF_cos = TRUE,
                 run_ARD = TRUE,
                 run_CUSP = TRUE,
                 run_signeR_code = TRUE,
                 run_sigProfiler = TRUE){
  
  registerDoParallel(ncores)
  # Create a directory
  out_dir <- paste0(simulation_dir, "/Scenario_", J, "_overdisp_", overdispersion, "_Knew_", K_new, "_theta_", theta)
  create_directory(out_dir)
  name_run <- paste(J, overdispersion, K_new, theta, sep = " ")
  # Simulate the data
  set.seed(10)
  if(generate_data){
    data_all <- lapply(1:nsims, function(i) {
      data <- simulate_data(J = J, cosmic_sig = cosmic_sig, 
                            K_new = K_new, alpha = 0.25, theta = theta, overdispersion = overdispersion)
    })
    saveRDS(data_all, file = paste0(out_dir, "/data.rds.gzip"), compress = "gzip")  
  } else {
    data_all <- readRDS(paste0(out_dir, "/data.rds.gzip"))
  }
  
  
  #---------------------------------- 1 - de novo CompressiveNMF
  nsamples <- 1000
  burnin <- 4000
  
  if (run_compNMF) {
    set.seed(10, kind = "L'Ecuyer-CMRG")
    print(paste0(name_run, "- Run compressive"))
    out_CompressiveNMF <- foreach(i = 1:nsims) %dopar% {
      out_name <- paste0(out_dir, "/CompressiveNMF")
      create_directory(out_name)
      out <- CompressiveNMF(X = data_all[[i]]$X, K = 20, epsilon = 0.001, nsamples = nsamples, burnin = burnin, nchains = 1, ncores = 1, progressbar = FALSE, verbose = FALSE)
      res <- Postprocess_Compressive(out, data = data_all[[i]])
      saveRDS(out, file = paste0(out_name, "/results_", i, ".rds.gzip"), compress = "gzip")
      res
    }
    saveRDS(out_CompressiveNMF, file = paste0(out_dir, "/CompressiveNMF.rds.gzip"), compress = "gzip")
  }
  
  
  #---------------------------------- 2 - CompressiveNMF with COSMIC data and de novo signatures
  if (run_compNMF_cos) {
    set.seed(10, kind = "L'Ecuyer-CMRG")
    print(paste0(name_run, "- Run compressive cosmic"))
    out_CompressiveNMF_cosmic <- foreach(i = 1:nsims) %dopar% {
      out_name <- paste0(out_dir, "/CompressiveNMF_cosmic")
      create_directory(out_name)
      out <- CompressiveNMF(X = data_all[[i]]$X, K = 15, use_cosmic = TRUE, betah_optimal = TRUE, epsilon = 0.001, nsamples = nsamples, burnin = burnin, nchains = 1, ncores = 1, progressbar = FALSE, verbose = FALSE)
      res <- Postprocess_Compressive(out, data = data_all[[i]])
      saveRDS(out, file = paste0(out_name, "/results_", i, ".rds.gzip"), compress = "gzip")
      res
    }
    saveRDS(out_CompressiveNMF_cosmic, file = paste0(out_dir, "/CompressiveNMF_cosmic.rds.gzip"), compress = "gzip")
  }
  
  #---------------------------------- 3 - ARD with SignatureAnalyzer
  if (run_ARD) {
    set.seed(10, kind = "L'Ecuyer-CMRG")
    print(paste0(name_run, "- Run ard"))
    out_ARD <- mclapply(
      X = c(1:nsims),
      FUN = function(i) {
        print(i)
        out_name <- paste0(out_dir, "/ARD")
        create_directory(out_name)
        # out <- NMF_l1_ARD(V = data_all[[i]]$X, a = 5, K = 25, normalize = TRUE)
        out <- SignatureAnalyzer(X = data_all[[i]]$X, Kmax = 25)
        res <- Postprocess_ARD(out, data = data_all[[i]])
        saveRDS(out, file = paste0(out_name, "/results_", i, ".rds.gzip"), compress = "gzip")
        res
      },
      mc.cores = 1,
      mc.preschedule = FALSE
    )
    saveRDS(out_ARD, file = paste0(out_dir, "/ARD.rds.gzip"), compress = "gzip")
  }
  
  #---------------------------------- 4 - PoissonCUSP
  if (run_CUSP) {
    set.seed(10, kind = "L'Ecuyer-CMRG")
    print(paste0(name_run, "- Run CUSP"))
    out_PoissonCUSP <- foreach(i = 1:nsims) %dopar% {
      out_name <- paste0(out_dir, "/PoissonCUSP")
        create_directory(out_name)
        out <- PoissonCUSP(
          X = data_all[[i]]$X, K = 20, a0 = 1, b0 = 1,
          nsamples = nsamples, burnin = burnin, alpha = 0.5,
          mu_inf = 0.01)
        res <- Postprocess_PoissonCUSP(out, data = data_all[[i]])
        saveRDS(out, file = paste0(out_name, "/results_", i, ".rds.gzip"), compress = "gzip")
        res
      }
    saveRDS(out_PoissonCUSP, file = paste0(out_dir, "/PoissonCUSP.rds.gzip"), compress = "gzip")
  }
  
  #---------------------------------- 5 - signeR
  if (run_signeR_code) {
    print(paste0(name_run, "- Run signeR"))
    set.seed(10, kind = "L'Ecuyer-CMRG")
    out_signeR <- foreach(i = 1:nsims) %dopar% {
      out_name <- paste0(out_dir, "/signeR")
      create_directory(out_name)
      out <- run_signeR(data_all[[i]]$X, Kmin = 2, Kmax = 20, estimate_hyper = FALSE, sequential = TRUE)
      res <- Postprocess_signeR(out, data = data_all[[i]])
      saveRDS(out, file = paste0(out_name, "/results_", i, ".rds.gzip"), compress = "gzip")
      res
    }
    saveRDS(out_signeR, file = paste0(out_dir, "/signeR.rds.gzip"), compress = "gzip")
  }
  
  #---------------------------------- 6 - sigprofiler
  if (run_sigProfiler) {
    print(paste0(name_run, "- Run sigpro"))
    set.seed(10, kind = "L'Ecuyer-CMRG")
    out_sigPro <- foreach(i = 1:nsims) %dopar% {
      out_dir_sigpro <- paste0(out_dir, "/sigpro/temp", i)
      out_name <- paste0(out_dir, "/SigPro")
      create_directory(out_name)
      out <- sigprofiler(
        X = data_all[[i]]$X, Kmin = 2, Kmax = 20,
        out_dir_sigpro = out_dir_sigpro,
        cores = 1)
      res <- Postprocess_SigProfiler(out, data = data_all[[i]])
      saveRDS(res, file = paste0(out_name, "/results_", i, ".rds.gzip"), compress = "gzip")
      res
    }
    saveRDS(out_sigPro, file = paste0(out_dir, "/sigProfiler.rds.gzip"), compress = "gzip")
  }
  
}

# Now, run the simulation
theta_list <- c(100)
overdispersion_list <- c(0, 0.15)
J_list <- c(50, 100, 200)
K_new_list <- c(2, 6)
simulation_dir <- "~/CompressiveNMF/output/main_simulation/"
create_directory(simulation_dir)

rerun <- FALSE # <----- Set to true to re-run. 
if(rerun){
  cat("start \n", file = "check.txt")
  # Run the simulation
  for(theta in theta_list){
    for(K_new in K_new_list){
      for(overd in overdispersion_list){
        for(J in J_list){
          try(main(J = J, K_new = K_new, theta = theta,
                   overdispersion = overd,
                   simulation_dir = simulation_dir, 
                   generate_data = TRUE,
                   run_compNMF =  TRUE,
                   run_compNMF_cos = TRUE,
                   run_ARD = TRUE,
                   run_CUSP = TRUE,
                   run_signeR_code = TRUE,
                   run_sigProfiler = TRUE),
              silent = FALSE, outFile = "log.txt")
        }
      }
    }
  }
}


#-------------- Unpack the output and save it
df_all <- data.frame()
for(theta in theta_list){
  for(K_new in K_new_list){
    for(overd in overdispersion_list){
      for(J in J_list){
        print(c(theta, K_new, overd, J))
        out_dir <- paste0(simulation_dir, "/Scenario_", J, "_overdisp_", overd, "_Knew_", K_new, "_theta_", theta)
        # Open results files
        out_CompressiveNMF <- open_rds_file(paste0(out_dir, "/CompressiveNMF.rds.gzip"))
        out_CompressiveNMF_cos <- open_rds_file(paste0(out_dir, "/CompressiveNMF_cosmic.rds.gzip")) 
        out_ARD <- open_rds_file(paste0(out_dir, "/ARD.rds.gzip"))
        out_PoissonCUSP <- open_rds_file(paste0(out_dir, "/PoissonCUSP.rds.gzip"))
        out_signeR <- open_rds_file(paste0(out_dir, "/signeR.rds.gzip"))
        out_sigPro <- open_rds_file(paste0(out_dir, "/sigProfiler.rds.gzip"))
        out_BayesNMF <- open_rds_file(paste0(out_dir, "/BayesNMF_brouwer.rds.gzip"))
        strange_sigpro <- unlist(lapply(out_sigPro, is.null))
        if(sum(strange_sigpro) > 0){
          out_sigPro <- out_sigPro[-which(strange_sigpro)] 
        }
        # Aggregate
        df <- aggregate_results(out_CompressiveNMF = out_CompressiveNMF,
                                out_CompressiveNMF_cosmic = out_CompressiveNMF_cos, 
                                out_PoissonCUSP = out_PoissonCUSP, 
                                out_ARD = out_ARD, 
                                out_signeR = out_signeR, 
                                out_sigPro = out_sigPro, 
                                out_BayesNMF = out_BayesNMF)
        df$J = J
        df$theta = theta
        df$K_new = K_new
        df$overd = overd 
        df_all <- rbind(df, df_all)
      }
    }
  }
}

# Save the output
write_csv(x = df_all, file =  paste0(simulation_dir, "/simulation_output.csv"))

#-----------------------------------------------------------  Extract the F1 range
df_F1 <- data.frame()
for(theta in theta_list){
  for(K_new in K_new_list){
    for(overd in 0.15){
      for(J in J_list){
        df_temp <- data.frame()
        print(c(theta, K_new, overd, J))
        out_dir <- paste0(simulation_dir, "/Scenario_", J, "_overdisp_", overd, "_Knew_", K_new, "_theta_", theta)
        # CompressiveNMF
        print("CompNMF")
        out_CompressiveNMF <- open_rds_file(paste0(out_dir, "/CompressiveNMF.rds.gzip"))
        df_temp <- rbind(df_temp, extract_F1_range(out_CompressiveNMF, method = "2.CompNMF"))
        # CompressiveNMF + cosmic
        print("CompNMFcos")
        out_CompressiveNMF_cos <- open_rds_file(paste0(out_dir, "/CompressiveNMF_cosmic.rds.gzip")) 
        df_temp <- rbind(df_temp, extract_F1_range(out_CompressiveNMF_cos, method = "1.CompNMFcos"))
        # ARD
        print("ARD")
        out_ARD <- open_rds_file(paste0(out_dir, "/ARD.rds.gzip"))
        df_temp <- rbind(df_temp, extract_F1_range(out_ARD, method = "5.ARD"))
        # CUSP
        print("CUSP")
        out_PoissonCUSP <- open_rds_file(paste0(out_dir, "/PoissonCUSP.rds.gzip"))
        df_temp <- rbind(df_temp, extract_F1_range(out_PoissonCUSP, method = "6.CUSP"))
        # SigneR
        print("signeR")
        out_signeR <- open_rds_file(paste0(out_dir, "/signeR.rds.gzip"))
        df_temp <- rbind(df_temp, extract_F1_range(out_signeR, method = "3.signeR"))
        # sigPro
        print("sigPro")
        out_sigPro <- open_rds_file(paste0(out_dir, "/sigProfiler.rds.gzip"))
        df_temp <- rbind(df_temp, extract_F1_range(out_sigPro, method = "4.SigPro"))
        # BayesNMF
        print("BayesNMF")
        out_BayesNMF <- open_rds_file(paste0(out_dir, "/BayesNMF_brouwer.rds.gzip"))
        df_temp <- rbind(df_temp, extract_F1_range(out_BayesNMF, method = "7.BayesNMF"))
        # Merge everything
        df_temp$J = J
        df_temp$theta = theta
        df_temp$K_new = K_new
        df_temp$overd = overd 
        df_F1 <- rbind(df_F1, df_temp)
      }
    }
  }
}

#write_csv(x = df_F1, file = paste0(simulation_dir, "/df_F1.csv"))
write_csv(x = df_F1, file = paste0(simulation_dir, "/df_F1_revision.csv"))










