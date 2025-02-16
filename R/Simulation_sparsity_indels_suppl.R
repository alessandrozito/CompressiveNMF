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

#################################################################
# Useful functions
#################################################################
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

open_rds_file <- function(file){
  if(file.exists(file)){
    out <- readRDS(file)
  } else {
    out <- NULL
  }
  return(out)
}

ReshapeResults_BayesNMF_brouwer <- function(i, out_dir, nChannels=96){
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
  R_all <- reshape_python_array(R, dimension = nChannels)
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
  select <- apply(sigs_norm, 2, function(x) cosine(x, rep(1, nChannels))) < 0.975
  
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

#################################################################
# Run the simulation
#################################################################

# Parameters of the simulation
ndatasets <- 20
J <- 100
K_new <- 2
overd_list <- c(0, 0.15)
theta <- 50
# Simulate the array of data
regenerate_data <- FALSE

main_dir <- "~/CompressiveNMF/output/sparse_indel_simulation/"
create_directory(main_dir)

if(regenerate_data){
  # Generate the data with increasing sparsity level
  for(overd in overd_list){
    data_all <- lapply(1:ndatasets, function(i) 
      simulate_data_Indels(J = J, 
                           K_new = K_new, 
                           cosmic_sig = c("ID1", "ID2", "ID8", "ID9"), 
                           alpha = 0.01, 
                           theta = 50, 
                           overdispersion = overd))
    directory <- paste0(main_dir, "Indel_theta_", theta, "_overd_", overd, "/")
    create_directory(directory)
    # Save data
    saveRDS(data_all, file = paste0(directory, "data.rds.gzip"), compress = "gzip")
  }
}

data <- simulate_data_Indels(J = J, 
                     K_new = K_new, 
                     cosmic_sig = c("ID1", "ID2", "ID8", "ID9"), 
                     alpha = 0.01, 
                     theta = 50, 
                     overdispersion = 0)
CompressiveNMF:::plot.ID.signature(data$Rmat)

#######################################################
# Part 1 - Run the simulation
#######################################################
run_compNMF <- FALSE
run_compNMF_cos <- FALSE
run_ARD <- FALSE
run_BayesNMF <- FALSE 
run_SigProfiler <- FALSE

# Number of samples
nsamples <- 1000
burnin <- 3000

registerDoParallel(min(c(ndatasets, parallel::detectCores() - 1)))
#---------------------------------- 1 - CompressiveNMF 
if (run_compNMF) {
  for(overd in overd_list){
    print(overd)
    # Load the data
    out_dir <-  paste0(main_dir, "Indel_theta_", theta, "_overd_", overd, "/")
    data_all <- readRDS(paste0(out_dir, "data.rds.gzip"))
    # Run the model in parallell
    out_CompressiveNMF <- foreach(i = 1:length(data_all)) %dopar% {
      out <- CompressiveNMF(X = data_all[[i]]$X, 
                            K = 20, 
                            epsilon = 0.001, 
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
  for(overd in overd_list){
    print(overd)
    # Load the data
    out_dir <-  paste0(main_dir, "Indel_theta_", theta, "_overd_", overd, "/")
    data_all <- readRDS(paste0(out_dir, "data.rds.gzip"))
    # Run the model in parallell
    out_CompressiveNMF_cos <- foreach(i = 1:length(data_all)) %dopar% {
      out <- CompressiveNMF(X = data_all[[i]]$X, 
                            K = 10,
                            use_cosmic = TRUE,
                            epsilon = 0.001, 
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
  for(overd in overd_list){
    print(overd)
    # Load the data
    out_dir <-  paste0(main_dir, "Indel_theta_", theta, "_overd_", overd, "/")
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
# Step 1 - run the python code

if(run_BayesNMF){

create_directories <- FALSE # <------ Set to TRUE if run for the first time

if(create_directories) {
  for(overd in overd_list){
    print(overd)
    # Set output directory
    out_dir <-  paste0(main_dir,  "Indel_theta_", theta, "_overd_", overd, "/")
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

# <------------- In the console, run the script ------------->
# python R/run_BayesNMF_python/sparse_Indels_simulation.py
# <---------------------------------------------------------->

# Convert the python results to R object
for(overd in overd_list){
  print(overd)
  out_dir <- paste0(main_dir,  "Indel_theta_", theta, "_overd_", overd, "/")
  for(i in 1:20){
    print(i)
    results <- ReshapeResults_BayesNMF_brouwer(i, out_dir, nChannels =83)
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

#---------------------------------- 5 - SigProfiler with R (use Pfrom SigProfilerExtractor import sigpro as sigython script)
# The version iof sigprofiler does not work for some reason, so I need to re-cast
# and run it directly in python
if(run_SigProfiler){
  
  create_directories <- TRUE # <------ Set to TRUE if run for the first time
  
  if(create_directories) {
    for(overd in overd_list){
      print(overd)
      # Set output directory
      out_dir <-  paste0(main_dir,  "Indel_theta_", theta, "_overd_", overd, "/")
      data <- readRDS(paste0(out_dir, "/data.rds.gzip"))
      create_directory(paste0(out_dir, "/SigProfiler/"))
      print(out_dir)
      # Re-save the data as .txt
      out_dir_data <- paste0(out_dir, "/SigProfiler/data/")
      out_dir_output <- paste0(out_dir, "/SigProfiler/output/")
      
      create_directory(out_dir_data)
      create_directory(out_dir_output)
      
      for(i in 1:length(data)){
        # Get the data
        X <- as.data.frame(data[[i]]$X)
        X <- as.data.frame(cbind("Mutation" = rownames(X), X))
        nm <- paste0(out_dir_data, "X_", i, ".txt")
        write_tsv(X, file = nm, )
        # Prepare the output folder
        nm_out <- paste0(out_dir_output, "sim_", i, "/")
        create_directory(nm_out)
      }
    }
  }
}  


# <------------- In the console, run the script ------------->
# python3 R/run_SigProfiler_python/SigProfiler_ID_sparse_simulation.py
# <---------------------------------------------------------->

# Convert the python results to R object
for(overd in overd_list){
  print(overd)
  out_dir <- paste0(main_dir,  "Indel_theta_", theta, "_overd_", overd, "/")
  data_all <- readRDS(paste0(out_dir, "/data.rds.gzip"))
  # Now, postprocess the output and save it
  res <- lapply(1:20, function(i){
    nm <- paste0(out_dir, "SigProfiler/output/sim_", i, 
                 "/ID83/Suggested_Solution/ID83_De-Novo_Solution/")
    # Read the Signatures results
    Signatures <- as.data.frame(suppressMessages(readr::read_table(file = paste0(nm, "Signatures/ID83_De-Novo_Signatures.txt"))))
    channels <- Signatures$MutationType
    Signatures <- as.matrix(Signatures[, -1])
    rownames(Signatures) <- channels
    # Read Activities results
    Weights <- suppressMessages(readr::read_table(file = paste0(nm, "Activities/ID83_De-Novo_Activities_refit.txt")))  
    Weights <- t(as.matrix(Weights[, -1]))
    
    resSigPro <- list("Signatures" = Signatures,
                      "Weights" = Weights,
                      "K" = ncol(Signatures),
                      "time" = 0)
    Postprocess_SigProfiler(resSigPro, data_all[[i]])
  })
  # Save the merged output
  saveRDS(res, file = paste0(out_dir, "/sigProfiler.rds.gzip"), compress = "gzip")
}



#######################################################
# Part 2 - Make plots
#######################################################

# Aggregate output
df_results <- data.frame()
for(overd in overd_list){
  print(overd)
  # Set output directory
  out_dir <-  paste0(main_dir,  "Indel_theta_", theta, "_overd_", overd, "/")
  #------------------------------------------------- CompressiveNMF
  out_CompressiveNMF <- readRDS(paste0(out_dir, "CompressiveNMF.rds.gzip"))
  df_temp <- data.frame(do.call("rbind", lapply(out_CompressiveNMF, 
                                                function(x) x$results)))
  df_temp$Method = "CompNMF" 
  df_temp$overd = overd
  df_results <- rbind(df_results, df_temp)
  #------------------------------------------------- CompressiveNMF + cosmic
  out_CompressiveNMF_cosmic <- readRDS(paste0(out_dir, "CompressiveNMF_cos.rds.gzip"))
  df_temp <- data.frame(do.call("rbind", lapply(out_CompressiveNMF_cosmic, 
                                                function(x) x$results)))
  df_temp$Method = "CompNMFcos" 
  df_temp$overd = overd
  df_results <- rbind(df_results, df_temp)
  #------------------------------------------------- ARD (tan and fevotte, 2011)
  out_ARD <- readRDS(paste0(out_dir, "ARD.rds.gzip"))
  df_temp <- data.frame(do.call("rbind", lapply(out_ARD, 
                                                function(x) add_ess(x$results))))
  df_temp$Method = "SignatureAnalyzer" 
  df_temp$overd = overd
  df_results <- rbind(df_results, df_temp)
  #------------------------------------------------- BayesNMF ARD (Brouwer et al. 2017)
  out_BayesNMF <- readRDS(paste0(out_dir, "BayesNMF_brouwer.rds.gzip"))
  df_temp <- data.frame(do.call("rbind", lapply(out_BayesNMF, 
                                                function(x) add_ess(x$results))))
  df_temp$Method = "BayesNMF" 
  df_temp$overd = overd
  df_results <- rbind(df_results, df_temp)
  #------------------------------------------------- SigProfilerExtractor 
  out_SigPro <- readRDS(paste0(out_dir, "sigProfiler.rds.gzip"))
  df_temp <- data.frame(do.call("rbind", lapply(out_SigPro, 
                                                function(x) add_ess(x$results))))
  df_temp$Method = "sigProfiler" 
  df_temp$overd = overd
  df_results <- rbind(df_results, df_temp)
}

# Make plots
df_results$Method <- factor(df_results$Method, levels = c("CompNMF", "CompNMFcos", "SignatureAnalyzer", "sigProfiler", "BayesNMF"))

pK <- df_results %>%
  dplyr::select(Method, overd, K) %>%
  mutate(overd = paste0("tau = ", overd)) %>%
  gather(key = "quantity", value = "value", -Method, -overd) %>%
  ggplot() +
  geom_hline(yintercept = 6, col = "gray50", linetype = "dashed")+
  geom_boxplot(aes(y = value, x = Method, color = Method, fill = Method), alpha = 0.7)+
  scale_color_manual(values = c("forestgreen", "blue", "firebrick", "darkorange", "black"))+
  scale_fill_manual(values = c("forestgreen", "blue", "firebrick", "darkorange", "black"))+
  ylim(c(3, 20))+
  theme_bw()+
  facet_wrap(~overd+quantity, scales = "free_y", ncol = 1)+
  theme(legend.position = "none", 
        axis.title = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1))
  
pSens <- df_results %>%
  dplyr::select(Method, overd, Sensitivity) %>%
  mutate(overd = paste0("tau = ", overd)) %>%
  gather(key = "quantity", value = "value", -Method, -overd) %>%
  ggplot() +
  geom_boxplot(aes(y = value, x = Method, color = Method, fill = Method), alpha = 0.7)+
  scale_color_manual(values = c("forestgreen", "blue", "firebrick", "darkorange", "black"))+
  scale_fill_manual(values = c("forestgreen", "blue", "firebrick", "darkorange", "black"))+
  theme_bw()+
  facet_wrap(~overd+quantity, scales = "free_y", ncol = 1)+
  theme(legend.position = "none", 
        axis.title = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1))

pPrec <- df_results %>%
  dplyr::select(Method, overd, Precision) %>%
  mutate(overd = paste0("tau = ", overd)) %>%
  gather(key = "quantity", value = "value", -Method, -overd) %>%
  ggplot() +
  geom_boxplot(aes(y = value, x = Method, color = Method, fill = Method), alpha = 0.7)+
  scale_color_manual(values = c("forestgreen", "blue", "firebrick", "darkorange", "black"))+
  scale_fill_manual(values = c("forestgreen", "blue", "firebrick", "darkorange", "black"))+
  theme_bw()+
  facet_wrap(~overd+quantity, scales = "free_y", ncol = 1)+
  theme(legend.position = "none", 
        axis.title = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1))


p1 <- ggarrange(pK, pSens, pPrec, nrow = 1, common.legend = TRUE, legend = "none")
ggsave(p1, filename = "figures/Sparse_indels_simulation_Sens_prec.pdf", height = 5.72, width = 5.59)


p2 <- df_results %>%
  dplyr::select(Method, overd, rmse_Counts, rmse_Weights, rmse_Signatures) %>%
  mutate(overd = paste0("tau = ", overd)) %>%
  gather(key = "quantity", value = "value", -Method, -overd) %>%
  mutate(quantity = case_when(quantity == "rmse_Counts" ~ "RMSE Counts",
                              quantity == "rmse_Signatures" ~ "RMSE Signatures",
                              quantity == "rmse_Weights" ~ "RMSE Loadings")) %>%
  ggplot() +
  geom_boxplot(aes(y = value, x = Method, color = Method, fill = Method), alpha = 0.7)+
  scale_color_manual(values = c("forestgreen", "blue", "firebrick", "darkorange", "black"))+
  scale_fill_manual(values = c("forestgreen", "blue", "firebrick", "darkorange", "black"))+
  theme_bw()+
  facet_wrap(~overd+quantity, scales = "free_y")+
  theme(legend.position = "none", 
        axis.title = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(p2, filename = "figures/Sparse_indels_simulation_RMSE.pdf", height = 5.72, width = 5.59)

ggarrange(p1,p2, nrow = 1, common.legend = TRUE, legend = "none")
ggsave(filename = "figures/Sparse_indels_simulation.pdf", height = 5.72, width = 6 * 2)

