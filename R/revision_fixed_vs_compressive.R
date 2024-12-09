# Simulation - impact of strength-matching vs fixed choice of hyperprior. 
# This will show the benefit of having the compressive property

library(CompressiveNMF)
library(tidyverse)
library(mcmcse)
library(foreach)
library(lsa)
library(doParallel)
library(ggpubr)
library(simplecolors)

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


# Hyperparameters of the simulation *(fixed vs compressive)
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
                      overdispersion = overd))
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
                                                nburn = 3000, 
                                                ncores = 1){
  
  subdir <- paste0(main_dir, "Scenario_", J, "_overd_", overd, "/")
  # Load the data
  data_all <- readRDS(paste0(subdir, "data.rds.gzip"))
  # Run all models in parallel
  registerDoParallel(ncores)
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
nsamples <- 1000
nburn <- 3000
ncores <- ndatasets

#-------------------------------------------- Run compressive case
run_compressive <- FALSE
run_correct <- FALSE
run_overd <- FALSE

if(run_compressive){
  for(J in J_range){
    
    # Correctly specified case
    if(run_correct) {
      overd <- 0
      print(paste0("Compressive: J - ", J, ", overd - ", overd))
      set.seed(10, kind = "L'Ecuyer-CMRG")
      run_simulation_fixed_vs_compressive(J = J, 
                                          overd = overd, 
                                          case = "compressive", 
                                          nsamples = nsamples, 
                                          nburn = nburn, 
                                          ncores = ncores)
    }
    
    # Incorrect case
    if(run_overd) {
      overd <- 0.15
      print(paste0("Compressive: J - ", J, ", overd - ", overd))
      set.seed(10, kind = "L'Ecuyer-CMRG")
      run_simulation_fixed_vs_compressive(J = J, 
                                          overd = overd, 
                                          case = "compressive", 
                                          nsamples = nsamples, 
                                          nburn = nburn, 
                                          ncores = ncores)
    }
  }
}

#-------------------------------------------- Run fixed case
run_fixed <- FALSE
run_correct <- FALSE
run_overd <- FALSE

if(run_fixed){
  for(J in J_range){
    
    # Correctly specified case
    if(run_correct) {
      overd <- 0
      print(paste0("Fixed: J - ", J, ", overd - ", overd))
      set.seed(10, kind = "L'Ecuyer-CMRG")
      run_simulation_fixed_vs_compressive(J = J, 
                                          overd = overd, 
                                          case = "fixed", 
                                          nsamples = nsamples, 
                                          nburn = nburn, 
                                          ncores = ncores)
    }
    
    # Incorrect case
    if(run_overd) {
      overd <- 0.15
      print(paste0("Fixed: J - ", J, ", overd - ", overd))
      set.seed(10, kind = "L'Ecuyer-CMRG")
      run_simulation_fixed_vs_compressive(J = J, 
                                          overd = overd, 
                                          case = "fixed", 
                                          nsamples = nsamples, 
                                          nburn = nburn, 
                                          ncores = ncores)
    }
  }
}


#-------------------------------------------- 
# Plot the functions
#-------------------------------------------- 

df_results <- data.frame()
main_dir <- "~/CompressiveNMF/output/compressive_vs_fixed_simulation/"
for(J in J_range){
  for(overd in overd_list){
    subdir <- paste0(main_dir, "Scenario_", J, "_overd_", overd, "/")
    # Read csv results
    res_compressive <- as.data.frame(read_csv(paste0(subdir, "results_compressive.csv"), show_col_types = FALSE))
    res_fixed <- as.data.frame(read_csv(paste0(subdir, "results_fixed.csv"), show_col_types = FALSE))
    # Append
    df_results <- rbind(df_results, res_compressive,  res_fixed)
  }
}

# #------------------- Number of signatures, precision and sensitivity
# p_K_prec_sens_0 <- df_results %>%
#   filter(overd == 0) %>%
#   mutate(Jlist = as.factor(J), 
#          ovd = paste0(ovd = paste0("tau = ", overd)), 
#          `Estimated K` = Kest) %>%
#   dplyr::select(Jlist, ovd, case, `Estimated K`, Precision, Sensitivity) %>%
#   gather(key = "quant", value = "value", -Jlist, -ovd, -case) %>%
#   mutate(Ktrue = case_when(quant == "Estimated K" ~ 10,
#                            TRUE ~ NA)) %>%
#   ggplot()+
#     theme_bw()+
#     geom_hline(aes(yintercept = Ktrue), color = "turquoise", linetype = "solid")+
#     geom_boxplot(aes(x = Jlist, y = value, fill = case, color = case)) +
#     facet_grid(quant~ovd, scales = "free", shrink = FALSE)+
#     scale_color_manual("Prior", values = c("firebrick", "blue"))+
#      scale_fill_manual("Prior", values = c("#F8766D", "#619CFF")) +
#     theme(axis.title.y = element_blank(), 
#           panel.spacing = unit(0, "lines"))+
#   xlab("Sample size J")
# 
# p_K_prec_sens_015 <- df_results %>%
#   filter(overd == 0.15) %>%
#   mutate(Jlist = as.factor(J), 
#          ovd = paste0(ovd = paste0("tau = ", overd)), 
#          `Estimated K` = Kest) %>%
#   dplyr::select(Jlist, ovd, case, `Estimated K`, Precision, Sensitivity) %>%
#   gather(key = "quant", value = "value", -Jlist, -ovd, -case) %>%
#   mutate(Ktrue = case_when(quant == "Estimated K" ~ 10,
#                            TRUE ~ NA)) %>%
#   ggplot()+
#   theme_bw()+
#   geom_hline(aes(yintercept = Ktrue), color = "turquoise", linetype = "solid")+
#   geom_boxplot(aes(x = Jlist, y = value, fill = case, color = case)) +
#   facet_grid(quant~ovd, scales = "free", shrink = FALSE)+
#   scale_color_manual("Prior", values = c("firebrick", "blue"))+
#   scale_fill_manual("Prior", values = c("#F8766D", "#619CFF")) +
#   theme(axis.title.y = element_blank(), 
#         panel.spacing = unit(0, "lines"))+
#   xlab("Sample size J")

p_K_prec_sens <- df_results %>%
  mutate(Jlist = as.factor(J), 
         ovd = paste0(ovd = paste0("tau = ", overd)), 
         `Estimated K` = Kest, 
         case = case_when(case == "compressive"~"Compressive",
                          TRUE~"Fixed-strength")) %>%
  dplyr::select(Jlist, ovd, case, `Estimated K`, Precision, Sensitivity) %>%
  gather(key = "quant", value = "value", -Jlist, -ovd, -case) %>%
  mutate(Ktrue = case_when(quant == "Estimated K" ~ 10,
                           TRUE ~ NA)) %>%
  ggplot()+
  theme_bw()+
  geom_hline(aes(yintercept = Ktrue), color = "turquoise", linetype = "solid")+
  geom_boxplot(aes(x = Jlist, y = value, fill = case, color = case)) +
  facet_wrap(ovd~quant, scales = "free_y", shrink = FALSE)+
  scale_color_manual("Prior", values = c("firebrick", "blue"))+
  scale_fill_manual("Prior", values = c("#F8766D", "#619CFF")) +
  theme(axis.title.y = element_blank(), 
        panel.spacing.y = unit(0, "lines"), 
        legend.position = "top")+
  xlab("Sample size J")
p_K_prec_sens

ggsave(plot = p_K_prec_sens, 
       filename = "figures/simulation_fixed_vs_compressive_K_prec_sens.pdf", 
       width = 8.41, height = 5.27)

#------------------- RMSE of all quantities
# p_rmse_0 <- df_results %>%
#   filter(overd == 0) %>%
#   mutate(Jlist = as.factor(J), 
#          ovd = paste0(ovd = paste0("tau = ", overd)), 
#          `RMSE Counts` = rmse_Counts, 
#          `RMSE Loadings` = rmse_Weights, 
#          `RMSE Signatures` = rmse_Signatures) %>%
#   dplyr::select(Jlist, ovd, case, `RMSE Counts`, `RMSE Loadings`, `RMSE Signatures`) %>%
#   gather(key = "quant", value = "value", -Jlist, -ovd, -case) %>%
#   ggplot()+
#   theme_bw()+
#   geom_boxplot(aes(x = Jlist, y = value, fill = case, color = case)) +
#   facet_grid(quant~ovd, scales = "free", shrink = FALSE)+
#   scale_color_manual("Prior", values = c("firebrick", "blue"))+
#   scale_fill_manual("Prior", values = c("#F8766D", "#619CFF")) +
#   theme(axis.title.y = element_blank(), 
#         panel.spacing = unit(0, "lines"))+
#   xlab("Sample size J")
# 
# p_rmse_015 <- df_results %>%
#   filter(overd == 0.15) %>%
#   mutate(Jlist = as.factor(J), 
#          ovd = paste0(ovd = paste0("tau = ", overd)), 
#          `RMSE Counts` = rmse_Counts, 
#          `RMSE Loadings` = rmse_Weights, 
#          `RMSE Signatures` = rmse_Signatures) %>%
#   dplyr::select(Jlist, ovd, case, `RMSE Counts`, `RMSE Loadings`, `RMSE Signatures`) %>%
#   gather(key = "quant", value = "value", -Jlist, -ovd, -case) %>%
#   ggplot()+
#   theme_bw()+
#   geom_boxplot(aes(x = Jlist, y = value, fill = case, color = case)) +
#   facet_grid(quant~ovd, scales = "free", shrink = FALSE)+
#   scale_color_manual("Prior", values = c("firebrick", "blue"))+
#   scale_fill_manual("Prior", values = c("#F8766D", "#619CFF")) +
#   theme(axis.title.y = element_blank(), 
#         panel.spacing = unit(0, "lines"))+
#   xlab("Sample size J")

p_rmse <- df_results %>%
  mutate(Jlist = as.factor(J), 
         ovd = paste0(ovd = paste0("tau = ", overd)), 
         `RMSE Counts` = rmse_Counts, 
         `RMSE Loadings` = rmse_Weights, 
         `RMSE Signatures` = rmse_Signatures, 
         case = case_when(case == "compressive"~"Compressive",
                          TRUE~"Fixed-strength")) %>%
  dplyr::select(Jlist, ovd, case, `RMSE Counts`, `RMSE Loadings`, `RMSE Signatures`) %>%
  gather(key = "quant", value = "value", -Jlist, -ovd, -case) %>%
  ggplot()+
  theme_bw()+
  geom_boxplot(aes(x = Jlist, y = value, fill = case, color = case)) +
  facet_wrap(ovd~quant, scales = "free_y")+
  scale_color_manual("Prior", values = c("firebrick", "blue"))+
  scale_fill_manual("Prior", values = c("#F8766D", "#619CFF")) +
  theme(axis.title.y = element_blank(), 
        panel.spacing.y = unit(0, "lines"), legend.position = "top")+
  xlab("Sample size J")

#ggpubr::ggarrange(p_K_prec_sens_0, p_K_prec_sens_015, p_rmse_0, p_rmse_015, common.legend = TRUE, legend = "top", nrow = 1)

#ggpubr::ggarrange(p_rmse, common.legend = TRUE, legend = "top", nrow = 1)
ggsave(plot = p_rmse, 
       filename = "figures/simulation_fixed_vs_compressive_rmse.pdf", 
       width = 8.41, height = 5.27)



#ggpubr::ggarrange(p_K_prec_sens, p_rmse, common.legend = TRUE, legend = "top", nrow = 1)
#ggsave(filename = "figures/simulation_fixed_vs_compressive.pdf", 
#       width = 9.30, height = 5.11)


#-------------------------------------------- 
# Show an individual case. How do the signatures and the loadings look like?
data <- readRDS("output/compressive_vs_fixed_simulation/Scenario_300_overd_0.15/data.rds.gzip")[[3]]

# Compressive case
set.seed(10)
out_comp <- CompressiveNMF(X = data$X,
                      a0 = ncol(data$X) + 1, 
                      b0 = epsilon * ncol(data$X), 
                      cutoff_excluded = 0,
                      K = K, 
                      epsilon = epsilon,
                      nsamples = nsamples, 
                      burnin = nburn, 
                      verbose = TRUE)

# Fixed case
set.seed(10)
out_fixed <- CompressiveNMF(X = data$X,
                        a0 = strength_fixed + 1, 
                        b0 = epsilon * strength_fixed, 
                        cutoff_excluded = 0,
                        K = K, 
                        epsilon = epsilon,
                        nsamples = nsamples, 
                        burnin = nburn, 
                        verbose = TRUE)


# Match them via Hungarian algorithm for better visualization
match_mut <- match_MutSign(R_hat = out_comp$Signatures, 
                           R_true = out_fixed$Signatures)
sig_compressive <- match_mut$R_hat
sig_fixed <- match_mut$R_true

# Put the "compressed out signatures" at the bottom. 
is_flat <- apply(sig_compressive, 2, function(x) cosine(x, rep(1, 96)) > 0.95)
sig_compressive <- cbind(sig_compressive[, !is_flat], sig_compressive[, is_flat])
sig_fixed <- cbind(sig_fixed[, !is_flat], sig_fixed[, is_flat])
match <- c(match_mut$match[!is_flat], match_mut$match[is_flat])

# Name the signatures
colnames(sig_compressive) <- LETTERS[1:K]
colnames(sig_fixed) <- LETTERS[1:K]

# Plot signatures side by side
p_sig <-ggarrange(CompressiveNMF:::plot.SBS.signature(sig_compressive) +
            theme(axis.text.x = element_blank(), 
                  plot.title = element_text(hjust = 0.5)) +
            ggtitle("Compressive"), 
          CompressiveNMF:::plot.SBS.signature(sig_fixed) +
            theme(axis.text.x = element_blank(), 
                  plot.title = element_text(hjust = 0.5))+
            ggtitle("Fixed-strength"), 
          nrow = 1)

ggsave(plot = p_sig , 
       filename = "figures/simulation_fixed_vs_compressive_signature_ex.pdf", 
       width = 8.41, height = 5.19)

# Plot the relevance weights
mu_comp <- data.frame(out_comp$mcmc_out[[1]]$Mu[, match])
colnames(mu_comp) <- LETTERS[1:K]
mu_fixed <- data.frame(out_fixed$mcmc_out[[1]]$Mu)
mu_fixed <- cbind(mu_fixed[, !is_flat], mu_fixed[, is_flat])
colnames(mu_fixed) <- LETTERS[1:K]

colors_palette_20 <- c("black",
                       rev(RColorBrewer::brewer.pal(n = 9, name = "Blues")), 
                       "gray",
                       RColorBrewer::brewer.pal(n = 9, name = "Reds"))
p_trace_comp <- mu_comp %>%
  mutate(iter = 1:nrow(mu_comp), 
         case = "Compressive") %>%
  gather(key = "sig", value = "mu", -iter, -case) %>%
  ggplot()+
  #theme_bw()+
  geom_line(aes(x = iter, y = mu, col = sig))+
  scale_color_manual(values = colors_palette_20)+
  xlab("Post burn-in MCMC iteration")+
  facet_grid(.~case, scales = "free_y")

p_trace_fixed <- mu_fixed %>%
  mutate(iter = 1:nrow(mu_fixed), 
         case = "Fixed-strength") %>%
  gather(key = "sig", value = "mu", -iter, -case) %>%
  ggplot()+
  #theme_bw()+
  geom_line(aes(x = iter, y = mu, col = sig))+
  scale_color_manual(values = colors_palette_20)+
  xlab("Post burn-in MCMC iteration")+
  facet_grid(.~case, scales = "free_y")

p_trace_all <- rbind(mu_comp %>%
        mutate(iter = 1:nrow(mu_comp), 
               case = "Compressive") %>%
        gather(key = "sig", value = "mu", -iter, -case), 
      mu_fixed %>%
        mutate(iter = 1:nrow(mu_fixed), 
               case = "Fixed-strength") %>%
        gather(key = "sig", value = "mu", -iter, -case)) %>%
  ggplot()+
  #theme_bw()+
  geom_line(aes(x = iter, y = mu, col = sig))+
  scale_color_manual("Signature", values = colors_palette_20)+
  xlab("Post burn-in MCMC iteration")+
  facet_grid(.~case, scales = "free_y") +
  theme(legend.position = "top")+
  guides(colour = guide_legend(nrow = 2))+
  ylab(expression(paste("Posterior samples of ", mu[k])))

ggsave(plot = p_trace_all , 
       filename = "figures/simulation_fixed_vs_compressive_traceplots.pdf", 
       width = 8.41, height = 4.06)
