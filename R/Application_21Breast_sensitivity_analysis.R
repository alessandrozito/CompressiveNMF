# This file runs a sensitivity analysis for the 21 breast cancer in 
# CompressiveNMF and CompressiveNMF + cosmic.

library(CompressiveNMF)
library(doParallel)
library(tidyverse)
library(foreach)
library(lsa)
library(ggpubr)


source("R/Postprocess_functions.R")

# Load the 21 breast cancer
load("data/21breast.rdata")


# Possible range of values for the sensitivity analysis of the 21 breast cancer
a_range <- c(1, 2, 0.5)
alpha_range <- c(0.5, 1, 2)
epsilon_range <- c(0.01, 0.1, 0.25)
K_range <- c(15, 20, 25, 40)
values_all <- expand_grid(a_range, alpha_range, epsilon_range, K_range)
  
# Parameters of the simulation
nsamples <- 2000
burnin <- 10000

rerun <- FALSE # <------- set to TRUE to re-run the study
if(rerun){
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
    alpha <- values_all$alpha_range[i]
    epsilon <- values_all$epsilon_range[i]
    
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
}

# Analyze the ouput
df_results <-read_csv("output/Application_21brca/Sensitivity_21brca_CompNMF.csv")

pK <- ggplot(df_results %>%
         mutate(a = paste0("a = ", a), 
                epsilon = paste0("eps = ", epsilon),
                alpha = paste0("alpha = ", alpha))) +
  geom_point(aes(x = as.factor(K), y = Kest, shape = epsilon, color = epsilon), 
             size = 2.5, stroke = .8) +
  scale_shape_manual(values = c(1,3,4))+
  facet_grid(a~alpha)+
  ylim(c(3,10))+
  theme_bw()+
  theme(panel.spacing = unit(0, "lines"))+
  ylab("N. estimated signatures") +
  xlab("K")
ggsave(pK, filename = "figures/Sensitivity_K_21breast.pdf", width = 6.85, height = 3.21)


p1 <- ggplot(df_results %>%
         mutate(a = paste0("a = ", a), 
                epsilon = paste0("eps = ", epsilon),
                alpha = paste0("alpha = ", alpha))) +
  geom_point(aes(x = as.factor(K), y = Precision, shape = epsilon, color = epsilon), 
             size = 2.5, stroke = .8) +
  scale_shape_manual(values = c(1,3,4))+
  facet_grid(a~alpha)+
  #ylim(c(3,10))+
  theme_bw()+
  ylab("Precision") +
  xlab("K")+
  ggtitle("Precision")

p2 <- ggplot(df_results %>%
         mutate(a = paste0("a = ", a), 
                epsilon = paste0("eps = ", epsilon),
                alpha = paste0("alpha = ", alpha))) +
  geom_point(aes(x = as.factor(K), y = Sensitivity, shape = epsilon, color = epsilon), 
             size = 2.5, stroke = .8) +
  scale_shape_manual(values = c(1,3,4))+
  facet_grid(a~alpha)+
  #ylim(c(3,10))+
  theme_bw()+
  ylab("Sensitivity") +
  xlab("K")+
  ggtitle("Sensitivity")

ggarrange(p1, p2, common.legend = TRUE)


p3 <- ggplot(df_results %>%
               mutate(a = paste0("a = ", a), 
                      epsilon = paste0("eps = ", epsilon),
                      alpha = paste0("alpha = ", alpha))) +
  geom_point(aes(x = as.factor(K), y = rmse_Signatures, shape = epsilon, color = epsilon), 
             size = 2.5, stroke = .8) +
  scale_shape_manual(values = c(1,3,4))+
  facet_grid(a~alpha)+
  #ylim(c(3,10))+
  theme_bw()+
  ylab("RMSE") +
  xlab("K")+
  theme(panel.spacing = unit(0, "lines"))+
  ggtitle("Signatures")

p4 <- ggplot(df_results %>%
               mutate(a = paste0("a = ", a), 
                      epsilon = paste0("eps = ", epsilon),
                      alpha = paste0("alpha = ", alpha))) +
  geom_point(aes(x = as.factor(K), y = rmse_Weights, shape = epsilon, color = epsilon), 
             size = 2.5, stroke = .8) +
  scale_shape_manual(values = c(1,3,4))+
  facet_grid(a~alpha)+
  #ylim(c(3,10))+
  theme_bw()+
  ylab("RMSE") +
  xlab("K")+
  theme(panel.spacing = unit(0, "lines"))+
  ggtitle("Loadings")

p_RMSE <- ggarrange(p3, p4, common.legend = TRUE, legend = "top")
ggsave(p_RMSE, filename = "figures/Sensitivity_RMSE_21breast.pdf", 
       width = 8.77, height = 4.82)




