#------------------------------------------------- Figure S7.6 in the supplement
#
# This file further test the sensitivity of the approach when using different 
# values of alpha and a. We pick K = 20 and epsilon = 0.001, just like the 
# simulation in the main manuscript. 

# This file runs the simulation to show the sensitivity of CompressiveNMF to 
# the values of epsilon, K, a and alpha.  
library(CompressiveNMF)
library(tidyverse)
library(mcmcse)
library(ggpubr)
library(grid)
library(foreach)
library(lsa)
library(doParallel)

source("R/Postprocess_functions.R")
source("R/simulate_data.R")

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

test_alpha_a <- function(a_range, 
                         alpha_range, 
                         J = 50, 
                         ndatasets = 20, 
                         overdispersion = 0, 
                         nsamples = 50,
                         K = 20, 
                         epsilon = 0.001,
                         burnin = 100, 
                         num_cores = NULL){
  values <- expand.grid("a" = a_range, "alpha" = alpha_range)
  
  # Create the cluster
  if(is.null(num_cores)){
    num_cores <- min(nrow(values), parallel::detectCores() - 1)
  } 
  #cl <- makeCluster(num_cores, )
  registerDoParallel(num_cores)
  # Simulate a dataset
  data_all <- lapply(1:ndatasets, 
                     function(x) simulate_data(K_new = 2, 
                                               J = J,
                                               overdispersion = overdispersion))
  
  # Parallelize
  results <- foreach(s = c(1:nrow(values)), .combine = "rbind") %dopar% {

    output <- vector(mode = "list", length = ndatasets)
    
    for(j in 1:ndatasets){
      data <- data_all[[j]]
      # Run compressiveNMF
      res <- CompressiveNMF(data$X, 
                                 ncores = 1, nchains = 1,
                                 K = K, 
                                 alpha = values$alpha[s],
                                 a = values$a[s],
                                 epsilon = epsilon,
                                 burnin =  burnin,
                                 nsamples = nsamples)
      # Get results
      Lambda <- get_Lambda_Comp(res)
      Lambda_true <- data$Rmat %*% data$Theta
      rmse_Lambda <- sqrt(mean((Lambda_true - Lambda)^2))
      rmse_Counts <- sqrt(mean((data$X - Lambda)^2))
      # Find estimated signatures
      id_selected <- res$RelWeights > 1.5 * epsilon
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
      rmse_Counts_filer <- sqrt(mean((data$X - R_hat %*% Theta_hat)^2))
      rmse_Lambda_filer <- sqrt(mean((Lambda_true - R_hat %*% Theta_hat)^2))
      output[[j]]  <- c("Kest" = K_selected,
                   "mean" = mean(res$RelWeights[!id_selected]), 
                   "lowCI" = quantile(res$RelWeights[!id_selected], probs = c(0.05)), 
                   "highCI" = quantile(res$RelWeights[!id_selected], probs = c(0.95)), 
                   sens_prec, 
                   rmse_R, 
                   rmse_Theta, 
                   "mean_cos_sim" = cos_sim, 
                   "rmse_Lambda" = rmse_Lambda, 
                   "rmse_Counts" = rmse_Counts,
                   "rmse_Lambda_filter" = rmse_Lambda_filer, 
                   "rmse_Counts_filter" = rmse_Counts_filer)
    }
    output<- data.frame(do.call("rbind", output))
    output$epsilon = epsilon
    output$Kused = K
    output$a = values$a[s]
    output$alpha = values$alpha[s]
    output
  }
  return(results)
}

# Set up the grid of values
a_range <- c(0.1, 0.5, 1, 2)
alpha_range <- c(0.01, 0.1, 0.5, 1, 2)

ndatasets <- 40
nsamples <- 500
burnin <- 3000
J <- 100

rerun <- FALSE
if(rerun) {
  set.seed(10, kind = "L'Ecuyer-CMRG")
  results_correct <- test_alpha_a(a_range = a_range, alpha_range = alpha_range,
                                  K = 20, epsilon = 0.001,
                                  J = J, 
                                  ndatasets = ndatasets, 
                                  overdispersion = 0,
                                  nsamples = nsamples, 
                                  burnin = burnin, 
                                  num_cores = 10)
  results_correct$overdispersion <- 0
  write_csv(results_correct, 
            file = "output/sensitivity_simulation/sensitivity_a_alpha_100_correct.csv")
  
  set.seed(10, kind = "L'Ecuyer-CMRG")
  results_misp <- test_alpha_a(a_range = a_range, alpha_range = alpha_range,
                               K = 20, epsilon = 0.001,
                               J = J, 
                               ndatasets = ndatasets, 
                               overdispersion = 0.15,
                               nsamples = nsamples, 
                               burnin = burnin, 
                               num_cores = 10)
  results_misp$overdispersion <- 0.15
  write_csv(results_misp, 
            file = "output/sensitivity_simulation/sensitivity_a_alpha_100_misp.csv")
  
}

################################################ 
# Figure S7.6 in the supplement
################################################ 
df_a_alpha <- rbind(read_csv("output/sensitivity_simulation/sensitivity_a_alpha_100_correct.csv"),
                 read_csv("output/sensitivity_simulation/sensitivity_a_alpha_100_misp.csv"))

p_K <- df_a_alpha %>%
  mutate(a = as.factor(a),
         overdispersion = paste0("tau = ", overdispersion)) %>%
  group_by(a, alpha, overdispersion) %>%
  summarise(`Estimated K` = mean(Kest)) %>%
  gather(key = "quant", value = "value", -a, -alpha, -overdispersion) %>%
  ggplot() +
  geom_hline(yintercept =  6, linetype = "dashed", color = "gray50") +
  geom_point(aes(x = alpha, y = value, color = a)) +
  geom_line(aes(x = alpha, y = value, color = a)) +
  theme_bw() +
  facet_wrap(overdispersion~quant, scales = "free_y", nrow = 2)+
  theme(axis.title = element_blank(), 
        panel.spacing.y = unit(0, "lines")) +
  ylim(c(4, 20))

p_load <- df_a_alpha %>%
  mutate(a = as.factor(a),
         overdispersion = paste0("tau = ", overdispersion)) %>%
  group_by(a, alpha, overdispersion) %>%
  summarise(`RMSE Loadings` = mean(rmse_Weights)) %>%
  gather(key = "quant", value = "value", -a, -alpha, -overdispersion) %>%
  ggplot() +
  geom_point(aes(x = alpha, y = value, color = a)) +
  geom_line(aes(x = alpha, y = value, color = a)) +
  theme_bw() +
  facet_wrap(overdispersion~quant, scales = "free_y", nrow = 2)+
  theme(axis.title = element_blank(), 
        panel.spacing.y = unit(0, "lines"))

p_sig <- df_a_alpha %>%
  mutate(a = as.factor(a),
         overdispersion = paste0("tau = ", overdispersion)) %>%
  group_by(a, alpha, overdispersion) %>%
  summarise(`RMSE Signatures` = mean(rmse_Signatures)) %>%
  gather(key = "quant", value = "value", -a, -alpha, -overdispersion) %>%
  ggplot() +
  geom_point(aes(x = alpha, y = value, color = a)) +
  geom_line(aes(x = alpha, y = value, color = a)) +
  theme_bw() +
  facet_wrap(overdispersion~quant, scales = "free_y", nrow = 2)+
  theme(axis.title = element_blank(), 
        panel.spacing.y = unit(0, "lines"))

figure <- ggpubr::ggarrange(p_K, p_load, p_sig, common.legend = TRUE, nrow = 1)

figure <- annotate_figure(figure, 
                          bottom = textGrob("Value of alpha", gp = gpar(cex = 1)))
ggsave(plot = figure, filename = "figures/sensitivty_a_alpha_new.pdf", 
       width = 7.58, height = 4.45)


