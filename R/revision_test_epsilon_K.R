# This file runs the simulation to show the sensitivity of CompressiveNMF to 
# the values of epsilon, K, a and alpha.  

library(CompressiveNMF)
library(tidyverse)
library(mcmcse)
library(foreach)
library(lsa)
library(doParallel)

source("R/Postprocess_functions.R")

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
  cosmic_data <- CompressiveNMF::COSMIC_SBS96_GRCh37_v3.4
  Rmat_random <- t(LaplacesDemon::rdirichlet(n = K_new, alpha = rep(alpha, 96)))
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

test_epsilon_k <- function(K_range, 
                           epsilon_range, 
                           J = 50, 
                           ndatasets = 20, 
                           overdispersion = 0, nsamples = 50, burnin = 100){
  values <- expand.grid("K" = K_range, "epsilon" = epsilon_range)
  registerDoParallel(min(nrow(values), parallel::detectCores() - 1))
  results <- foreach(s = c(1:nrow(values)), .combine = "rbind") %dopar% {
    # Simulate a dataset
    data_all <- lapply(1:ndatasets, function(x) simulate_data(K_new = 2, 
                                                              J = J, 
                                                              overdispersion = overdispersion))
    res <- lapply(data_all, FUN = function(data) CompressiveNMF(data$X, ncores = 1, nchains = 1,
                                                                K = values$K[s],  
                                                                epsilon = values$epsilon[s],
                                                                cutoff_excluded = 0,
                                                                burnin =  burnin,
                                                                nsamples = nsamples))
    
    output <- lapply(1:ndatasets, function(i) {
      Lambda <- get_Lambda_Comp(res[[i]])
      Lambda_true <- data_all[[i]]$Rmat %*% data_all[[i]]$Theta
      rmse_Lambda <- sqrt(mean((Lambda_true - Lambda)^2))
      rmse_Counts <- sqrt(mean((data_all[[i]]$X - Lambda)^2))
      # Find estimated signatures
      id_selected <- res[[i]]$RelWeights > 1.5 * values$epsilon[s]
      K_selected <- sum(id_selected)
      # Compare true vs estimated
      R_true <- apply(data_all[[i]]$Rmat, 2, function(x) x / sum(x))
      R_hat <- res[[i]]$Signatures[, id_selected]
      matchedSign <- match_MutSign(R_true = R_true, R_hat = R_hat)
      cos_sim <- mean(get_cosine_similarity(matchedSign))
      # Step 4 - calculate the RMSE between Theta and the rest
      Theta_hat <- res[[i]]$Weights[id_selected, ]
      rmse_R <- compute_RMSE_Signature(R_hat = matchedSign$R_hat, R_true = matchedSign$R_true)  
      rmse_Theta <- compute_RMSE_Theta(Theta_true = data_all[[i]]$Theta, Theta_hat = Theta_hat, matchedSign$match)  
      # Step 5 - calculate the sensitivity and precision
      sens_prec  <- Compute_sensitivity_precision(R_hat = R_hat, R_true)
      
      results <- c("Kest" = K_selected,
                   "mean" = mean(res[[i]]$RelWeights[!id_selected]), 
                   "lowCI" = quantile(res[[i]]$RelWeights[!id_selected], probs = c(0.05)), 
                   "highCI" = quantile(res[[i]]$RelWeights[!id_selected], probs = c(0.95)), 
                   sens_prec, 
                   rmse_R, 
                   rmse_Theta, 
                   "mean_cos_sim" = cos_sim, 
                   "rmse_Lambda" = rmse_Lambda, 
                   "rmse_Counts" = rmse_Counts)
    })
    output<- data.frame(do.call("rbind", output))
    output$epsilon = values$epsilon[s]
    output$Kused = values$K[s]
    output
  }
  return(results)
}

epsilon_range <- c(0.001, 0.01, 0.1, 0.25, 0.5, 1, 1.5, 2)
K_range <- c(5, 10, 20, 30, 40, 50)


ndatasets <- 20
nsamples <- 1000
burnin <- 3000
J <- 100

set.seed(10, kind = "L'Ecuyer-CMRG")
results_correct <- test_epsilon_k(K_range = K_range, 
                                  epsilon_range = epsilon_range, 
                                  J = J, 
                                  ndatasets = ndatasets, 
                                  overdispersion = 0,
                                  nsamples = nsamples, burnin = burnin)
results_correct$overdispersion <- 0

set.seed(10, kind = "L'Ecuyer-CMRG")
results_misp <- test_epsilon_k(K_range = K_range,
                               epsilon_range = epsilon_range,
                               J = J,
                               ndatasets = ndatasets,
                               overdispersion = 0.15,
                               nsamples = nsamples, burnin = burnin)
results_misp$overdispersion <- 0.15

write_csv(rbind(results_correct, results_misp), 
          file = "output/sensitivity_simulation/sensitivity_epsilon_K_100.csv")


df_Keps <-read_csv("output/sensitivity_simulation/sensitivity_epsilon_K_100.csv")

# 
# dodge_eps <- data.frame("Kused" = as.factor(c(1:5) * 10), 
#                       "dodge" = c(-0.02, -0.01, 0, 0.01, 0.02))
# 
# ggplot(df_Keps %>%
#        mutate(Kused = as.factor(Kused)) %>%
#        left_join(dodge_eps, by = "Kused")) +
# geom_abline(slope = 1, intercept = 0, col = "grey") +
# geom_point(aes(x = epsilon, y = mean, color = Kused), 
#            position = position_dodge(width = 0.05))+
# theme_bw()
# 
# 
# ## Plot of excluded columns
df_Keps %>%
  filter(Kused >5) %>%
  mutate(Kused = as.factor(Kused),
         tau = as.factor(overdispersion)) %>%
  group_by(epsilon, Kused, tau) %>%
  summarise(mean_eps_excluded = mean(mean)) %>%
  #left_join(dodge_eps, by = "Kused") %>%
  ggplot() +
  geom_abline(slope = 1, intercept = 0, col = "grey") +
  geom_segment(x = 1, xend = 1, y = 0, yend = Inf, col = "grey50", linetype = "dashed")+
  geom_point(aes(x = epsilon, y = mean_eps_excluded, color = tau))+
  theme_bw() +
  facet_grid(~ Kused)+
  theme(aspect.ratio = 1)

# 
# 
## RMSE as a function of K and epsilon
df_Keps %>%
  mutate(epsilon = as.factor(epsilon),
         overdispersion = paste0("overd = ", overdispersion)) %>%
  group_by(epsilon, Kused, overdispersion) %>%
  summarise(rmse = mean(rmse_Counts)) %>%
  ggplot() +
  geom_point(aes(x = Kused, y = rmse, color = epsilon)) +
  geom_line(aes(x = Kused, y = rmse, color = epsilon)) +
  theme_bw() +
  facet_wrap(~overdispersion, scales = "free")

## Number of factors as a function of K and epsilon
df_Keps %>%
  mutate(epsilon = as.factor(epsilon),
         overdispersion = paste0("overd = ", overdispersion)) %>%
  group_by(epsilon, Kused, overdispersion) %>%
  summarise(Kest = mean(Kest)) %>%
  ggplot() +
  geom_point(aes(x = Kused, y = Kest, color = epsilon)) +
  geom_line(aes(x = Kused, y = Kest, color = epsilon)) +
  theme_bw() +
  facet_wrap(~overdispersion, scales = "free")
# 
# 
# 
df_Keps %>%
mutate(epsilon = as.factor(epsilon)) %>%
group_by(epsilon, Kused, overdispersion) %>%
summarise(rmse_Weights = mean(rmse_Weights)) %>%
ggplot()+
geom_line(aes(x = Kused, y = rmse_Weights, color = epsilon))+
  facet_wrap(~overdispersion, scales = "free")

df_Keps %>%
  mutate(epsilon = as.factor(epsilon)) %>%
  group_by(epsilon, Kused, overdispersion) %>%
  summarise(rmse = mean(rmse_Signatures)) %>%
  ggplot()+
  geom_line(aes(x = Kused, y = rmse, color = epsilon))+
  facet_wrap(~overdispersion, scales = "free")


# 
# # Let's now plot the rmse 
# ggplot(results_correct %>%
#        mutate(Kused = as.factor(Kused))) +
#        #geom_abline(slope = 1, intercept = 0, col = "grey") +
#        geom_point(aes(x= Kused, y = rmse_Lambda, color = as.factor(epsilon)))+
#        geom_line(aes(x= Kused, y = rmse_Lambda, color = as.factor(epsilon)))+
#      theme_bw()
#      
#      
#      ggplot(results_correct %>%
#               mutate(Kused = as.factor(Kused),
#                      epsilon = as.factor(epsilon))) +
#        #geom_abline(slope = 1, intercept = 0, col = "grey") +
#        geom_boxplot(aes(x= Kused, y = rmse_Signatures, color = epsilon))+
#        theme_bw()
#      
#      
#      ggplot(results_misp %>%
#               mutate(Kused = as.factor(Kused))) +
#        geom_abline(slope = 1, intercept = 0, col = "grey") +
#        geom_point(aes(x= epsilon, y = mean, color = Kused))+
#        theme_bw()
#      
#      ggplot(results_correct %>%
#               mutate(Kused = as.factor(Kused),
#                      epsilon = as.factor(epsilon))) +
#        #geom_abline(slope = 1, intercept = 0, col = "grey") +
#        geom_boxplot(aes(y= Kest, x= epsilon, color = epsilon))+
#        facet_wrap(~Kused, nrow = 1)+
#        theme_bw()
#      
#      
# 
# 
