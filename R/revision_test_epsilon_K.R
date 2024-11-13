# This file runs the simulation to show the sensitivity of CompressiveNMF to 
# the values of epsilon, K, a and alpha.  

library(CompressiveNMF)
library(tidyverse)
library(mcmcse)
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


test_epsilon_k <- function(K_range, 
                           epsilon_range, 
                           J = 50, 
                           ndatasets = 20, 
                           overdispersion = 0, 
                           nsamples = 50, 
                           burnin = 100){
  values <- expand.grid("K" = K_range, "epsilon" = epsilon_range)
  
  # Create the cluster
  num_cores <- min(nrow(values), parallel::detectCores() - 1)
  #cl <- makeCluster(num_cores, )
  registerDoParallel(num_cores)
  # Simulate a dataset
  data_all <- lapply(1:ndatasets, function(x) 
    simulate_data(K_new = 2, 
                  J = J, 
                  overdispersion = overdispersion))
  # Parallelize
  results <- foreach(s = c(1:nrow(values)), .combine = "rbind") %dopar% {
    
    #res <- vector(mode = "list", length = ndatasets)
    output <- vector(mode = "list", length = ndatasets)
    
    for(j in 1:ndatasets){
      data <- data_all[[j]]
      # Run the method
      res <- CompressiveNMF(data$X, ncores = 1, nchains = 1,
                     K = values$K[s],  
                     epsilon = values$epsilon[s],
                     cutoff_excluded = 0,
                     burnin =  burnin,
                     nsamples = nsamples)
      # Evaluate results
      Lambda <- get_Lambda_Comp(res)
      Lambda_true <- data$Rmat %*% data$Theta
      rmse_Lambda <- sqrt(mean((Lambda_true - Lambda)^2))
      rmse_Counts <- sqrt(mean((data$X - Lambda)^2))
      # Find estimated signatures
      id_selected <- res$RelWeights > 1.5 * values$epsilon[s]
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
      output[[j]] <- c("Kest" = K_selected,
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
    
    # output <- lapply(1:ndatasets, function(i) {
    #   Lambda <- get_Lambda_Comp(res[[i]])
    #   Lambda_true <- data_all[[i]]$Rmat %*% data_all[[i]]$Theta
    #   rmse_Lambda <- sqrt(mean((Lambda_true - Lambda)^2))
    #   rmse_Counts <- sqrt(mean((data_all[[i]]$X - Lambda)^2))
    #   # Find estimated signatures
    #   id_selected <- res[[i]]$RelWeights > 1.5 * values$epsilon[s]
    #   K_selected <- sum(id_selected)
    #   # Compare true vs estimated
    #   R_true <- apply(data_all[[i]]$Rmat, 2, function(x) x / sum(x))
    #   R_hat <- res[[i]]$Signatures[, id_selected]
    #   matchedSign <- match_MutSign(R_true = R_true, R_hat = R_hat)
    #   cos_sim <- mean(get_cosine_similarity(matchedSign))
    #   # Step 4 - calculate the RMSE between Theta and the rest
    #   Theta_hat <- res[[i]]$Weights[id_selected, ]
    #   rmse_R <- compute_RMSE_Signature(R_hat = matchedSign$R_hat, R_true = matchedSign$R_true)  
    #   rmse_Theta <- compute_RMSE_Theta(Theta_true = data_all[[i]]$Theta, Theta_hat = Theta_hat, matchedSign$match)  
    #   # Step 5 - calculate the sensitivity and precision
    #   sens_prec  <- Compute_sensitivity_precision(R_hat = R_hat, R_true)
    #   rmse_Counts_filer <- sqrt(mean((data_all[[i]]$X - R_hat %*% Theta_hat)^2))
    #   rmse_Lambda_filer <- sqrt(mean((Lambda_true - R_hat %*% Theta_hat)^2))
    #   results <- c("Kest" = K_selected,
    #                "mean" = mean(res[[i]]$RelWeights[!id_selected]), 
    #                "lowCI" = quantile(res[[i]]$RelWeights[!id_selected], probs = c(0.05)), 
    #                "highCI" = quantile(res[[i]]$RelWeights[!id_selected], probs = c(0.95)), 
    #                sens_prec, 
    #                rmse_R, 
    #                rmse_Theta, 
    #                "mean_cos_sim" = cos_sim, 
    #                "rmse_Lambda" = rmse_Lambda, 
    #                "rmse_Counts" = rmse_Counts,
    #                "rmse_Lambda_filter" = rmse_Lambda_filer, 
    #                "rmse_Counts_filter" = rmse_Counts_filer)
    # })
    output<- data.frame(do.call("rbind", output))
    output$epsilon = values$epsilon[s]
    output$Kused = values$K[s]
    output
  }
  #stopCluster(cl)
  return(results)
}

#epsilon_range <- c(0.001, 0.01, 0.1, 0.25, 0.5, 1, 1.5, 2)
epsilon_range <- c(0.001, 0.01, 0.1, 0.25, 1, 2)
K_range <- c(5, 10, 20, 30, 40, 50)


ndatasets <- 40
nsamples <- 500
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
write_csv(results_correct, 
          file = "output/sensitivity_simulation/sensitivity_epsilon_K_100_correct.csv")

set.seed(10, kind = "L'Ecuyer-CMRG")
results_misp <- test_epsilon_k(K_range = K_range,
                               epsilon_range = epsilon_range,
                               J = J,
                               ndatasets = ndatasets,
                               overdispersion = 0.15,
                               nsamples = nsamples, burnin = burnin)
results_misp$overdispersion <- 0.15
write_csv(results_misp, 
          file = "output/sensitivity_simulation/sensitivity_epsilon_K_100_misp.csv")


#write_csv(rbind(results_correct, results_misp), 
#          file = "output/sensitivity_simulation/sensitivity_epsilon_K_100.csv")

df_Keps <- rbind(read_csv("output/sensitivity_simulation/sensitivity_epsilon_K_100_correct.csv"),
                    read_csv("output/sensitivity_simulation/sensitivity_epsilon_K_100_misp.csv"))


#df_Keps <- read_csv("output/sensitivity_simulation/sensitivity_epsilon_K_100.csv")

# # 
# # dodge_eps <- data.frame("Kused" = as.factor(c(1:5) * 10), 
# #                       "dodge" = c(-0.02, -0.01, 0, 0.01, 0.02))
# # 
# # ggplot(df_Keps %>%
# #        mutate(Kused = as.factor(Kused)) %>%
# #        left_join(dodge_eps, by = "Kused")) +
# # geom_abline(slope = 1, intercept = 0, col = "grey") +
# # geom_point(aes(x = epsilon, y = mean, color = Kused), 
# #            position = position_dodge(width = 0.05))+
# # theme_bw()
# # 
# # 
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
# # 
# # # 
# # # 
# # ## RMSE as a function of K and epsilon
pK_rmse <- df_Keps %>%
  mutate(epsilon = as.factor(epsilon),
         overdispersion = paste0("overd = ", overdispersion)) %>%
  group_by(epsilon, Kused, overdispersion) %>%
  summarise(rmse = mean(rmse_Counts)) %>%
  ggplot() +
  geom_point(aes(x = Kused, y = rmse, color = epsilon)) +
  geom_line(aes(x = Kused, y = rmse, color = epsilon)) +
  theme_bw() +
  facet_wrap(~overdispersion, scales = "free")+
  ylab("RMSE for the counts")
# 
# # ## Number of factors as a function of K and epsilon
pK_sig <- df_Keps %>%
  mutate(epsilon = as.factor(epsilon),
         overdispersion = paste0("overd = ", overdispersion)) %>%
  group_by(epsilon, Kused, overdispersion) %>%
  summarise(Kest = mean(Kest)) %>%
  ggplot() +
  geom_hline(yintercept =  6, linetype = "dashed", color = "gray50") +
  geom_point(aes(x = Kused, y = Kest, color = epsilon)) +
  geom_line(aes(x = Kused, y = Kest, color = epsilon)) +
  theme_bw() +
  facet_wrap(~overdispersion, scales = "free")+
  ylab("Estimated n. of signatures")+
  xlab("Value of K in the model")

ggpubr::ggarrange(pK_rmse, pK_sig, common.legend = TRUE)
# 
# # # 
# # # 
# # # 
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

# # 
# # 
# # # 
# # # # Let's now plot the rmse 
# # # ggplot(results_correct %>%
# # #        mutate(Kused = as.factor(Kused))) +
# # #        #geom_abline(slope = 1, intercept = 0, col = "grey") +
# # #        geom_point(aes(x= Kused, y = rmse_Lambda, color = as.factor(epsilon)))+
# # #        geom_line(aes(x= Kused, y = rmse_Lambda, color = as.factor(epsilon)))+
# # #      theme_bw()
# # #      
# # #      
# # #      ggplot(results_correct %>%
# # #               mutate(Kused = as.factor(Kused),
# # #                      epsilon = as.factor(epsilon))) +
# # #        #geom_abline(slope = 1, intercept = 0, col = "grey") +
# # #        geom_boxplot(aes(x= Kused, y = rmse_Signatures, color = epsilon))+
# # #        theme_bw()
# # #      
# # #      
# # #      ggplot(results_misp %>%
# # #               mutate(Kused = as.factor(Kused))) +
# # #        geom_abline(slope = 1, intercept = 0, col = "grey") +
# # #        geom_point(aes(x= epsilon, y = mean, color = Kused))+
# # #        theme_bw()
# # #      
# # #      ggplot(results_correct %>%
# # #               mutate(Kused = as.factor(Kused),
# # #                      epsilon = as.factor(epsilon))) +
# # #        #geom_abline(slope = 1, intercept = 0, col = "grey") +
# # #        geom_boxplot(aes(y= Kest, x= epsilon, color = epsilon))+
# # #        facet_wrap(~Kused, nrow = 1)+
# # #        theme_bw()
# # #      
# # #      
# # # 
# # # 
