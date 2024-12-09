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


simulate_data_sparse <- function(J = 100, 
                                 cosmic_sig = c("SBS1", "SBS2", "SBS5", "SBS13"), 
                                 K_new = 0, 
                                 alpha = 0.25, 
                                 theta = 100, 
                                 overdispersion = 0, 
                                 pi0 = 0.1){
  
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
  
  # Zero-out the weights at random with probability pi0
  is_Theta_zero <- matrix(rbinom(K * J, 1, pi0), ncol = J, nrow = K)
  if(sum(is_Theta_zero) > 0) {
    array_zeros <- which(is_Theta_zero == 1, arr.ind = TRUE)
    for(r in 1:nrow(array_zeros)){
      Theta[array_zeros[r, 1], array_zeros[r, 2]] <- rgamma(1, 0.1, 0.5)
    }
  }
  
  # Generate the counts
  Lambda <- Rmat %*% Theta
  X <- matrix(rnbinom(length(Lambda), size = 1/overdispersion, mu = c(Lambda)), nrow = 96, ncol = J)
  rownames(X) <- rownames(Rmat) 
  
  # Check if all patients have at least one mutation. If not, 
  # substitute with a 1 in the entry has the highest lambda value
  mean_X <- colMeans(X)
  if(any(mean_X == 0)){
    id_X_zero <- which(mean_X == 0)
    for(j in id_X_zero){
      id_max_lambda <- which.max(Lambda[, j])
      X[id_max_lambda, j] <- 1
    }
  }
  
  return(list(X = X, Rmat = Rmat, Theta = Theta))
}

# 
# set.seed(10)
# library(CompressiveNMF)
# data <- simulate_data_sparse(K_new = 2, pi0 = 0.9)
# plot(table(data$X))
# sum(data$X == 0)/(96 * 100)
# sum(colMeans(data$X) > 0)
# 
# out <- CompressiveNMF(X = data$X, K = 20, nsamples = 1000, burnin = 2000)
# print(out)
# plot(out)
# 
# 
# Rmat <- data$Rmat
# Theta <- data$Theta
# X <- data$X
# Lambda <- Rmat %*% Theta
# 
# id_zero_x <- which(colSums(X) == 0)
# Lambda[, id_zero_x]
# 
# 
# # Patients must have at least one active signatur
# sum(colMeans(X) > 0)
# 
# array_zeros <- which(is_Theta_zero == 1, arr.ind = TRUE)
# 
# array_zeros
# 
# # Plot the matrix now based on the sparsity level
# library(pheatmap)
# data <- simulate_data_sparse(K_new = 2, pi0 = 0.75)
# pheatmap(data$X,
#          color=colorRampPalette(brewer.pal(9,"Blues")[c(1,8)])(30),
#          cluster_cols = FALSE,
#          cluster_rows = FALSE,
#          annotation_names_row = FALSE,
#          show_rownames = FALSE,
#          show_colnames = FALSE,
#          legend = FALSE, silent = FALSE,
#          #cellwidth = 1,
#          #cellheight = 1,
#          main = "temp",
#          fontsize = 9,
#          border_color=FALSE,
#          annotation_legend=FALSE,
#          gaps_row = FALSE,
#          gaps_col = FALSE)
# 
# load("data/MutMatrix_IFM2009.rdata")
# MutMatrix_IFM2009$SBS
# pheatmap(MutMatrix_IFM2009$ID,
#          color=colorRampPalette(brewer.pal(9,"Blues")[c(1,8)])(30),
#          cluster_cols = FALSE,
#          cluster_rows = FALSE,
#          annotation_names_row = FALSE,
#          show_rownames = FALSE,
#          show_colnames = FALSE,
#          legend = FALSE, silent = FALSE,
#          #cellwidth = 1,
#          #cellheight = 1,
#          main = "temp",
#          fontsize = 9,
#          border_color=FALSE,
#          annotation_legend=FALSE,
#          gaps_row = FALSE,
#          gaps_col = FALSE)
# 
# data <- simulate_data_sparse(K_new = 2, pi0 = 0.8)
# X_temp <- data$X
# #rownames(X_temp) <- 1:nrow(X_temp)
# longData <- reshape2::melt(MutMatrix_IFM2009$DBS)
# colnames(longData) <- c("Mutation", "Patient", "value")
# longData[longData == 0] <- NA
# ggplot(longData, aes(x = Patient, y = Mutation)) +
#   geom_raster(aes(fill=value), hjust = 0, vjust = 0) +
#   #scale_fill_gradient2(low="yellow",midpoint = 0.5, mid ="orange",
#   #                    high="red3", limits = c(0,1)) +
#   #scale_y_reverse()+
#   theme_minimal()+
#   ylab("Mutation channel") +
#   scale_fill_distiller("Mutations", palette = "YlOrRd", direction = 1)+
#   #theme_minimal() +
#   theme(aspect.ratio = 1,
#         #axis.title=element_blank(),
#         axis.text.y=element_blank(),
#         axis.ticks=element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank()) #+
#   #facet_wrap(.~network, ncol = ncol, nrow = nrow)
# 
# # Simulate the data now!
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
