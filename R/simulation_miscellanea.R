data <- simulate_data(K_new = 5, K_cosmic = 5, J = 100, signal = 50, sd_perturb = 0.005,
                      alpha_new = 0.1)
j <- 2
plot_cosmic_signature(sign = colnames(data$Rmat)[j], sign_inferred = data$Rmat[, j])
k <- 7
plot_signature(data$Rmat[, k])


find_cosine_diff <- function(alpha_new = 0.1, CS = NULL){
  # Generate the signature
  Rnew <- rgamma(n = 96, alpha_new, 0.5)
  Rnew <- Rnew / sum(Rnew)
  # Calculate difference with whole cosmic database
  if(is.null(CS)){
    load("data/Cosmic_data.rdata")
    CS <- cosmic_data[, -c(1,2,3)]
  }
  quantile(apply(CS, 2, function(x) cosine(x, Rnew)))
}

find_cosine_diff(0.1)



# simulate_data <- function(K_new = NULL, K_cosmic = NULL, J = 100, signal = 100,
#                           perturb_cosmic = TRUE, prob_perturbation = 0.05,
#                           size_perturbation = 0.005, alpha_new = 0.5) {
#   I <- 96
#   Rmat <- NULL
#   # Add a perturbed version of the cosmic signatures
#   if (!is.null(K_cosmic)) {
#     load("data/Cosmic_data.rdata")
#     signature_names <- colnames(cosmic_data)[-c(1, 2, 3)]
#     sampled_sign <- sample(signature_names, size = K_cosmic)
#     Rcosmic <- as.matrix(cosmic_data[, sampled_sign])
#     if (perturb_cosmic) {
#       perturbations <- matrix(size_perturbation * rbinom(I * K_cosmic, size = 1, prob = prob_perturbation), nrow = I)
#       Rcosmic <- apply(Rcosmic + perturbations, 2, function(x) x / sum(x))
#     }
#     Rmat <- cbind(Rmat, Rcosmic)
#   }
#   # Add additional extra signatures
#   if (!is.null(K_new)) {
#     Rnew <- matrix(rgamma(n = I * K_new, alpha_new, 0.5), nrow = I, ncol = K_new)
#     Rnew <- apply(Rnew, 2, function(x) x / sum(x))
#     colnames(Rnew) <- paste0("SBSnew", 1:K_new)
#     Rmat <- cbind(Rmat, Rnew)
#   }
#   # Sample the weights
#   Theta <- signal * matrix(rgamma(n = ncol(Rmat) * J, 1, 0.5), nrow = ncol(Rmat), ncol = J)
#   # Generate the counts
#   Lambda <- Rmat %*% Theta
#   X <- matrix(rpois(n = length(Lambda), c(Lambda)), nrow = I, ncol = J)
#   return(list(X = X, Rmat = Rmat, Theta = Theta))
# }


#---------------------------------- 6 - SparseSignatures (?)
# data <- data_all[[1]]
# x <- t(data$X)
# ss_try <- SparseSignatures::nmfLasso(x = x, K = 8, normalize_counts = FALSE, nmf_runs = 5,
#                                      lambda_rate_alpha = 0, lambda_rate_beta = 0.01)
# 
# ss_cross <- SparseSignatures::nmfLassoCV(x = x, K = 5:10, normalize_counts = FALSE,
#                                          nmf_runs = 5, lambda_values_alpha = 0,
#                                          cross_validation_repetitions = 10, iterations = 20)
# 
# cv_out <- ss_cross
# 
# cv_mses <- cv_out$grid_search_mse[1, , ]
# cv_means_mse <- matrix(
#   sapply(cv_mses, FUN = mean),
#   nrow = dim(cv_mses)[1]
# )
# dimnames(cv_means_mse) <- dimnames(cv_mses)
# min_ii <- which(cv_means_mse == min(cv_means_mse), arr.ind = TRUE)
# min_Lambda <- rownames(cv_means_mse)[min_ii[1]]
# min_K <- colnames(cv_means_mse)[min_ii[2]]
# cat("Minimum MSE at:", min_Lambda, "and", min_K, "\n")
# 
# 
# # Check lambda values
# lambda_test_values <- c(0.01, 0.02, 0.05, 0.1, 0.2, 0.225, 0.25, 0.275, 0.3)
# results <- lambdaRangeBetaEvaluation(
#   x = x,
#   K = 10, normalize_counts = FALSE, iterations = 35, nmf_runs = 5,
#   lambda_values = lambda_test_values,
#   num_processes = 20
# )
# 
# for (i in 1:length(lambda_test_values)) {
#   print(colnames(results)[[i]])
#   print(results[[i]]$loglik_progression)
#   plot(results[[i]]$loglik_progression)
# }
# 
# plot(c(-14367716, -56298111, -56342555, -59080222, rep(-59080222, 10)))
# 
# 
# ss_match <- match_MutSign(R_true = data$Rmat, R_hat = t(ss_try$beta))
# j <- 1
# plot_signature(signature = ss_match$R_true[, j], sign_inferred = ss_match$R_hat[,j])
# 
# 
# Lambda_true <-  data$Rmat %*% data$Theta
# Lambda_ss <- t(ss_try$alpha %*% ss_try$beta)
# 
# plot(Lambda_ss, Lambda_true)

#---------------------------------- 7 - SigProfiler - the package of Alexandrov et al. (?)


library(SigProfilerExtractorR)

path_to_example_data <- importdata("matrix")
data <- path_to_example_data # here you can provide the path of your own data

sigprofilerextractor("matrix", 
                     "example_output", 
                     data, minimum_signatures = 1,
                     maximum_signatures=4,
                     nmf_replicates=5,
                     min_nmf_iterations = 1000,
                     max_nmf_iterations =100000,
                     nmf_test_conv = 1000,
                     nmf_tolerance = 1e-08)
sigprofilerextractor

xx <- read_table("/home/alezito/.local/share/r-miniconda/envs/r-reticulate/lib/python3.6/site-packages/SigProfilerExtractor/data/Samples.txt")
sigpro2("matrix",  "example_output", data) 

sigpro2 <- function (input_type, output, input_data, reference_genome = "GRCh37", 
                     opportunity_genome = "GRCh37", cosmic_version = 3.3, context_type = "default", 
                     exome = F, minimum_signatures = , maximum_signatures = 25, 
                     nmf_replicates = 100, resample = T, batch_size = 1, cpu = -1, 
                     gpu = F, nmf_init = "random", precision = "single", matrix_normalization = "gmm", 
                     seeds = "random", min_nmf_iterations = 10000, max_nmf_iterations = 1e+06, 
                     nmf_test_conv = 10000, nmf_tolerance = 1e-15, nnls_add_penalty = 0.05, 
                     nnls_remove_penalty = 0.01, initial_remove_penalty = 0.05, 
                     collapse_to_SBS96 = T, clustering_distance = "cosine", export_probabilities = T, 
                     make_decomposition_plots = T, stability = 0.8, min_stability = 0.2, 
                     combined_stability = 1, allow_stability_drop = F, get_all_signature_matrices = F) 
{
  sys <- reticulate::import("sys")
  sigpro <- reticulate::import("SigProfilerExtractor.sigpro")
  minimum_signatures = as.integer(minimum_signatures)
  maximum_signatures = as.integer(maximum_signatures)
  nmf_replicates = as.integer(nmf_replicates)
  min_nmf_iterations = as.integer(min_nmf_iterations)
  max_nmf_iterations = as.integer(max_nmf_iterations)
  nmf_tolerance = as.numeric(nmf_tolerance)
  nmf_test_conv = as.integer(nmf_test_conv)
  nnls_add_penalty = as.numeric(nnls_add_penalty)
  nnls_remove_penalty = as.numeric(nnls_remove_penalty)
  initial_remove_penalty = as.numeric(initial_remove_penalty)
  stability = as.numeric(stability)
  min_stability = as.numeric(min_stability)
  combined_stability = as.numeric(combined_stability)
  batch_size = as.integer(batch_size)
  cpu = as.integer(cpu)
  cosmic_version = as.numeric(cosmic_version)
  sigpro$sigProfilerExtractor(input_type, output, input_data, 
                              reference_genome = reference_genome, opportunity_genome = opportunity_genome, 
                              cosmic_version = cosmic_version, context_type = context_type, 
                              exome = exome, minimum_signatures = minimum_signatures, 
                              maximum_signatures = maximum_signatures, nmf_replicates = nmf_replicates, 
                              resample = resample, batch_size = batch_size, cpu = cpu, 
                              gpu = gpu, nmf_init = nmf_init, precision = precision, 
                              matrix_normalization = matrix_normalization, seeds = seeds, 
                              min_nmf_iterations = min_nmf_iterations, max_nmf_iterations = max_nmf_iterations, 
                              nmf_test_conv = nmf_test_conv, nmf_tolerance = nmf_tolerance, 
                              nnls_add_penalty = nnls_add_penalty, nnls_remove_penalty = nnls_remove_penalty, 
                              initial_remove_penalty = initial_remove_penalty,  
                              clustering_distance = clustering_distance, export_probabilities = export_probabilities, 
                              make_decomposition_plots = make_decomposition_plots, 
                              stability = stability, min_stability = min_stability, 
                              combined_stability = combined_stability,
                              get_all_signature_matrices = get_all_signature_matrices)
  sys$stdout$flush()
}

out_CompressiveNMF_cosmic <- readRDS("output/Simulation1/CompressiveNMF_cosmic.rds.gzip")
df_res <- aggregate_results(out_CompressiveNMF = out_CompressiveNMF, 
                            out_CompressiveNMF_cosmic = out_CompressiveNMF_cosmic, out_PoissonCUSP = out_PoissonCUSP)
df_res <- aggregate_results(out_CompressiveNMF = out_CompressiveNMF,
                            out_CompressiveNMF_cosmic = out_CompressiveNMF_cosmic,
                            out_PoissonCUSP = out_PoissonCUSP, 
                            out_ARD = out_ARD,
                            out_signeR = out_signeR)


df_res[df_res$Method == "SignatureAnalyzer", ]$K
df_res[df_res$Method == "CompressiveNMF_cosmic", ]$K
df_res[df_res$Method == "CompressiveNMF", ]$K
df_res[df_res$Method == "PoissonCUSP", ]$K
df_res[df_res$Method == "signeR", ]$K


df_res$diffK <- abs(10 - df_res$K)
df_res %>%
  group_by(Method) %>%
  summarise_all(mean)


ggplot(df_res) +
  geom_boxplot(aes(x = Method, y = cos_sim))

ggplot(df_res) +
  geom_boxplot(aes(x = Method, y = rmse_Lambda))

ggplot(df_res) +
  geom_boxplot(aes(x = Method, y = rmse_Counts))

# Did we retrieve the correct signatures?
s <- 2
match <- out_CompressiveNMF_cosmic[[s]]$signatures$match
mu <- out_CompressiveNMF_cosmic[[s]]$Mu_hat
names(mu) <- colnames(S)
(true_sign <- colnames(out_CompressiveNMF_cosmic[[s]]$signatures$R_true))
(inferred_sign <- names(mu)[mu > 0.025][match])

rbind(true_sign, inferred_sign)

j <- 1
plot_cosmic_signature(true_sign[j], sign_inferred = out_CompressiveNMF_cosmic[[s]]$signatures$R_hat[, j])

j <- 3
p1 <- plot_signature(signature = out_CompressiveNMF_cosmic[[s]]$signatures$R_true[, j], sign_inferred = out_CompressiveNMF_cosmic[[s]]$signatures$R_hat[, j])

p2 <- plot_cosmic_signature("SBS93",  sign_inferred = out_CompressiveNMF_cosmic[[s]]$signatures$R_hat[, j])

ggarrange(p1,p2, ncol = 1)
plot_cosmic_signature("SBS5")
plot_cosmic_signature("SBS41")


id <- sample(1:ncol(S), 2)
cosine(S[ , id[1]], S[ , id[2]])


j <- 1
plot_signature(signature = out_CompressiveNMF_cosmic[[s]]$signatures$R_true[, j], sign_inferred = out_CompressiveNMF_cosmic[[s]]$signatures$R_hat[, j])







#
#
#
#
#
# j <- 5
# plot_cosmic_signature(sign = colnames(post_ard$signatures$R_true)[j],
#                       sign_inferred = post_ard$signatures$R_hat[, j])
#
#
#
#
#
# # SignatureAnalyzer
# set.seed(40)
# data <- simulate_mutations_from_cosmic(signal = 50, K0 = 8, J = 50)
# ard <- NMF_l1_ARD(V = data$X, K = 20, normalize = TRUE)
# (post_ard <- Postprocess_ARD(resARD =  ard, data = data))
#
#
#
# out1 <- CompressiveNMF(X = data$X, K = 20, nsamples = 1000, burnin = 500, alpha = 0.4)
# (post <- Postprocess_Compressive(resComp = out1, data = data, fast = FALSE, cutoff = 0.025))
#
# out2 <- CompressiveNMF(X = data$X, S = data$Rmat, K = 10, nsamples = 1000, burnin = 500, alpha = 0.4)
# (post2 <- Postprocess_Compressive(resComp = out2, data = data, fast = FALSE, cutoff = 0.025))
#
#
# S <- as.matrix(cosmic_data[, -c(1,2,3)])
# out3 <- CompressiveNMF(X = data$X, S = S, K = 3, nsamples = 500, burnin = 250, alpha = 0.4)
# (post3 <- Postprocess_Compressive(resComp = out3, data = data, fast = FALSE, cutoff = 0.025))
#
#
# j <- 8
# plot_cosmic_signature(sign = colnames(post3$signatures$R_true)[j], sign_inferred = post3$signatures$R_hat[, j])
#
#
# res_mu <- colMeans(out3$Mu)
# names(res_mu)[1:ncol(S)] <- colnames(S)
# res_mu[res_mu > 0.025]
# sort(colnames(data$Rmat))
# sort(names(res_mu[res_mu > 0.025]))
#
# post_ard$results
# post$results
# post3$results
#
# colMeans(out1$Mu)
#
#
# plot_cosmic_signature("SBS3")
# plot_cosmic_signature("SBS40c")
#
#
# # Simulation
# j <- 8
# plot_cosmic_signature(sign = colnames(post_ard$signatures$R_true)[j], sign_inferred = post_ard$signatures$R_hat[, j])
#
# p_list <- list()
# for(i in 1:16){
#   if(i%%2 == 1){
#     p_list[[i]] <- plot_cosmic_signature(sign = colnames(data$Rmat)[round((i + 1)/2)], add_aetiology = FALSE)
#   } else {
#     p_list[[i]] <- plot_signature(signature = ard$W[, round(i/2)]) #plot_cosmic_signature(sign = colnames(data$Rmat)[i], add_aetiology = FALSE)
#   }
# }
#
# for(i in 9:16){
#   p_list[[i]] <- plot_signature(signature = ard$W[, i-8])
# }
#
# pdf(file = "Plot.pdf", width = 20, height = 26)
# ggarrange(plotlist = p_list, ncol = 2, nrow = 8, legend = "none")
# dev.off()
#
#
#
#
#
# out1 <- CompressiveNMF(X = data$X, K = 20, nsamples = 2000, burnin = 1000, alpha = 0.25)
# post <- Postprocess_Compressive(resComp = out1, data = data, fast = FALSE)
#
#
#
#
#
#
#
# library(tictoc)
#
#
# tic()
# out1 <- CompressiveNMF(X = data$X, K = 15, nsamples = 2000, burnin = 1000, alpha = 0.4)
# toc()
#
#
# tic()
# out2 <- CompressiveNMF_fast(X = data$X, K = 20, nsamples = 3000, burnin = 1000, alpha = 0.4)
# toc()
#
# sort(round(colMeans(out1$Mu), 5))
# sort(round(colMeans(out2$Mu), 5))
# Postprocess_Compressive(out2, data, fast = TRUE)
# post <- Postprocess_Compressive(resComp = out1, data = data, fast = FALSE)
# #post <- Postprocess_Compressive(resComp = out2, data = data, fast = TRUE)
#
#
# ard <- NMF_l1_ARD(V = data$X, K = 20, normalize = TRUE)
# post_ard <- Postprocess_ARD(resARD =  ard, data = data)
#
#
# j <- 3
# plot_cosmic_signature(colnames(post_ard$signatures$R_true)[j], sign_inferred = post_ard$signatures$R_hat[, j])
#
#
#
# # Use signeR now
#
# out_SignR <- signeR(M = data$X,  samples = "col", nlim = c(3,12))
#
# post_signeR <- match_MutSign(R_true = data$Rmat, R_hat = out_SignR$Phat)
#
#
# #plot_cosmic_signature(colnames(post$signatures$R_true)[j], sign_inferred = sign_inferred)
# p1 <- plot_cosmic_signature(colnames(post$signatures$R_true)[j], sign_inferred = post$signatures$R_hat[, j])
# p2 <- plot_cosmic_signature(colnames(post$signatures$R_true)[j], sign_inferred = post_ard$signatures$R_hat[, j])
# p3 <- plot_cosmic_signature(colnames(post_signeR$signatures$R_true)[j], sign_inferred = post_signeR$R_hat[, j])
# ggarrange(p1, p2, p3, ncol = 1)
#
#
#
#
# M <- data$X
# S <- matrix(rgamma(96 * 3, shape = 1, rate = 1), 96)
# S <- apply(S, 2, function(x) x/sum(x))
# out_SignR <- signeR(M=M,  P = S, fixedP = TRUE, samples = "col", nsig = 8)
#
#
#
#
#
# out_SignR$Phat
#
#
# library(SigProfilerExtractorR)
#
#
# out500 <- readRDS("output/Simulation_500.rds.gzip")
# path_to_example_table = importdata("matrix")
# data = path_to_example_table
# install("GRCh37", rsync=FALSE, bash=TRUE)
# read_table(data)
#
#
#
#
#
#

# Test why SignatureAnalyzer takes so long
data <- data_all[[27]]
data$X

set.seed(20, kind = "L'Ecuyer-CMRG")
try <- NMF_l1_ARD(V = data$X, K = 25)














