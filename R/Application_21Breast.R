# This file analyzes the 21breast cancer dataset from NikZainal et al (2016)
# The dataset is available in the package signeR
library(CompressiveNMF)
# Load the packages
library(tidyverse)
library(coda)
source("R/SignatureAnalyzer.R")
source("R/SigProfilerExtractor.R")
source("R/signeR.R")
source("R/PoissonCUSP.R")
source("R/CompressiveNMF.R")
source("R/plot_signatures.R")
source("R/Postprocess_functions.R")
source("R/plot_signatures.R")


# Load the dataset
mut <- t(read.table(system.file("extdata","21_breast_cancers.mutations.txt", package="signeR"), header=TRUE, check.names=FALSE))
transformed_names <- apply(stringr::str_split(rownames(mut), ":", simplify = TRUE), 1, function(x) {
  y <- str_split(x[2], "", simplify = TRUE)
  y[2] <- paste0("[", x[1],"]")
  paste0(y, collapse = "")
})

X <- mut[order(transformed_names), ]
rownames(X) <- sort(transformed_names)

################################################################################
# Part 0. analyze the data
################################################################################
base_subs <- c("C>A", "C>T", "C>G", "T>A", "T>G", "T>C")
mut_all <- sapply(base_subs, function(x) colSums(X[grepl(x, rownames(X)), ]))
mean(colSums(mut_all[rownames(mut_all) != "PD4120a", ]))
sum(mut_all[rownames(mut_all) == "PD4120a", ])

mean(colSums(mut_all[rownames(mut_all) != "PD4120a", ])[1:3])
mean(colSums(mut_all[rownames(mut_all) != "PD4120a", ])[4:6])

################################################################################
# Part 1. Model fitting
################################################################################

rerun <- FALSE
if (rerun) {
  set.seed(10)
  
  #-------------------------------------------------------------------------------
  # SignatureAnalyzer
  #-------------------------------------------------------------------------------
  
  #---------------------------------------------------------------- L1W.L2H
  # This is the version described in the paper. The L1KL version yielded very similar results.
  out_ARD_pcawg <- SignatureAnalyzer(X, method = "L1W.L2H")
  saveRDS(out_ARD_pcawg, file = "output/Application_21brca/ARD_pcawg.rds.gzip", compress = "gzip")
  
  #---------------------------------------------------------------- L1KL
  out_ARD_kl <- SignatureAnalyzer(X, method = "L1KL")
  saveRDS(out_ARD_kl, file = "output/Application_21brca/ARD_kl.rds.gzip", compress = "gzip")
  
  #-------------------------------------------------------------------------------
  # CompressiveNMF
  #-------------------------------------------------------------------------------
  
  #--------------------------- CompressiveNMF
  out_CompNMF <- CompressiveNMF(X = X, K = 15, alpha = 0.5, a = 1, burnin = 10000, epsilon = 0.01,
                                nsamples = 2000, progressbar = FALSE, ncores = 4, nchains = 4)
  print(out_CompNMF)
  saveRDS(out_CompNMF, file = "output/Application_21brca/CompressiveNMF.rds.gzip", compress = "gzip")
  
  # To save space, we remove the chains that are not used
  out_CompNMF_bestchain <- out_CompNMF
  out_CompNMF_bestchain$mcmc_out[-out_CompNMF_bestchain$selected_chain] <- NULL
  out_CompNMF_bestchain$selected_chain <- 1
  saveRDS(out_CompNMF_bestchain, file = "output/Application_21brca/CompressiveNMF_bestchain.rds.gzip", compress = "gzip")
  
  #-------------- CompressiveNMF + cosmic
  out_CompNMF_cosmic_all <- CompressiveNMF(X = X, use_cosmic = TRUE,
                                           K = 10, epsilon = 0.01,
                                           alpha = 0.5, a = 1, burnin = 10000, 
                                           nsamples = 2000, progressbar = FALSE, ncores = 4, nchains = 4, swap_prior = TRUE)
  print(out_CompNMF_cosmic_all)
  saveRDS(out_CompNMF_cosmic_all_nochains, file = "output/Application_21brca/CompressiveNMF_cosmic_all_best.rds.gzip")
  
  out_CompNMF_cosmic_all <- readRDS("output/Application_21brca/CompressiveNMF_cosmic_all.rds.gzip")
  out_CompNMF_cosmic_all
  
  # To save space, we remove the chains that are not used and we remove the redundant signatures
  out_CompNMF_cosmic_all_bestchain <- out_CompNMF_cosmic_all
  out_CompNMF_cosmic_all_bestchain$mcmc_out[-out_CompNMF_cosmic_all_bestchain$selected_chain] <- NULL
  out_CompNMF_cosmic_all_bestchain$selected_chain <- 1
  
  chain <- out_CompNMF_cosmic_all_bestchain$mcmc_out[[1]]
  nonzero_sign <- which(colMeans(chain$Mu) > 0.05)
  out_CompNMF_cosmic_all_bestchain$mcmc_out[[1]]$Signatures <- chain$Signatures[, , nonzero_sign] 
  out_CompNMF_cosmic_all_bestchain$mcmc_out[[1]]$Weights <- chain$Weights[, nonzero_sign , ] 
  out_CompNMF_cosmic_all_bestchain$mcmc_out[[1]]$Mu <- chain$Mu[, nonzero_sign] 
  
  saveRDS(out_CompNMF_cosmic_all_bestchain, file = "output/Application_21brca/CompressiveNMF_cosmic_all_bestchain.rds.gzip", compress = "gzip")
  
  #-------------------------------------------------------------------------------
  # SigneR
  #-------------------------------------------------------------------------------
  out_signeR <- run_signeR(X, Kmin = 2, Kmax = 15, estimate_hyper = TRUE)
  saveRDS(out_signeR, file = "output/Application_21brca/signeR.rds.gzip", compress = "gzip")
  
  #-------------------------------------------------------------------------------
  # PoissonCUSP
  #-------------------------------------------------------------------------------
  set.seed(42)
  out_CUSP <- PoissonCUSP(X, K = 15, nsamples = 2000, burnin = 10000)
  saveRDS(out_CUSP, file = "output/Application_21brca/PoissonCUSP.rds.gzip", compress = "gzip")
  
  #-------------------------------------------------------------------------------
  # SigProfiler
  #-------------------------------------------------------------------------------
  out_sigPro <- sigprofiler(X, "sigprofiler_out/brca21/", Kmin = 2, Kmax = 15, cores = 20)
  saveRDS(out_sigPro, file = "output/Application_21brca/sigPro.rds.gzip", compress = "gzip")
  
}


################################################################################
# Part 2. Analyze the output
################################################################################

# Reload all the models
out_ARD_pcawg <- readRDS("output/Application_21brca/ARD_pcawg.rds.gzip")
out_CompNMF <- readRDS("output/Application_21brca/CompressiveNMF_bestchain.rds.gzip")
out_CompNMF_cosmic_all <- readRDS("output/Application_21brca/CompressiveNMF_cosmic_all_bestchain.rds.gzip")
out_signeR <- readRDS("output/Application_21brca/signeR.rds.gzip")
out_sigPro <- readRDS("output/Application_21brca/sigPro.rds.gzip")
out_CUSP <- readRDS("output/Application_21brca/PoissonCUSP.rds.gzip")


#-------------------------------------------------- Figure 3 - Panel A
df_matrix <- data.frame()

# Compressive NMF
df_temp <- match_to_cosmic_uncertainty_CompNMF(out_CompNMF)
df_temp$avg_theta <- rowMeans(out_CompNMF$Weights)
df_temp$method <- "CompNMF"
df_matrix <- rbind(df_matrix, df_temp)

# Compressive NMF + cosmic
df_temp <- match_to_cosmic_uncertainty_CompNMF(out_CompNMF_cosmic_all)
df_temp$avg_theta <- rowMeans(out_CompNMF_cosmic_all$Weights)
df_temp$method <- "CompNMF + cosmic"
df_matrix <- rbind(df_matrix, df_temp)

# PoissonCUSP
pos_cusp <- get_posterior_CUSP(out_CUSP)
df_temp <- match_to_cosmic_uncertainty_PoissonCUSP(out_CUSP)
df_temp$avg_theta <- rowMeans(pos_cusp$Theta_hat[pos_cusp$nspike < 0.05, ])
df_temp$method <- "PoissonCUSP"
df_matrix <- rbind(df_matrix, df_temp)

# SigProfiler
df_temp <- match_to_cosmic(out_sigPro$Signatures)
df_temp$lowCI <- NA
df_temp$highCI <- NA
df_temp$avg_theta <- rowMeans(out_sigPro$Weights)
df_temp$method <- "SigProfiler"
df_matrix <- rbind(df_matrix, df_temp)

# SigneR
df_temp <- match_to_cosmic_uncertainty_signeR(out_signeR)
df_temp$avg_theta <- rowMeans(out_signeR$Ehat)
df_temp$method <- "signeR"
df_matrix <- rbind(df_matrix, df_temp)

# SigAnalyzer
df_temp <- match_to_cosmic(out_ARD_pcawg$Signature.norm)
df_temp$lowCI <- NA
df_temp$highCI <- NA
df_temp$avg_theta <- rowMeans(out_ARD_pcawg$Exposure)
df_temp$method <- "SignatureAnalyzer"
df_matrix <- rbind(df_matrix, df_temp)

col_values <- c("blue", "skyblue1", "darkorange", "#FFD166", "red", "brown", "#999999")


#------------------------------------------------------------------ # Figure 3 (A)
df_matrix %>%
  mutate(best_cosmic = fct_relevel(best_cosmic,
                                   "SBS1", "SBS2", "SBS3", "SBS8", "SBS9", "SBS13", "SBS34", "SBS40a", "SBS96", "SBS98"), 
         method = fct_relevel(method,
                              "CompNMF + cosmic", "CompNMF", "signeR", "SigProfiler", "SignatureAnalyzer", "PoissonCUSP")) %>%
  ggplot() + 
  geom_point(aes(x = best_cosmic, y = cosine_sim, color = method))+
  geom_errorbar(aes(x = best_cosmic, ymin=lowCI, ymax=highCI, color = method), width = 0.5)+
  coord_cartesian(ylim = c(0.45, 1.0)) +
  scale_fill_manual(values = col_values) +
  scale_color_manual(values = col_values) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, color = "gray35",
                                   vjust = 1,hjust = 1), 
        axis.title.x = element_blank(),
        strip.text = element_text(colour = "black", size = 12),
        legend.position = "none") +
  facet_wrap(~method, ncol = 1) +
  ylab("Cosine Similarity with closest cosmic signature")
ggsave("figures/uncertainty_all_methods.pdf", height = 6.44, width = 2.50)


#------------------------------------------------------------ Figure 3 - Panel B,C
# Plot the signature of CompNMF + cosmic.
sigMat <- out_CompNMF_cosmic_all$Signatures
# Get the posterior samples for all signatures
chain <- out_CompNMF_cosmic_all$mcmc_out[[out_CompNMF_cosmic_all$selected_chain]]
nonzero_sign <- which(colMeans(chain$Mu) > 0.05)
sigChain <- chain$Signatures[, , nonzero_sign]
lowCI <- apply(sigChain, c(2,3), function(x) quantile(x, 0.05))
highCI <- apply(sigChain, c(2,3), function(x) quantile(x, 0.95))
colnames(lowCI) <- colnames(highCI) <- colnames(sigMat)
rownames(sigMat) <- rownames(lowCI) <- rownames(highCI) <- CompressiveNMF::COSMIC_SBS96_GRCh37_v3.4$Channel
pCompNMF_cos <- CompressiveNMF:::plot.SBS.signature(sigMat, lowCI = lowCI, highCI = highCI) + theme(axis.text.x = element_blank())

# Plot the signature of CompNMF.
sigMat <- out_CompNMF$Signatures
# Get the posterior samples for all signatures
chain <- out_CompNMF$mcmc_out[[out_CompNMF$selected_chain]]
nonzero_sign <- which(colMeans(chain$Mu) > 0.05)
df_sim <- match_to_cosmic(sigMat) %>% arrange(best_cosmic)
colnames(sigMat) <- paste0("Sig", c("D", "F", "A", "B", "E", "C"))
sigChain <- chain$Signatures[, , nonzero_sign]
lowCI <- cbind(apply(sigChain, c(2,3), function(x) quantile(x, 0.05)), matrix(0, ncol = 2, nrow = 96))
highCI <- cbind(apply(sigChain, c(2,3), function(x) quantile(x, 0.95)), matrix(0, ncol = 2, nrow = 96))
sigMat <- cbind(sigMat, matrix(0, ncol = 2, nrow = 96))
colnames(lowCI) <- colnames(highCI) <- colnames(sigMat)
rownames(sigMat) <- rownames(lowCI) <- rownames(highCI) <- CompressiveNMF::COSMIC_SBS96_GRCh37_v3.4$Channel
pCompNMF <- CompressiveNMF:::plot.SBS.signature(sigMat, lowCI = lowCI, highCI = highCI)  + theme(axis.text.x = element_blank())

ggpubr::ggarrange(pCompNMF_cos, pCompNMF, nrow = 1)
ggsave("figures/sig_CompNMF_and_CompNMF_cos_21breast.pdf", width = 9.72, height = 6.44)
#--------------------------------------------------------------------- RMSE

# CompNMF
round(sqrt(mean((X - out_CompNMF$Signatures %*% out_CompNMF$Weights)^2)), 2)
# CompNMF + cosmic 
round(sqrt(mean((X - out_CompNMF_cosmic_all$Signatures %*% out_CompNMF_cosmic_all$Weights)^2)), 2)
# signeR
round(sqrt(mean((X - out_signeR$Phat %*% out_signeR$Ehat)^2)), 2)
# PoissonCUSP
round(sqrt(mean((X - pos_cusp$R_hat %*% pos_cusp$Theta_hat)^2)), 2)
# SignatureAnalyzer
round(sqrt(mean((X - out_ARD_pcawg$Signature.norm %*% out_ARD_pcawg$Exposure)^2)), 2)
# SigProfiler
round(sqrt(mean((X - out_sigPro$Signatures %*% out_sigPro$Weights)^2)), 2)


#--------------------------------------------------------------------- 
# Table 8.1
# Effective Sample size for the Bayesian methods
df_ess <- rbind(get_ESS_Comp(out_CompNMF), 
                get_ESS_Comp(out_CompNMF_cosmic_all), 
                get_ESS_signeR(out_signeR, which_samples = 1:2000), 
                get_ESS_PoissonCUSP(out_CUSP))
rownames(df_ess) <- c("CompNMF", "CompNMF+comsic", "signeR", "PoissonCUSP")
round(df_ess, 1)

#----------------------------------------------------- Figures in the supplement

#-------------------------------------------------- CompNMF cosmic all
sigMat <- out_CompNMF_cosmic_all$Signatures
Wmat <- out_CompNMF_cosmic_all$Weights
Wmat <- t(apply(Wmat, 2, function(x) x/sum(x)))
rownames(Wmat) <- colnames(X)
colnames(Wmat)[8] <- "SBSnew1"
plot_weights(Wmat)
ggsave("figures/Weight_compNMF_cosmic.pdf", height = 3.13, width = 5.05)
plot_weights(Wmat)
plot_matrix_signature_v2(sigMat) + theme(aspect.ratio = 0.6) + labs(title = "CompNMF + cosmic")
ggsave("figures/sig_suppl_CompNMF_cosmic.pdf", height = 9.08, width = 9.08)

#-------------------------------------------------- CompNMF 
sigMat <- out_CompNMF$Signatures
Wmat <- out_CompNMF$Weights
Wmat <- t(apply(Wmat, 2, function(x) x/sum(x)))
rownames(Wmat) <- colnames(X)
df_sim <- match_to_cosmic(sigMat) %>% arrange(best_cosmic)
Wmat <- Wmat[, df_sim$signature]
sigMat <- sigMat[, df_sim$signature]
colnames(Wmat) <- paste0("Sig", toupper(letters)[1:ncol(Wmat)], " (", paste(df_sim$best_cosmic, round(df_sim$cosine_sim, 2)), ")")
colnames(sigMat) <- paste0("Sig", toupper(letters)[1:ncol(Wmat)])

# Plot signatures
plot_matrix_signature_v2(sigMat) + theme(aspect.ratio = 0.6) + labs(title = "CompNMF")
ggsave("figures/sig_suppl_CompNMF.pdf", height = 9.08, width = 9.08)
# Plot weights
plot_weights(Wmat, col_palette = c("darkblue", "#2473D0", "#A4C1ED","#F2DAAC", "#E82C36", "grey45"))
ggsave("figures/Weight_compNMF.pdf", height = 3.13, width = 5.05)


#-------------------------------------------------- signeR 
sigMat = out_signeR$Phat
Wmat <- out_signeR$Ehat
Wmat <- t(apply(Wmat, 2, function(x) x/sum(x)))
rownames(Wmat) <- colnames(X)
df_sim <- match_to_cosmic(sigMat) %>% arrange(best_cosmic)
sigMat <- sigMat[, df_sim$signature]
Wmat <- Wmat[, df_sim$signature]
colnames(Wmat) <- paste0("Sig", toupper(letters)[1:ncol(Wmat)], " (", paste(df_sim$best_cosmic, round(df_sim$cosine_sim, 2)), ")")
colnames(sigMat) <- paste0("Sig", toupper(letters)[1:ncol(Wmat)])

# Plot signatures
plot_matrix_signature_v2(sigMat) + theme(aspect.ratio = 0.6) + labs(title = "signeR")
ggsave("figures/sig_suppl_signeR.pdf", height = 9.08, width = 9.08)
# Plot weights
plot_weights(Wmat, col_palette = c("darkblue", "#2473D0", "#A4C1ED", "#F2DAAC", "#A69576"))+
  guides(fill = guide_legend(nrow = 3))
ggsave("figures/Weight_signeR.pdf", height = 3.13, width = 5.05)


#-------------------------------------------------------- PoissonCUSP
sigMat <- pos_cusp$R_hat[, pos_cusp$nspike<0.05]; colnames(sigMat) <- paste0("temp", c(1:ncol(sigMat)))
Wmat <- pos_cusp$Theta_hat[pos_cusp$nspike<0.05, ]; rownames(Wmat) <- paste0("temp", c(1:ncol(sigMat)))
Wmat <- t(apply(Wmat, 2, function(x) x/sum(x)))
rownames(Wmat) <- colnames(X)
df_sim <- match_to_cosmic(sigMat) %>% arrange(best_cosmic)
sigMat <- sigMat[, df_sim$signature]
Wmat <- Wmat[, df_sim$signature]
colnames(Wmat) <- paste0("Sig", toupper(letters)[1:ncol(Wmat)], " (", paste(df_sim$best_cosmic, round(df_sim$cosine_sim, 2)), ")")
colnames(sigMat) <- paste0("Sig", toupper(letters)[1:ncol(Wmat)])

# Plot weights
plot_weights(Wmat, col_palette = c("darkblue", "#2473D0", "#A4C1ED", "#F2DAAC", "red", "darkred", "grey45"))+
  guides(fill = guide_legend(nrow = 3))
ggsave("figures/Weight_PoissonCUSP.pdf", height = 3.13, width = 5.05) 
# Plot signatures
plot_matrix_signature_v2(sigMat) + theme(aspect.ratio = .6) + labs(title = "PoissonCUSP")
ggsave("figures/sig_suppl_PoissonCUSP.pdf", height = 9.08, width = 9.08)

#-------------------------------------------------------- SignatureAnalyzer
sigMat <- out_ARD_pcawg$Signature.norm
Wmat <- out_ARD_pcawg$Exposure
Wmat <- t(apply(Wmat, 2, function(x) x/sum(x)))
rownames(Wmat) <- colnames(X)
df_sim <- match_to_cosmic(sigMat) %>% arrange(best_cosmic)
sigMat <- sigMat[, df_sim$signature]
Wmat <- Wmat[, df_sim$signature]
colnames(Wmat) <- paste0("Sig", toupper(letters)[1:ncol(Wmat)], " (", paste(df_sim$best_cosmic, round(df_sim$cosine_sim, 2)), ")")
colnames(sigMat) <- paste0("Sig", toupper(letters)[1:ncol(Wmat)])

# Plot weights
plot_weights(Wmat, col_palette = c("darkblue", "#2473D0", "#A4C1ED", "#A4C1ED7D", "#F2DAAC"))+
  guides(fill = guide_legend(nrow = 3))
ggsave("figures/Weight_SigAnalyzer.pdf", height = 3.13, width = 5.05)
# Plot signatures
plot_matrix_signature_v2(sigMat) + theme(aspect.ratio = .6) + labs(title = "SignatureAnalyzer")
ggsave("figures/sig_suppl_SigAnalyzer.pdf", height = 9.08, width = 9.08)

#-------------------------------------------------------- SigProfiler
sigMat <- out_sigPro$Signatures
Wmat <- out_sigPro$Weights
Wmat <- t(apply(Wmat, 2, function(x) x/sum(x)))
rownames(Wmat) <- colnames(X)
df_sim <- match_to_cosmic(sigMat) %>% arrange(best_cosmic)
sigMat <- sigMat[, df_sim$signature]
Wmat <- Wmat[, df_sim$signature]
colnames(Wmat) <- paste0("Sig", toupper(letters)[1:ncol(Wmat)], " (", paste(df_sim$best_cosmic, round(df_sim$cosine_sim, 2)), ")")
colnames(sigMat) <- paste0("Sig", toupper(letters)[1:ncol(Wmat)])
# Plot weights
plot_weights(Wmat, col_palette = c("#A4C1ED", "#F2DAAC", "#F29D52")) +
  guides(fill = guide_legend(nrow = 3))
ggsave("figures/Weight_SigProfiler.pdf", height = 3.13, width = 5.05)
# Plot signatures
plot_matrix_signature_v2(sigMat) + theme(aspect.ratio = .6) + labs(title = "SigProfiler")
ggsave("figures/sig_suppl_SigProfiler.pdf", height = 9.08, width = 9.08)



#-------------------------------------------------------- Add BayesNMF from brouwer
dir_files <- "~/CompressiveNMF/output/Application_21brca/BayesNMF_brouwer/"
ReshapeResults_BayesNMF_application <- function(dir_files, J = 21){
  
  # R matrix  
  R <- as.matrix(read_table(paste0(dir_files, "Signatures.txt"), 
                            col_names = FALSE, show_col_types = FALSE))
  R_all <- reshape_python_array(R, dimension = 96)
  Rmean <- apply(R_all, c(2, 3), mean)
  lowCI <- apply(R_all, c(2, 3), function(x) quantile(x, 0.05))
  highCI <- apply(R_all, c(2, 3), function(x) quantile(x, 0.95))
  # Theta
  Theta <- as.matrix(read_table(paste0(dir_files, "Loadings.txt"), 
                                col_names = FALSE, show_col_types = FALSE))
  Theta_all <- reshape_python_array(Theta, dimension = J)
  Theta_mean <- apply(Theta_all, c(2, 3), mean)
  
  # Relevance weights
  Lambda <- as.matrix(read_table(paste0(dir_files, "lambda.txt"), 
                                 col_names = FALSE, show_col_types = FALSE))
  
  # time for execution
  time_exec <- as.matrix(read_table(paste0(dir_files, "times.txt"), 
                                    col_names = FALSE, show_col_types = FALSE))
  time_exec <- time_exec[nrow(time_exec), 1]
  
  # Calculate the effective sample sizes
  EffectiveSigs <- apply(R_all, c(2,3), function(x) coda::effectiveSize(x))
  EffectiveTheta <- apply(Theta_all, c(2,3), function(x) coda::effectiveSize(x))
  EffectiveLambda <- coda::effectiveSize(Lambda)
  
  # Select the number of signatures by excluding the ones that are plainly flat
  norm_weight <- colSums(Rmean)
  sigs_norm <- sapply(1:ncol(Rmean), function(i) Rmean[, i]/norm_weight[i])
  Theta_norm <- t(sapply(1:ncol(Rmean), function(i) Theta_mean[, i]*norm_weight[i]))
  select <- apply(sigs_norm, 2, function(x) cosine(x, rep(1, 96))) < 0.975
  
  lowCI <- apply(lowCI,2,  function(x) x/sum(x))
  highCI <- apply(highCI,2,  function(x) x/sum(x))
  return(list(
    Signatures = sigs_norm,
    CIsigs = list("lowCI" = lowCI, 
                  "highCI" = highCI),
    Theta = Theta_norm, 
    RelWeights = colMeans(Lambda),
    EffectiveSigs = EffectiveSigs,
    EffectiveTheta = EffectiveTheta, 
    EffectiveLambda = EffectiveLambda,
    time = unname(time_exec)
  ))
}













