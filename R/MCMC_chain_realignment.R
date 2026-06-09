# Three examples of label switching
library(CompressiveNMF)
# Load the packages
library(tidyverse)
library(patchwork)
library(coda)
source("~/CompressiveNMF/R/SignatureAnalyzer.R")
source("~/CompressiveNMF/R/SigProfilerExtractor.R")
source("~/CompressiveNMF/R/signeR.R")
source("~/CompressiveNMF/R/PoissonCUSP.R")
source("~/CompressiveNMF/R/CompressiveNMF.R")
source("~/CompressiveNMF/R/plot_signatures.R")
source("~/CompressiveNMF/R/Postprocess_functions.R")
source("~/CompressiveNMF/R/plot_signatures.R")


#-----------------------------------------------------------------------------
detect_label_switch <- function(out_CompNMF, X, 
                                     percSamples = 0.01, 
                                     mu_threshold = 0.005, 
                                     post_start = 100,
                                burnin_line = NULL) {
  
  # 1. Extract MCMC chain and dimensions
  mcmc <- out_CompNMF$mcmc_out[[out_CompNMF$selected_chain]]
  n_iter <- nrow(mcmc$Signatures)
  final_samples <- round(percSamples * n_iter)
  chain_idx <- 1:(n_iter - final_samples)
  
  # 2. Estimate Mu and identify active signatures
  MuEst <- colMeans(mcmc$Mu[-chain_idx, , drop = FALSE])
  id_sigs <- which(MuEst > mu_threshold)
  
  if (length(id_sigs) == 0) stop("No signatures passed the Mu threshold.")
  
  # 3. Estimate Signatures
  SigsEst <- apply(mcmc$Signatures[-chain_idx, , id_sigs, drop = FALSE], c(2, 3), mean)
  rownames(SigsEst) <- rownames(X)
  colnames(SigsEst) <- paste0("Sig", sprintf("%02d", seq_along(id_sigs)))
  
  # 4. Calculate Cosine Similarities (compact lapply replaces the for-loop)
  df_res <- lapply(seq_along(id_sigs), function(i) {
    sims <- cosine(t(mcmc$Signatures[, , id_sigs[i]]), as.matrix(SigsEst))
    colnames(sims) <- colnames(SigsEst)
    
    as.data.frame(sims) %>%
      mutate(iteration = row_number(),
             Signature = paste0("Sig", sprintf("%02d", i))) %>%
      pivot_longer(cols = -c(iteration, Signature), names_to = "Label")
  }) %>% bind_rows()
  
  # 5. Build Plot 1 (Cosine Similarities)
  p1 <- ggplot(df_res, aes(x = iteration, y = value, color = Label)) +
    geom_line() +
    facet_wrap(~ Signature) + 
    ylim(c(0, 1)) + 
    theme_bw() +
    ylab("Cosine sim. with estimate") +
    geom_vline(xintercept = n_iter - final_samples, linetype = "dashed", color = "red")
  if(!is.null(burnin_line)){
    p1 <- p1 + geom_vline(xintercept = burnin_line, linetype = "dashed", color = "blue")
  }
  # 6. Build Plot 2 (Logposterior)
  # subsetting dynamically instead of hardcoding 1000:20000
  idx_post <- post_start:n_iter 
  df_post <- data.frame(iteration = idx_post, logpost = mcmc$logposterior[idx_post])
  
  p2 <- ggplot(df_post, aes(x = iteration, y = logpost)) +
    geom_line() +
    theme_bw() +
    geom_vline(xintercept = n_iter - final_samples, linetype = "dashed", color = "red") +
    labs(x = "MCMC iteration", y = "Logposterior")
  if(!is.null(burnin_line)){
    p2 <- p2 + geom_vline(xintercept = burnin_line, linetype = "dashed", color = "blue")
  }
  
  # 7. Assemble final layout 
  # Note: plot_SBS_signature must be loaded in your environment
  final_plot <- p1 + (plot_SBS_signature(SigsEst) / p2 + plot_layout(heights = c(2, 1)))
  
  # Return both the plot and the extracted signatures
  return(list(
    plot = final_plot,
    SigsEst = SigsEst
  ))
}

plot_SBS_signature2 <- function(signatures,
         lowCI = NULL,
         highCI = NULL,
         palette = c("#40BDEE", "#020202", "#E52925", "#CCC9CA", "#A3CF62", "#ECC5C5")) {
  #load("~/CompressiveNMF/data/Cosmic_data.rdata")
  signatures <- as.matrix(signatures)
  names_sig <- rownames(signatures)
  df_plot <- data.frame(signatures) %>%
    dplyr::mutate(Channel = names_sig,
                  Triplet = apply(stringr::str_split(names_sig, "", simplify = TRUE), 1,
                                  function(x) paste0(x[c(1,3,7)], collapse = "")),
                  Mutation = apply(stringr::str_split(names_sig, "", simplify = TRUE), 1,
                                   function(x) paste0(x[c(3,4,5)], collapse = "")),
                  Mutation = as.factor(Mutation)) %>%
    tidyr::gather(key = "Sig", value = "Prob", -Channel, -Triplet, -Mutation) %>%
    dplyr::mutate(Sig = factor(Sig, levels = colnames(signatures)))
  
  if(!is.null(lowCI) & !is.null(highCI)){
    df_plot <- df_plot %>%
      dplyr::left_join(data.frame(lowCI) %>%
                         dplyr::mutate(Channel = names_sig,
                                       Triplet = apply(stringr::str_split(names_sig, "", simplify = TRUE), 1,
                                                       function(x) paste0(x[c(1,3,7)], collapse = "")),
                                       Mutation = apply(stringr::str_split(names_sig, "", simplify = TRUE), 1,
                                                        function(x) paste0(x[c(3,4,5)], collapse = "")),
                                       Mutation = as.factor(Mutation)) %>%
                         tidyr::gather(key = "Sig", value = "lowCI", -Channel, -Triplet, -Mutation),
                       by = c("Channel", "Triplet", "Mutation", "Sig")) %>%
      dplyr::left_join(data.frame(highCI) %>%
                         dplyr::mutate(Channel = names_sig,
                                       Triplet = apply(stringr::str_split(names_sig, "", simplify = TRUE), 1,
                                                       function(x) paste0(x[c(1,3,7)], collapse = "")),
                                       Mutation = apply(stringr::str_split(names_sig, "", simplify = TRUE), 1,
                                                        function(x) paste0(x[c(3,4,5)], collapse = "")),
                                       Mutation = as.factor(Mutation)) %>%
                         tidyr::gather(key = "Sig", value = "highCI", -Channel, -Triplet, -Mutation),
                       by = c("Channel", "Triplet", "Mutation", "Sig"))
  }
  
  p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = Triplet, y = Prob, fill = Mutation))+
    ggplot2::geom_bar(stat = "identity", width = 0.7) +
    ggplot2::facet_grid(Sig~Mutation, scales = "free")+
    ggplot2::theme_minimal()+
    ggplot2::scale_fill_manual(values = palette)+
    ggplot2::theme(
      legend.position = "none",
      axis.title = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 90, color = "gray35",
                                          vjust = .5, size = 6.5, margin = ggplot2::margin(t = -4)),
      panel.grid = ggplot2::element_blank(),
      panel.spacing.x = ggplot2::unit(0, "lines"),
      panel.spacing.y = ggplot2::unit(0,"lines"))
  
  if(!is.null(lowCI) & !is.null(highCI)){
    p <- p +
      ggplot2::geom_linerange(ggplot2::aes(x = Triplet, ymin =lowCI, ymax = highCI), color = "grey65")
  }
  
  return(p)
  
}

#-----------------------------------------------------------------------------

#==============================================================
# Example 1 - Simulated data, misspecified setting
#==============================================================
data_sim <- readRDS("~/CompressiveNMF/output/main_simulation/Scenario_100_overdisp_0.15_Knew_2_theta_100/data.rds.gzip")
Xsim <- data_sim[[1]]$X

set.seed(42)
out_Sim <- CompressiveNMF(X = Xsim, K = 20, epsilon = 0.001,
                          alpha = 0.5, 
                          a = 1, burnin = 0,
                          nsamples = 10000,
                          progressbar = TRUE,
                          ncores = 1,
                          nchains = 1)

pSimulation <- plot_compNMF_diagnostics(out_Sim, Xsim, mid_dashed_line = 5000)

#==============================================================
# Example 2 - 21 breast cancer application
#==============================================================

#----------------------------------------------------- Load the dataset
mut <- t(read.table(system.file("extdata","21_breast_cancers.mutations.txt", package="signeR"), header=TRUE, check.names=FALSE))
transformed_names <- apply(stringr::str_split(rownames(mut), ":", simplify = TRUE), 1, function(x) {
  y <- str_split(x[2], "", simplify = TRUE)
  y[2] <- paste0("[", x[1],"]")
  paste0(y, collapse = "")
})

X <- mut[order(transformed_names), ]
rownames(X) <- sort(transformed_names)


#----------------------------------------------------- Informative prior
set.seed(42)
out_CompNMF_cosmic <- CompressiveNMF(X = X, use_cosmic = TRUE,
                                     K = 10, epsilon = 0.001,
                                     alpha = 0.5, a = 1, burnin = 0,
                                     nsamples = 10000, progressbar = TRUE,
                                     ncores = 1, nchains = 1, swap_prior = FALSE)
dimnames(out_CompNMF_cosmic$mcmc_out[[1]]$Signatures) <- list(NULL, rownames(X), colnames(out_CompNMF_cosmic$mcmc_out[[1]]$init_pars$R))


#----------------------------------------------------- Fully unsupervised
set.seed(42)
out_CompNMF <- CompressiveNMF(X = X, use_cosmic = FALSE,
                                     K = 20, epsilon = 0.001,
                                     alpha = 0.5, a = 1, burnin = 0,
                                     nsamples = 20000, progressbar = TRUE,
                                     ncores = 1, nchains = 1, swap_prior = FALSE)

p21Breast <- detect_label_switch(out_CompNMF, X, burnin_line = 10000)
ggsave(plot = p21Breast$plot, filename = "~/CompressiveNMF/figures/Plot_LabelSwitching_21Breast.pdf", 
       width = 13.10, height = 4.86)

#==============================================================
# Prior swap example in the two cases (careful elicitation vs vague prior)
#==============================================================

data_sim2 <- readRDS("~/CompressiveNMF/output/main_simulation/Scenario_100_overdisp_0.15_Knew_6_theta_100/data.rds.gzip")
Xsim2 <- data_sim2[[1]]$X

# Informative prior without elicitation
sigs <- c("SBS1", "SBS2", "SBS13", "SBS3", "SBS5", "SBS8", "SBS34", "SBS40a")
S <- CompressiveNMF::COSMIC_v3.4_SBS96_GRCh37[, sigs]

betah_notune <- rep(100, ncol(S))
names(betah_notune) <- colnames(S)

#---- No tune: the prior is not tuned carefully, and we get a mismatch
set.seed(10)
out_CompNMF_notune <- CompressiveNMF(X = Xsim2, S = S, betah = betah_notune,
                                     K = 10,  epsilon = 0.001,
                                     alpha = 0.5, a = 1, burnin = 2000,
                                     nsamples = 1000, progressbar = TRUE,
                                     ncores = 1, nchains = 1, swap_prior = FALSE)
Sigs_notune <- out_CompNMF_notune$Signatures
rownames(Sigs_notune) <- rownames(Xsim2)
plot_SBS_signature(Sigs_notune)
match_to_cosmic(Sigs_notune)

Lambda_notune <- out_CompNMF_notune$Signatures %*% out_CompNMF_notune$Weights
sqrt(mean((Xsim2 - Lambda_notune)^2))


#---- Tune: the prior is tuned in an informative manner. This leads to better results
betah_tune <- CompressiveNMF::Betah_SBS96_GRCh37_v3.4[colnames(S)]

set.seed(10)
out_CompNMF_tune <- CompressiveNMF(X = Xsim2, S = S, betah = betah_tune,
                                     K = 10,  epsilon = 0.001,
                                     alpha = 0.5, a = 1, burnin = 2000,
                                     nsamples = 1000, progressbar = TRUE,
                                     ncores = 1, nchains = 1, swap_prior = FALSE)
Sigs_tune <- out_CompNMF_tune$Signatures
rownames(Sigs_tune) <- rownames(Xsim2)
plot_SBS_signature(Sigs_tune)
match_to_cosmic(Sigs_tune)

Lambda_tune <- out_CompNMF_tune$Signatures %*% out_CompNMF_tune$Weights
sqrt(mean((Xsim2 - Lambda_tune)^2))

#---- Swap: the prior is tuned in an informative manner and we perform a swap at 2/3 of burnin. 
set.seed(10)
out_CompNMF_swap <- CompressiveNMF(X = Xsim2, S = S, betah = betah_tune,
                                   K = 10,  epsilon = 0.001,
                                   alpha = 0.5, a = 1, burnin = 2000,
                                   nsamples = 1000, progressbar = TRUE,
                                   ncores = 1, nchains = 1, swap_prior = TRUE)
Sigs_swap <- out_CompNMF_swap$Signatures
rownames(Sigs_swap) <- rownames(Xsim2)
plot_SBS_signature(Sigs_swap)
match_to_cosmic(Sigs_swap)

Lambda_swap <- out_CompNMF_swap$Signatures %*% out_CompNMF_swap$Weights
sqrt(mean((Xsim2 - Lambda_swap)^2))


#---- Make the final plot now 
match_swap_notune <- RcppHungarian::HungarianSolver(1 - cosine(Sigs_swap, Sigs_notune))
match_swap_tune <- RcppHungarian::HungarianSolver(1 - cosine(Sigs_swap, Sigs_tune))

# Final plot of the signatures
plot_SBS_signature2(Sigs_notune[, match_swap_notune$pairs[, 2]])
plot_SBS_signature2(Sigs_tune[, match_swap_tune$pairs[, 2]])
plot_SBS_signature2(Sigs_swap)



