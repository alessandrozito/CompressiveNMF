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



#----- Old code


dimnames(out_CompNMF$mcmc_out[[1]]$Signatures) <- list(NULL, rownames(X), colnames(out_CompNMF$mcmc_out[[1]]$init_pars$R))

burnin <- 1:(nrow(out_CompNMF$mcmc_out[[1]]$Signatures) - 250)
MuEst <- colMeans(out_CompNMF$mcmc_out[[1]]$Mu[-burnin, ])
id_sigs <- which(MuEst > 0.005)

SigsEst <- apply(out_CompNMF$mcmc_out[[1]]$Signatures[-burnin, , ], c(2, 3), mean)[, id_sigs]
rownames(SigsEst) <- rownames(X)
colnames(SigsEst) <- paste0("Sig", sprintf("%02d", 1:ncol(SigsEst)))

df_res <- data.frame()
for(i in 1:length(id_sigs)){
  res <- cosine(t(out_CompNMF$mcmc_out[[1]]$Signatures[, , id_sigs[i]]), 
                as.matrix(SigsEst))
  colnames(res) <- colnames(SigsEst)
  df_tmp <- as.data.frame(res) %>%
    mutate(iteration = row_number()) %>%
    pivot_longer(cols = -iteration)
  df_tmp$Signature <- paste0("Sig", sprintf("%02d", i))
  df_tmp$Label <- df_tmp$name
  df_res <- rbind(df_res, df_tmp)
}

p1 <- ggplot(df_res, aes(x = iteration, y = value, color = Label))+
  geom_line() +
  facet_wrap(~ Signature) + 
  ylim(c(0, 1)) + 
  theme_bw() +
  ylab("Cosine sim. with estimate") +
  geom_vline(xintercept = 10000, linetype = "dashed") +
  geom_vline(xintercept = 20000 - 250, linetype = "dashed", color = "red")

p2 <- ggplot() +
  geom_line(aes(x = 1000:20000, y = out_CompNMF$mcmc_out[[1]]$logposterior[1000:20000])) +
  theme_bw() +
  geom_vline(xintercept = 10000, linetype = "dashed") +
  geom_vline(xintercept = 20000 - 250, linetype = "dashed", color = "red") +
  ylab("Logposterior")+
  xlab("MCMC iteration")


p1 + (plot_SBS_signature(SigsEst) / p2 + plot_layout(heights = c(2, 1)))



# I also want to plot these
as.data.frame(out_CompNMF$mcmc_out[[1]]$Mu) %>%
  mutate(iteration = row_number()) %>%
  pivot_longer(cols = -iteration) %>%
  left_join(data.frame(name = names(id_sigs), 
                       signature = colnames(SigsEst)), by = "name") %>%
  mutate(signature = case_when(!is.na(signature)~ signature, 
                          TRUE ~ name)) %>%
  ggplot(aes(x = iteration, y = value, color = signature)) +
  geom_line() 

p3 <- as.data.frame(out_CompNMF$mcmc_out[[1]]$Mu) %>%
  mutate(iteration = row_number()) %>%
  pivot_longer(cols = -iteration) %>%
  left_join(data.frame(name = names(id_sigs), 
                       signature = colnames(SigsEst)), by = "name") %>%
  mutate(signature = case_when(!is.na(signature) ~ signature, 
                               TRUE ~ name)) %>%
  
  # 1. Create the color grouping
  mutate(color_group = if_else(signature %in% c("Sig01", "Sig02", "Sig03", 
                                                "Sig04", "Sig05", "Sig06"), 
                               signature, "Other")) %>%
  
  # 2. Sort so "Other" is plotted first (FALSE/0 comes before TRUE/1)
  # This ensures the lightgray lines stay in the background
  arrange(color_group != "Other") %>%
  
  # 3. Use group = signature (so lines stay separate) and color = color_group
  ggplot(aes(x = iteration, y = value, group = signature, color = color_group)) +
  geom_line() +
  theme_bw()+
  # 4. Manually assign colors
  scale_color_manual(
    values = c(
      "Sig01" = "#F8766D", # Standard ggplot default colors
      "Sig02" = "#B79F00",
      "Sig03" = "#00BA38",
      "Sig04" = "#00BFC4",
      "Sig05" = "#619CFF",
      "Sig06" = "#F564E3",
      "Other" = "gray"
    ),
    # Optional: Removes "Other" from your legend so it looks cleaner
    breaks = c("Sig01", "Sig02", "Sig03", "Sig04", "Sig05", "Sig06") 
  )



df_res <- data.frame()
for(i in 1:length(id_sigs)){
  res <- cosine(t(out_CompNMF$mcmc_out[[1]]$Signatures[, , id_sigs[i]]), 
                as.matrix(SigsEst))
  maxRes <- apply(res, 1, which.max)
  cosine_sig <- sapply(1:nrow(res), function(j) res[j, maxRes[j]])# res[, i]
  df_res <- rbind(df_res, data.frame("iteration" = 1:nrow(res), "signature" = i, 
                                     "best_match" = maxRes, "cosine" = cosine_sig))
}

df_res$signature <- as.factor(df_res$signature)
df_res$best_match <- as.factor(df_res$best_match)

ggplot(df_res, aes(x = iteration, y = signature)) +
  geom_point(aes(color = best_match, size = cosine)) +
  theme_bw() + 
  scale_size_continuous(range = c(0.5, 5), breaks = seq(0.5, 1, by = 0.1)) +
  labs(x = "Iteration", y = "Signature")

yy <- sapply(1:nrow(res), function(j) res[j, maxRes[j]])

colnames(SigsEst) <- paste0("Sig", sprintf("%02d", 1:ncol(SigsEst)))
df_res <- data.frame()
for(i in 1:length(id_sigs)){
  res <- cosine(t(out_CompNMF$mcmc_out[[1]]$Signatures[, , id_sigs[i]]), 
                as.matrix(SigsEst))

  df_tmp <- as.data.frame(res) %>%
    mutate(iteration = row_number()) %>%
    pivot_longer(cols = -iteration)
  df_tmp$Signature <- paste0("Sig", sprintf("%02d", i))
  
  df_res <- rbind(df_res, df_tmp)
}

ggplot(df_res, aes(x = iteration, y = value, color = name))+
  geom_line() +
  facet_wrap(~ Signature) + 
  ylim(c(0, 1)) + 
  theme_bw()

plot(out_CompNMF$mcmc_out[[1]]$logposterior[1000:20000], type = "l")

#-----------------------------------------------------



iters <- 4500:5000
chain_cosines <- apply(out_CompNMF_cosmic$mcmc_out[[1]]$Signatures[iters, ,], 1, function(R) match_to_cosmic(R))
df_chain <- data.frame(do.call("rbind", chain_cosines))
df_chain$iteration <- rep(iters, each = nrow(chain_cosines[[1]]))

chain_mu <- out_CompNMF_cosmic$mcmc_out[[1]]$Mu[iters, ]
df_mu <- as.data.frame(chain_mu) %>%
  mutate(iteration = iters) %>%
  pivot_longer(cols = -iteration,
               names_to = "signature",
               values_to = "mu")

df_chain_all <- df_chain %>%
  left_join(df_mu, by = c("signature", "iteration")) %>%
  mutate(present = mu > 0.1)

ggplot(df_chain_all) +
  geom_point(aes(x = iteration, y = signature, color = best_cosmic))

df_chain_all$signature <- factor(df_chain_all$signature, levels = colnames(out_CompNMF_cosmic$mcmc_out[[1]]$init_pars$R))
df_chain_all$best_cosmic <- factor(df_chain_all$best_cosmic, levels = colnames(CompressiveNMF::COSMIC_v3.4_SBS96_GRCh37))

ggplot(df_chain_all, aes(x = iteration, y = signature)) +
  # Layer 1: Draw colored points ONLY for rows where mu >= 0.1
  geom_point(
    data = filter(df_chain_all, present == TRUE),
    aes(color = best_cosmic),
    shape = 16,
    size = 1) +
  # Layer 2: Draw a black 'x' ONLY for rows where mu < 0.1
  geom_point(
    size = 1,
    data = filter(df_chain_all, present == FALSE),
    shape = 4,        # 4 is an 'x'
    color = "black"   # Force the color to be black
  ) +
  theme_bw()

df_chain_filtered <- df_chain_all %>%
  group_by(signature) %>%
  # mean(present) calculates the fraction of iterations where present == TRUE
  filter(mean(present) > 0.75) %>%
  ungroup()

# Plot using the filtered dataset
ggplot(df_chain_filtered, aes(x = iteration, y = signature)) +
  # Layer 1: Draw colored points ONLY for rows where mu >= 0.1
  geom_point(
    data = filter(df_chain_filtered, present == TRUE),
    aes(color = best_cosmic, size = cosine_sim),
    shape = 16
  ) +
  # Layer 2: Draw a black 'x' ONLY for rows where mu < 0.1
  geom_point(
    data = filter(df_chain_filtered, present == FALSE),
    shape = 4,
    color = "black"
  ) +
  theme_bw() +
  labs(
    title = "Signatures present in >25% of iterations",
    x = "Iteration",
    y = "Signature"
  )

# One can use MatchAlign if label signature is detected. 





