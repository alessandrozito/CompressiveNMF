#------------------------------------------------- Computation time + ESS (supplement)
#
#   (1) Computation time of CompressiveNMF as a function of the number of
#       mutation categories M and the number of samples N.
#
#   (2) Effective sample size (ESS) of the CompressiveNMF MCMC algorithm. We use
#       the ESS of the relevance weight (mu) chain of the selected signatures.
#
#-------------------------------------------------

library(CompressiveNMF)
library(tidyverse)
library(ggpubr)
library(grid)
library(foreach)
library(doParallel)
library(coda)

source("R/Postprocess_functions.R")

create_directory <- function(dir) {
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
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

#------------------------------------------------- Data generator (from Jeff)
simulate_data2 <- function(I = 96, J = 100, K_new = 6, alpha = 0.25,
                           theta = 100, overdispersion = 0){
  # Generate the signatures
  Rmat <- t(LaplacesDemon::rdirichlet(n = K_new, alpha = rep(alpha, I)))
  colnames(Rmat) <- paste0("SBSnew", 1:K_new)
  # Generate the exposures
  K <- ncol(Rmat)
  exposures <- rgamma(K, theta, 1)
  Theta <- matrix(rgamma(K * J, 0.5, 0.5), ncol = J, nrow = K)
  Theta <- apply(Theta, 2, function(x) x * exposures)
  rownames(Theta) <- colnames(Rmat)
  # Generate the counts
  Lambda <- Rmat %*% Theta
  X <- matrix(rnbinom(length(Lambda), size = 1 / overdispersion, mu = c(Lambda)),
              nrow = I, ncol = J)
  rownames(X) <- rownames(Rmat)
  return(list(X = X, Rmat = Rmat, Theta = Theta))
}

#------------------------------------------------- ESS of the MCMC algorithm
# We report the ESS of the relevance weights (mu) for the selected signatures,
# taken from the best chain, exactly as in print.CompressiveNMF(). We also
# return the ESS of the log-posterior as an overall mixing summary.
compute_ESS_Comp <- function(res){
  id_best <- res$selected_chain
  Mu_all  <- res$mcmc_out[[id_best]]$Mu
  # Keep only the signatures retained in the final solution
  sel <- colnames(Mu_all) %in% names(res$RelWeights)
  Mu_chain <- as.matrix(Mu_all[, sel, drop = FALSE])

  if(ncol(Mu_chain) == 0){
    ess_mu <- NA
  } else {
    ess_mu <- tryCatch(coda::effectiveSize(Mu_chain),
                       error = function(e) NA)
  }
  ess_lp <- tryCatch(coda::effectiveSize(res$mcmc_out[[id_best]]$logposterior),
                     error = function(e) NA)

  c("ess_mu_min"    = min(ess_mu),
    "ess_mu_median" = median(ess_mu),
    "ess_mu_mean"   = mean(ess_mu),
    "ess_logpost"   = as.numeric(ess_lp))
}

#------------------------------------------------- Main experiment
# grid_MN: a data.frame with columns M and N. For each (M, N) row we simulate
# `nreps` replicate datasets and run the sampler once on each.
run_timing_experiment <- function(grid_MN,
                                  nreps = 10,
                                  K = 20,
                                  epsilon = 0.001,
                                  alpha = 0.5,
                                  a = 1,
                                  nsamples = 500,
                                  burnin = 3000,
                                  overdispersion = 0,
                                  num_cores = NULL,
                                  seed = 10){

  # Expand the grid so that each (M, N) appears `nreps` times
  jobs <- grid_MN[rep(1:nrow(grid_MN), each = nreps), , drop = FALSE]
  jobs$rep <- rep(1:nreps, times = nrow(grid_MN))
  rownames(jobs) <- NULL

  # Pre-simulate every dataset up front (done serially in the main process so
  # the data are fully reproducible and independent of the parallel backend)
  set.seed(seed, kind = "L'Ecuyer-CMRG")
  data_all <- lapply(1:nrow(jobs), function(i){
    simulate_data2(I = jobs$M[i], J = jobs$N[i], overdispersion = overdispersion)
  })

  # Create the cluster
  if(is.null(num_cores)){
    num_cores <- min(nrow(jobs), parallel::detectCores() - 1)
  }
  registerDoParallel(num_cores)

  # Parallelize the runs (one single-chain MCMC per replicate)
  results <- foreach(i = c(1:nrow(jobs)), .combine = "rbind") %dopar% {

    data <- data_all[[i]]
    # Run CompressiveNMF once, single chain, and record its wall-time
    res <- CompressiveNMF(data$X,
                          ncores = 1, nchains = 1,
                          K = K,
                          alpha = alpha,
                          a = a,
                          epsilon = epsilon,
                          burnin = burnin,
                          nsamples = nsamples)

    # Mixing: ESS of the relevance weights (and of the log-posterior)
    ess <- compute_ESS_Comp(res)
    # Number of signatures retained
    K_selected <- sum(res$RelWeights > 1.5 * epsilon)

    data.frame("M"             = jobs$M[i],
               "N"             = jobs$N[i],
               "rep"           = jobs$rep[i],
               "Kused"         = K,
               "Kselected"     = K_selected,
               "time"          = res$time,                 # seconds
               "ess_mu_min"    = ess[["ess_mu_min"]],
               "ess_mu_median" = ess[["ess_mu_median"]],
               "ess_mu_mean"   = ess[["ess_mu_mean"]],
               "ess_logpost"   = ess[["ess_logpost"]],
               "ess_mean_per_sec" = ess[["ess_mu_mean"]] / res$time)
  }
  rownames(results) <- NULL
  return(results)
}

#------------------------------------------------- Settings
M_values <- c(12, 24, 48, 96, 192)   # vary M, fix N
N_fixed  <- 120
N_values <- c(30, 60, 120, 240, 480)  # vary N, fix M
M_fixed  <- 96

nreps    <- 10
nsamples <- 1000
burnin   <- 2000
K        <- 20
epsilon  <- 0.001

create_directory("output/computation_time")

#------------------------------------------------- Run the simulation
rerun <- TRUE
if(rerun) {

  # (a) Computation time as a function of M, with N fixed at 120
  set.seed(10, kind = "L'Ecuyer-CMRG")
  results_M <- run_timing_experiment(grid_MN = data.frame(M = M_values, N = N_fixed),
                                     nreps = nreps, K = K, epsilon = epsilon,
                                     nsamples = nsamples, burnin = burnin,
                                     overdispersion = 0, num_cores = 20)
  results_M$experiment <- "M"
  write_csv(results_M, file = "output/computation_time/time_ess_vary_M.csv")

  # (b) Computation time as a function of N, with M fixed at 96
  set.seed(10, kind = "L'Ecuyer-CMRG")
  results_N <- run_timing_experiment(grid_MN = data.frame(M = M_fixed, N = N_values),
                                     nreps = nreps, K = K, epsilon = epsilon,
                                     nsamples = nsamples, burnin = burnin,
                                     overdispersion = 0, num_cores = 20)
  results_N$experiment <- "N"
  write_csv(results_N, file = "output/computation_time/time_ess_vary_N.csv")

}

################################################
# Figures
################################################
results_M <- read_csv("output/computation_time/time_ess_vary_M.csv")
results_N <- read_csv("output/computation_time/time_ess_vary_N.csv")

# Helper: boxplots over the replicates with the medians connected by a line
boxplot_with_median_line <- function(df, xvar, yvar, xlab, ylab){
  xv <- df[[xvar]]
  w  <- 0.6 * min(diff(sort(unique(xv))))   # box width = 60% of the smallest gap
  ggplot(df, aes(x = .data[[xvar]], y = .data[[yvar]])) +
    geom_boxplot(aes(group = .data[[xvar]]), width = w,
                 outlier.size = 0.6, fill = "gray90") +
    stat_summary(fun = median, geom = "line", aes(group = 1),
                 color = "steelblue", linewidth = 0.7) +
    stat_summary(fun = median, geom = "point", color = "steelblue", size = 1.6) +
    scale_x_continuous(breaks = sort(unique(xv))) +
    theme_bw() +
    labs(x = xlab, y = ylab)
}

#------------------------------------------------- Figure 1: computation time
results_M$timemin <- log(results_M$time/60)
results_M$logM <- log(results_N$N)
p_time_M <- boxplot_with_median_line(results_M, "logM", "timemin",
                                     xlab = "Log n. of mutation categories (M)",
                                     ylab = "Log computation time (min)") +
  ggtitle("N = 120")# +
  #ylim(c(0, 20))

results_N$timemin <- log(results_N$time/60)
results_N$logN <- log(results_N$N)
p_time_N <- boxplot_with_median_line(results_N, "logN", "timemin",
                                     xlab = "Log n. of samples (N)",
                                     ylab = "Log computation time (min)") +
  ggtitle("M = 96") #+
  #ylim(c(0, 45))

fig_time <- ggpubr::ggarrange(p_time_M, p_time_N, nrow = 1)
ggsave(plot = fig_time, filename = "figures/computation_time_M_N.pdf",
       width = 7.58, height = 3.4)

#------------------------------------------------- Figure 2: ESS of the MCMC
p_ess_M <- boxplot_with_median_line(results_M, "M", "ess_mu_median",
                                    xlab = "Number of mutation categories (M)",
                                    ylab = "ESS of relevance weights") +
  ggtitle("N = 120")

p_ess_N <- boxplot_with_median_line(results_N, "N", "ess_mu_median",
                                    xlab = "Number of samples (N)",
                                    ylab = "ESS of relevance weights") +
  ggtitle("M = 96")

fig_ess <- ggpubr::ggarrange(p_ess_M, p_ess_N, nrow = 1)
fig_ess <- annotate_figure(fig_ess,
                           top = textGrob(paste0("ESS out of ", nsamples,
                                                 " posterior samples"),
                                          gp = gpar(cex = 1)))
ggsave(plot = fig_ess, filename = "figures/ess_CompressiveNMF.pdf",
       width = 7.58, height = 3.4)
