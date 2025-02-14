# This file samples form the inverse Kummer distribution. It uses the ratio of 
# uniform method.
library(rust)
library(bridgesampling)
library(tidyverse)

rInvKummer <- function(nsamples = 1, alpha, beta, gamma, delta) {
  # Function with the log_pdf of the kummer distribution
  log_pdf_Kummer <- function(x, alpha, beta, gamma, delta){
    (alpha - 1) * log(x) - beta * x - gamma * log(1 + delta * x)
  }
  
  # Function to sample from the Kummer distribution
  rKummer_ru <- function(nsamples = 1, alpha, beta,gamma, delta){
    x2 <- suppressWarnings(rust::ru(logf = log_pdf_Kummer, d = 1, n = nsamples, 
                   lower = 1e-12, init = 1,
                   alpha = alpha, beta = beta, gamma = gamma, delta = delta))
    return(c(x2$sim_vals))
  }
  
  nsamples <- floor(nsamples)
  samples <- 1/rKummer_ru(nsamples, alpha, beta, gamma, delta)
  return(samples)
}


dInvKummer <- function(x, alpha, beta, gamma, delta, nsamples_birdge = 1e5){
  # Pdf of the Inverse Kummer distribution
  log_pdf_InvKummer <- function(x, alpha, beta, gamma, delta){
    (- (alpha - gamma) - 1) * log(x) - gamma * log(1 + x/delta) - beta/x
  }
  # Sample from the distribution
  samples <- rInvKummer(nsamples_birdge, alpha, beta, gamma, delta)
  # Calculate the normalizing constant via bridgesampling
  mat_samples <- as.matrix(samples)
  colnames(mat_samples) <- "x"
  normconst <- bridgesampling::bridge_sampler(mat_samples,
                                              log_posterior = function(pars, data)
                                                log_pdf_InvKummer(x = pars,
                                                                  alpha = data$alpha,
                                                                  beta = data$beta,
                                                                  gamma = data$gamma,
                                                                  delta = data$delta),
                                              data = list("alpha" = alpha, "beta" = beta, "gamma" = gamma, "delta" = delta),
                                              lb = c("x" = 0),
                                              ub = c("x" = Inf), silent = T)$logml
  
  # Return the value of the density
  logpdf <- log_pdf_InvKummer(x, alpha, beta, gamma, delta) - normconst
  return(exp(logpdf))
}


