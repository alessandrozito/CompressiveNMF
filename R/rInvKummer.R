# This file samples form the inverse Kummer distribution. It uses the ratio of 
# uniform method.
library(rust)

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

