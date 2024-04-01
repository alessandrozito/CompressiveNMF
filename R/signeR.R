
# This function runs the method in signeR
library(signeR)
library(future)
run_signeR <- function(X, Kmin, Kmax, estimate_hyper = TRUE, sequential = FALSE){
  if(sequential){
    paral <- "sequential"
  } else {
    paral <- future::plan(multisession, workers = availableCores() - 1)
  }
  t_start <- Sys.time()
  out <- signeR(M = X, samples = "col", nlim = c(Kmin, Kmax),
                estimate_hyper = estimate_hyper, parallelization = paral)
  time <- difftime(Sys.time(), t_start, units = "secs")[[1]]
  out$time <- time
  return(out)
}

