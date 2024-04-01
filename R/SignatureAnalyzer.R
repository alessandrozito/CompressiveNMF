# Run SignatureAnalyzer with the sigminer implementation


library(sigminer)
SignatureAnalyzer <- function(X, Kmax = 25, method = "L1KL"){
  t_start <- Sys.time()
  out <- sigminer::sig_auto_extract(nmf_matrix = t(X), K0 = Kmax, method = method, cores = 20)
  t_end = Sys.time() - t_start
  out$W <- out$Signature.norm
  out$H <- out$Exposure
  time <- difftime(Sys.time(), t_start, units = "secs")[[1]]
  out$time <- time
  return(out)
} 


