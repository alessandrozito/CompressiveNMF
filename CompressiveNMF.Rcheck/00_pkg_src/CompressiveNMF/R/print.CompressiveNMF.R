plot_mu_chain <- function(mu_chain){
  ylim = c(min(mu_chain), max(mu_chain))
  plot(mu_chain[, 1], ylim = ylim, type = "l", ylab = "mu")
  if(ncol(mu_chain)>1){
    for(i in 2:ncol(mu_chain)){
      lines(mu_chain[, i], col = i)
    }
  }
}

plot_logpost_chain <- function(mcmc_out){
  nchains <- length(mcmc_out)
  ylims <- c(min(unlist(lapply(mcmc_out, function(x) min(x$logpost)))),
             max(unlist(lapply(mcmc_out, function(x) max(x$logpost)))))
  plot(mcmc_out[[1]]$logposterior, type = "l", ylab = "logposterior", ylim = ylims)
  if(nchains > 1){
    for(i in 2:nchains){
      lines(mcmc_out[[i]]$logposterior, col = i)
    }
  }
}

match_to_cosmic <- function(R){
  cosmic_data <- CompressiveNMF::COSMIC_SBS96_GRCh37_v3.4
  df_res <- data.frame()
  for(k in 1:ncol(R)){
    signature <- R[, k]
    dist <- apply(cosmic_data[, -c(1,2,3)], 2, function(x) lsa::cosine(signature, x))
    df_res <- rbind(df_res,
                    data.frame("signature" = colnames(R)[k],
                               "best_cosmic" = names(which.max(dist)),
                               "cosine_sim" = round(max(dist), 3)))
  }
  return(df_res)
}


# Function to print results
#' @export
print.CompressiveNMF <- function(object, ...){
  nchains <- length(object$mcmc_out)
  post <-lapply(object$mcmc_out, function(x) postprocess_mcmc_out(x, 0.05))
  logposterior <- sapply(1:nchains, function(i) post[[i]]$logpost)
  id_best <- which.max(logposterior)
  cat("Average log posteriors \n")
  print(cbind(chain = 1:nchains, logpost = logposterior, K = sapply(1:nchains, function(i) length(post[[i]]$RelWeights))))
  cat("Best solution is from chain" , which.max(logposterior), "\n")
  cat("Effective sizes of relevance weights \n")
  # Plot best mu
  par(mfrow = c(1,2))
  Mu_chain <- as.matrix(object$mcmc_out[[id_best]]$Mu[, post[[id_best]]$nonzero_sign])
  print(coda::effectiveSize(Mu_chain))
  plot_mu_chain(Mu_chain)
  # Plot logposteriors
  plot_logpost_chain(object$mcmc_out)
  par(mfrow = c(1,1))
}


