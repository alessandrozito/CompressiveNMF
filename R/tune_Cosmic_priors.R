library(tidyverse)
library(LaplacesDemon)
library(pracma)
library(foreach)
library(doParallel)
library(sigminer)
library(CompressiveNMF)

# Tuning functions
MedianCosineSim <- function(betah, sig_vec, nsamples){
  med <- median(apply(rdirichlet(nsamples, alpha = betah * sig_vec), 1, function(x) cosine(sig_vec, x)))
  return(med)
}

tune_betah <- function(sig_vec, betas = logseq(10, 5000, 1000), target_sim = 0.975, nsamples = 1000){
  for(i in 1:length(betas)){
    sim <- MedianCosineSim(betah = betas[i], sig_vec = sig_vec, nsamples = nsamples)
    diff <- sim - target_sim
    if(diff > 0){
      break
    }
  }
  return(betas[i])
}

betas <- pracma::logseq(10, 5000, 1000)
registerDoParallel(23)
nsamples <- 100

# Tune the Indel priors
SigMat <- CompressiveNMF::COSMIC_v3.4_ID83_GRCh37

set.seed(10, kind = "L'Ecuyer-CMRG")
optimal_betah <- foreach(i = 1:ncol(SigMat), .combine = "rbind") %dopar% {
  best_beta <- tune_betah(sig_vec = SigMat[, i],
                          betas = betas,
                          target_sim = 0.975,
                          nsamples = nsamples)
  data.frame("Signature" = colnames(SigMat)[i], "betah" = best_beta)
}

Betah_ID83_GRCh37_v3.4 <- optimal_betah$betah
names(Betah_ID83_GRCh37_v3.4) <- optimal_betah$Signature
save(Betah_ID83_GRCh37_v3.4, file = "data/Betah_ID83_GRCh37_v3.4.rdata")

# Tune the DBS priors
SigMat <- CompressiveNMF::COSMIC_v3.4_DBS78_GRCh37

set.seed(10, kind = "L'Ecuyer-CMRG")
optimal_betah <- foreach(i = 1:ncol(SigMat), .combine = "rbind") %dopar% {
  best_beta <- tune_betah(sig_vec = SigMat[, i],
                          betas = betas,
                          target_sim = 0.975,
                          nsamples = nsamples)
  data.frame("Signature" = colnames(SigMat)[i], "betah" = best_beta)
}

Betah_DBS78_GRCh37_v3.4 <- optimal_betah$betah
names(Betah_DBS78_GRCh37_v3.4) <- optimal_betah$Signature
save(Betah_DBS78_GRCh37_v3.4, file = "data/Betah_DBS78_GRCh37_v3.4.rdata")



