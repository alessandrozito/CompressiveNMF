# Load the libraries
library(CompressiveNMF)
library(tidyverse)
library(lsa)
library(foreach)
library(doParallel)

# Source the files
source("R/SignatureAnalyzer.R")
source("R/SigProfilerExtractor.R")
source("R/signeR.R")
source("R/SigProfilerExtractor.R")
source("R/PoissonCUSP.R")
source("R/Postprocess_functions.R")

match_to_cosmic <- function(R){
  I <- nrow(R)
  if(I == 96){
    Rcosmic <- CompressiveNMF::COSMIC_v3.4_SBS96_GRCh38
  } else if (I == 83){
    Rcosmic <- CompressiveNMF::COSMIC_v3.4_ID83_GRCh37
  } else if (I == 78) {
    Rcosmic <- CompressiveNMF::COSMIC_v3.4_DBS78_GRCh37
  }
  #load("~/CompressiveNMF/data/Cosmic_data.rdata")
  df_res <- data.frame()
  if(is.null(colnames(R))){
    colnames(R) <- paste0("Sig", 1:ncol(R))
  }
  for(k in 1:ncol(R)){
    signature <- R[, k]
    dist <- apply(Rcosmic, 2, function(x) cosine(signature, x))
    df_res <- rbind(df_res, 
                    data.frame("signature" = colnames(R)[k],
                               "best_cosmic" = names(which.max(dist)),
                               "cosine_sim" = round(max(dist), 3)))
  }
  return(df_res)
}

# Load the data
load("data/MutMatrix_IFM2009.rdata")
Xsbs <- MutMatrix_IFM2009$SBS
Xdbs <- MutMatrix_IFM2009$DBS
Xid <- MutMatrix_IFM2009$ID

########################################################
# PART 1 - Test the models
########################################################
nsamples <- 1000
burnin <- 4000
nchains <- 4
ncores <- 4

set.seed(42)
#----------------------------------------------------------- CompressiveNMF
# SBS data
out_CompNMF_SBS <- CompressiveNMF(X = Xsbs, K = 30, 
                                  nsamples = nsamples, burnin = burnin, 
                                  nchains = nchains, ncores = ncores)
# DBS data
out_CompNMF_DBS <- CompressiveNMF(X = Xdbs, K = 30, 
                                  nsamples = nsamples, burnin = burnin, 
                                  nchains = nchains, ncores = ncores, alpha = 0.01)
plot(out_CompNMF_DBS)
# Indels data
out_CompNMF_ID <- CompressiveNMF(X = Xid, K = 30, 
                                 nsamples = nsamples, burnin = burnin, 
                                 nchains = nchains, ncores = ncores, alpha = 0.5)
plot(out_CompNMF_ID)
#print(out_CompNMF_ID)

#----------------------------------------------------------- CompressiveNMF + cosmic
# SBS data
out_CompNMFcos_SBS <- CompressiveNMF(X = Xsbs, K = 15, use_cosmic = TRUE,
                                  nsamples = nsamples, burnin = burnin, 
                                  nchains = nchains, ncores = ncores)
# DBS data
out_CompNMFcos_DBS <- CompressiveNMF(X = Xdbs, K = 15, use_cosmic = TRUE,
                                  nsamples = nsamples, burnin = burnin, 
                                  nchains = nchains, ncores = ncores, alpha = 0.1)
#plot(out_CompNMFcos_DBS)
# Indels data
out_CompNMFcos_ID <- CompressiveNMF(X = Xid, K = 15, use_cosmic = TRUE,
                                 nsamples = nsamples, burnin = burnin, 
                                 nchains = nchains, ncores = ncores, alpha = 0.5)
#plot(out_CompNMFcos_ID)

#----------------------------------------------------------- ARD L1KL
out_ARDl1_SBS <- sigminer::sig_auto_extract(t(Xsbs), method = "L1KL")
out_ARDl1_DBS <- sigminer::sig_auto_extract(t(Xdbs), method = "L1KL")
out_ARDl1_ID <- sigminer::sig_auto_extract(t(Xid), method = "L1KL")

#----------------------------------------------------------- ARD L2
out_ARDl2_SBS <- sigminer::sig_auto_extract(t(Xsbs))
out_ARDl2_DBS <- sigminer::sig_auto_extract(t(Xdbs))
out_ARDl2_ID <- sigminer::sig_auto_extract(t(Xid))



match_to_cosmic(out_ARDl2_DBS$Signature.norm)
show_sig_profile(out_ARDl2_DBS$Signature.norm, mode = "DBS", style = 'cosmic')


match_to_cosmic(out_ARDl1_DBS$Signature.norm)
show_sig_profile(out_ARDl1_DBS$Signature.norm, mode = "DBS", style = 'cosmic')


match_to_cosmic(out_CompNMF_ID$Signatures)
plot(out_CompNMF_ID)

match_to_cosmic(out_CompNMFcos_ID$Signatures)
plot(out_CompNMFcos_ID)

CompressiveNMF:::plot.DBS.signature(out_ARDl2_DBS$Signature.norm)
CompressiveNMF:::plot_weights(t(out_ARDl2_DBS$Exposure.norm))


out_CompNMFcos_DBS_2 <- CompressiveNMF(X = Xdbs, K = 0, use_cosmic = TRUE, swap_prior = FALSE,
                                     nsamples = 100, burnin = 500, betah = 1e6, betah_optimal = FALSE,
                                     nchains = 1, ncores = 1, alpha = 0.5, a = 1)

plot(out_CompNMFcos_DBS_2, "loadings")
plot(out_CompNMFcos_DBS_2, "signatures")
print(out_CompNMFcos_DBS_2)
match_to_cosmic(out_CompNMFcos_DBS_2$Signatures)

plot(out_CompNMFcos_DBS_2$Signatures %*% out_CompNMFcos_DBS_2$Weights, Xdbs)

plot(out_ARDl1_DBS$Signature.norm %*% out_ARDl1_DBS$Exposure, Xdbs)


plot(out_CompNMF_SBS$Signatures %*% out_CompNMF_SBS$Weights, Xsbs)


#------- meeting with Jeff soon


# DBS data
nchains <- 1
ncores <- 1
nsamples <- 10
burnin <- 20
# Compressive case
set.seed(10)
out_CompNMF_DBS <- CompressiveNMF(X = Xdbs, 
                                  K = 20, 
                                  nsamples = nsamples, 
                                  burnin = burnin, 
                                  nchains = 4, 
                                  ncores = 4, 
                                  alpha = 0.5, 
                                  epsilon = 0.01,
                                  a = 1, cutoff_excluded = 0)
plot(out_CompNMF_DBS)
CompressiveNMF:::plot.DBS.signature(out_CompNMF_DBS$Signatures[, out_CompNMF_DBS$RelWeights > 1.5 * 0.01])
# Fixed strength
set.seed(10)
N_fixed <- 10
out_CompNMF_DBS_fixed <- CompressiveNMF(X = Xdbs, 
                                        K = 20, 
                                        a0 = N_fixed + 1,
                                        b0 = 0.01 * N_fixed,
                                        nsamples = nsamples, 
                                        burnin = burnin, 
                                        nchains = nchains, 
                                        ncores = ncores, 
                                        alpha = 0.5, 
                                        a = 1, cutoff_excluded = 0)
plot(out_CompNMF_DBS_fixed)
print(out_CompNMF_DBS_fixed)
CompressiveNMF:::plot.DBS.signature(out_CompNMF_DBS_fixed$Signatures[, out_CompNMF_DBS_fixed$RelWeights > 1.5 * 0.01])

# SignatureAnalyzer
out_ARDl1_DBS <- sigminer::sig_auto_extract(t(Xdbs), method = "L1KL")

# SigneR. # <--- does not work with DBS!!!!!!!!!
set.seed(10)
out_signeR <- signeR(Xdbs, samples = "cols", nlim = c(1, 10), estimate_hyper = FALSE)

# Sigprofiler # also apparently not working well
sigminer::sigprofiler_extract(nmf_matrix = t(Xdbs),
                              output = "temp_files/sigProDBS", 
                              range = 1:8, init_method = "random", 
                              py_path = reticulate::conda_list()[2, 2], 
                              refit_plot = FALSE, refit = FALSE, nrun = 5L,
                              sigprofiler_version = "1.1.23",is_exome = FALSE,
                              cores = 8)














