# Indel signatures in Multiple myeloma dataset IFM2009. 
# We consider data on indels, which is a novel analysis.

library(CompressiveNMF)
library(tidyverse)
library(lsa)

# Load the Indel data
df_MM <- read_tsv(file = "data/Indels_mutliple_myeloma.txt")
Xid <- as.matrix(df_MM[, -1])
rownames(Xid) <- df_MM$Mutation

# Plot the total distribution of indels
CompressiveNMF:::plot.ID.signature(rowSums(Xid))

# Run all models
nsamples <- 2000
burnin <- 10000

rerun <- TRUE
if(rerun){
  #--------------------- CompNMF 
  set.seed(42)
  out_CompNMF <- CompressiveNMF(X = Xid, 
                                K = 20, 
                                epsilon = 0.01,
                                use_cosmic = FALSE,
                                nchains = 4, ncores = 4,
                                nsamples = nsamples,
                                burnin = burnin, 
                                a = 0.5, alpha = 0.1)
  saveRDS(out_CompNMF, "output/Application_ID_mm/CompressiveNMF.rds.gzip", compress = "gzip")
  
  #--------------------- CompNMF + cosmic
  set.seed(42)
  out_CompNMF_cosmic <- CompressiveNMF(X = Xid, 
                                       K = 10, 
                                       epsilon = 0.01,
                                       use_cosmic = TRUE,
                                       nchains = 4, ncores = 4,
                                       nsamples = nsamples,
                                       burnin = burnin, 
                                       a = 0.5, alpha = 0.1)
  saveRDS(out_CompNMF_cosmic, "output/Application_ID_mm/CompressiveNMF_cosmic.rds.gzip", compress = "gzip")
  
  #--------------------- SignatureAnalyzer 
  set.seed(42)
  out_ARD <- sigminer::sig_auto_extract(t(Xid), cores = 10, K0 = 25)
  saveRDS(out_ARD, "output/Application_ID_mm/ARD.rds.gzip", compress = "gzip")
  
  set.seed(42)
  out_ARD_kl <- sigminer::sig_auto_extract(t(Xid), cores = 10, K0 = 25, method = "L1KL")
  saveRDS(out_ARD_kl, "output/Application_ID_mm/ARD_kl.rds.gzip", compress = "gzip")
  
  #--------------------- SigProfiler Extractor
  # In python3 console, run the following lines (after installing SigProfilerExtractor)
  # from SigProfilerExtractor import sigpro as sig
  # 
  # # Run SigProfilerExtractor on the MM data
  # project_name = "output/Application_ID_mm/sigprofiler_out/"  # Output directory name
  # input_type = "matrix"             # Specify the input type
  # input_data = "data/Indels_mutliple_myeloma.txt" 
  # 
  # # Call the sigprofiler function
  # sig.sigProfilerExtractor(input_type, project_name, input_data, 
  #                          minimum_signatures=2, maximum_signatures=15, 
  #                          nmf_replicates=10)
  
  #--------------------- SigneR returns an error.I think it is because there is 
  # a bug and it expects only SBS mutations
  # There is an error in signeR. I would just ignore it.
  # Xtemp <- rbind(Xid, matrix(0, nrow = 13, ncol = 216))
  # resSignr <- signeR(t(Xtemp), nsig = 4)
  # resSignr <- signeR(t(Xid), nsig = 4)
}

if(FALSE){
# Summarize the outputs
out_CompNMF <- readRDS("output/Application_ID_mm/CompressiveNMF.rds.gzip")
out_CompNMF_cosmic <- readRDS("output/Application_ID_mm/CompressiveNMF_cosmic.rds.gzip")
out_ARD <- readRDS("output/Application_ID_mm/ARD.rds.gzip")
out_ARD_kl <- readRDS("output/Application_ID_mm/ARD_kl.rds.gzip")


# Sigprofiler
sig_sigPro <- as.data.frame(read_tsv("output/Application_ID_mm/sigprofiler_out/ID83/Suggested_Solution/ID83_De-Novo_Solution/Signatures/ID83_De-Novo_Signatures.txt"))
rownames(sig_sigPro) <- sig_sigPro$MutationType
sig_sigPro <- as.matrix(sig_sigPro[, -1])
Theta_sigPro <- as.data.frame(read_tsv("output/Application_ID_mm/sigprofiler_out/ID83/Suggested_Solution/ID83_De-Novo_Solution/Activities/ID83_De-Novo_Activities_refit.txt"))
rownames(Theta_sigPro) <-Theta_sigPro$Samples
Theta_sigPro <- as.matrix(t(Theta_sigPro[, -1]))
out_SigPro <- list("Signatures" = sig_sigPro,
                   "Loadings" = Theta_sigPro)

plot(out_CompNMF_cosmic)
plot(out_CompNMF)
CompressiveNMF:::plot.ID.signature(out_ARD$Signature.norm)
CompressiveNMF:::plot.ID.signature(out_ARD_kl$Signature.norm)
CompressiveNMF:::plot.ID.signature(out_SigPro$Signatures)


cosmic <- CompressiveNMF::COSMIC_v3.4_ID83_GRCh37

sort(apply(CompressiveNMF::COSMIC_v3.4_ID83_GRCh37, 2, function(x) cosine(x, out_ARD$Signature.norm[, 6])))


# Plot signatures in order
sigComp <- out_CompNMF_cosmic$Signatures[, c(1, 2, 4, 3)]
colnames(sigComp) <- paste0("Sig", 1:ncol(sigComp))
CompressiveNMF:::plot.ID.signature(sigComp) 

sigARD <- out_ARD$Signature.norm[, c(2, 4, 1, 3,5,6)]
colnames(sigARD) <- paste0("Sig", 1:ncol(sigARD))
CompressiveNMF:::plot.ID.signature(sigARD)

sigSigPro <- out_SigPro$Signatures[, c(4, 3, 1, 2)]
colnames(sigSigPro) <- paste0("Sig", 1:ncol(sigSigPro))
CompressiveNMF:::plot.ID.signature(sigSigPro)


}
