# PCAWG MAF file
library(maftools)
library(sigminer)
library(CompressiveNMF)
library(tidyverse)
library(lsa) 
library(ggpubr)
source("R/Postprocess_functions.R")

#------------------------------------------------  Useful functions
match_signatures <- function(R, Rtarget){
  # Make sure that the order of the rows matches!
  R <- R[rownames(Rtarget), ]
  df_res <- data.frame()
  if(is.null(colnames(R))){
    colnames(R) <- paste0("Sig", 1:ncol(R))
  }
  for(k in 1:ncol(R)){
    signature <- R[, k]
    dist <- apply(Rtarget, 2, function(x) cosine(signature, x))
    df_res <- rbind(df_res, 
                    data.frame("signature" = colnames(R)[k],
                               "best_cosmic" = names(which.max(dist)),
                               "cosine_sim" = round(max(dist), 3)))
  }
  return(df_res)
}

match_signatures_Hungarian <- function(R, Rtarget){
  k1 <- ncol(R); k2 <- ncol(Rtarget)
  CosMat <- matrix(1, k1, k2)
  for (i in 1:k1) {
    for (j in 1:k2) {
      CosMat[i, j] <- 1 - cosine(R[, i], Rtarget[, j])
    }
  }
  match <- RcppHungarian::HungarianSolver(CosMat)$pairs[, 2]
  return(match[match!=0])
}

plot_sorted_loadings <- function(Wmat, order_loadings = 1:nrow(Wmat)) {
  Wmat <- Wmat[order_loadings, ]
  pW <- Wmat %>%
    as.data.frame() %>%
    mutate(patient = rownames(Wmat))%>%
    gather(key = "Signature", value = "weight", -patient) %>%
    mutate(Signature = as.factor(Signature),
           patient = as.factor(patient))%>%
    mutate(patient = ordered(patient, unique(patient))) %>%
    ggplot() +
    theme_minimal() +
    geom_bar(aes(x=patient, y = weight, fill = Signature), stat = "identity") +
    scale_x_discrete(position = "top") +
    theme(axis.text.x =  element_blank(),
          axis.title.x = element_blank(),
          legend.position = "right",
          panel.grid = element_blank())+
    ylab("Loadings")
  
  return(pW)
}

loadings_palette <- c("darkblue", "#2473D0", "#A4C1ED", "#F2DAAC", "red", "darkred", "grey45", "grey80", "grey10", "grey90")
#------------------------------------------------ 

# Load the Indel data from Pancreatic AdenoCarcinoma
df_Panc <- read_tsv(file = "data/Indels_PancAdenoCA.tsv")
Xid <- as.matrix(df_Panc[, -1])
rownames(Xid) <- df_Panc$Mutation

# Plot the total distribution of indels
CompressiveNMF:::plot.ID.signature(rowSums(Xid))

# Run all models
nsamples <- 2000
burnin <- 10000

rerun <- FALSE
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
                                a = 1)
  # Effective sample sizes
  ESS_compNMF <- get_ESS_Comp(resComp = out_CompNMF)  
  # Calculate posterior credible intervals
  lowCI_sigs <- apply(out_CompNMF$mcmc_out[[out_CompNMF$selected_chain]]$Signatures[, , colnames(out_CompNMF$Signatures)],
                 c(2,3), function(x) quantile(x, 0.05))
  highCI_sigs <- apply(out_CompNMF$mcmc_out[[out_CompNMF$selected_chain]]$Signatures[, , colnames(out_CompNMF$Signatures)],
                  c(2,3), function(x) quantile(x, 0.95))
  lowCI_theta <- apply(out_CompNMF$mcmc_out[[out_CompNMF$selected_chain]]$Weights[, colnames(out_CompNMF$Signatures), ],
                      c(2,3), function(x) quantile(x, 0.05))
  highCI_theta <- apply(out_CompNMF$mcmc_out[[out_CompNMF$selected_chain]]$Weights[, colnames(out_CompNMF$Signatures), ],
                       c(2,3), function(x) quantile(x, 0.95))
  # Remove MCMC output to save space, and add CI for theta and sigma
  for(i in 1:length(out_CompNMF$mcmc_out)){
    out_CompNMF$mcmc_out[[i]]$Signatures <- NULL
    out_CompNMF$mcmc_out[[i]]$Weights <- NULL
    out_CompNMF$mcmc_out[[i]]$Ysums <- NULL  
  }
  out_CompNMF$CI <- list(lowCI_sigs = lowCI_sigs, 
                         highCI_sigs = highCI_sigs,
                         lowCI_theta = lowCI_theta, 
                         highCI_theta = highCI_theta)
  out_CompNMF$ESS <- ESS_compNMF
  saveRDS(out_CompNMF, "output/Application_ID_panc/CompressiveNMF_nochains.rds.gzip", compress = "gzip")
  
  #--------------------- CompNMF + cosmic
  set.seed(42)
  out_CompNMF_cosmic <- CompressiveNMF(X = Xid, 
                                       K = 10, 
                                       epsilon = 0.01,
                                       use_cosmic = TRUE,
                                       nchains = 4, ncores = 4,
                                       nsamples = nsamples,
                                       burnin = burnin, 
                                       a = 1)
  # calculate effective sample sizes
  ESS_compNMFcos <- get_ESS_Comp(resComp = out_CompNMF_cosmic)  
  # Calculate posterior credible intervals
  lowCI_sigs <- apply(out_CompNMF_cosmic$mcmc_out[[out_CompNMF_cosmic$selected_chain]]$Signatures[, , colnames(out_CompNMF_cosmic$Signatures)],
                      c(2,3), function(x) quantile(x, 0.05))
  highCI_sigs <- apply(out_CompNMF_cosmic$mcmc_out[[out_CompNMF_cosmic$selected_chain]]$Signatures[, , colnames(out_CompNMF_cosmic$Signatures)],
                       c(2,3), function(x) quantile(x, 0.95))
  lowCI_theta <- apply(out_CompNMF_cosmic$mcmc_out[[out_CompNMF_cosmic$selected_chain]]$Weights[, colnames(out_CompNMF_cosmic$Signatures), ],
                       c(2,3), function(x) quantile(x, 0.05))
  highCI_theta <- apply(out_CompNMF_cosmic$mcmc_out[[out_CompNMF_cosmic$selected_chain]]$Weights[, colnames(out_CompNMF_cosmic$Signatures), ],
                        c(2,3), function(x) quantile(x, 0.95))
  # Remove MCMC output to save space, and add CI for theta and sigma
  for(i in 1:length(out_CompNMF_cosmic$mcmc_out)){
    out_CompNMF_cosmic$mcmc_out[[i]]$Signatures <- NULL
    out_CompNMF_cosmic$mcmc_out[[i]]$Weights <- NULL
    out_CompNMF_cosmic$mcmc_out[[i]]$Ysums <- NULL  
  }
  out_CompNMF_cosmic$CI <- list(lowCI_sigs = lowCI_sigs, 
                         highCI_sigs = highCI_sigs,
                         lowCI_theta = lowCI_theta, 
                         highCI_theta = highCI_theta)
  out_CompNMF_cosmic$ESS <- ESS_compNMFcos
  saveRDS(out_CompNMF_cosmic, "output/Application_ID_panc/CompressiveNMF_cosmic_nochains.rds.gzip", compress = "gzip")
  
  #--------------------- SignatureAnalyzer 
  set.seed(42)
  out_ARD <- sigminer::sig_auto_extract(t(Xid), cores = 10, K0 = 25)
  saveRDS(out_ARD, "output/Application_ID_panc/ARD.rds.gzip", compress = "gzip")
  
  set.seed(42)
  out_ARD_kl <- sigminer::sig_auto_extract(t(Xid), cores = 10, K0 = 25, method = "L1KL")
  saveRDS(out_ARD_kl, "output/Application_ID_panc/ARD_kl.rds.gzip", compress = "gzip")
  
  #--------------------- SigProfiler Extractor
  # In python3 console, run the following lines (after installing SigProfilerExtractor)
  # from SigProfilerExtractor import sigpro as sig
  # 
  # from SigProfilerExtractor import sigpro as sig
  # 
  # # Run SigProfilerExtractor on the data
  # project_name = "output/Application_ID_panc/sigprofiler_out/"  # Output directory name
  # input_type = "matrix"             # Specify the input type
  # input_data = "data/Indels_PancAdenoCA.tsv"
  # 
  # # Call the sigprofiler function
  # sig.sigProfilerExtractor(input_type, project_name, input_data,
  #                          minimum_signatures=2, maximum_signatures=15,
  #                          nmf_replicates=25)
  
}

out_CompNMF <- readRDS("output/Application_ID_panc/CompressiveNMF.rds.gzip")
out_CompNMF_cosmic <- readRDS("output/Application_ID_panc/CompressiveNMF_cosmic.rds.gzip")

# Load the outputs 
out_CompNMF <- readRDS("output/Application_ID_panc/CompressiveNMF_nochains.rds.gzip")
out_CompNMF_cosmic <- readRDS("output/Application_ID_panc/CompressiveNMF_cosmic_nochains.rds.gzip")
out_ARD <- readRDS("output/Application_ID_panc/ARD.rds.gzip")
out_ARD_kl <- readRDS("output/Application_ID_panc/ARD_kl.rds.gzip")

# Sigprofiler
sig_sigPro <- as.data.frame(read_tsv("output/Application_ID_panc/sigprofiler_out/ID83/Suggested_Solution/ID83_De-Novo_Solution/Signatures/ID83_De-Novo_Signatures.txt"))
rownames(sig_sigPro) <- sig_sigPro$MutationType
sig_sigPro <- as.matrix(sig_sigPro[, -1])
Theta_sigPro <- as.data.frame(read_tsv("output/Application_ID_panc/sigprofiler_out/ID83/Suggested_Solution/ID83_De-Novo_Solution/Activities/ID83_De-Novo_Activities_refit.txt"))
rownames(Theta_sigPro) <-Theta_sigPro$Samples
Theta_sigPro <- as.matrix(t(Theta_sigPro[, -1]))
out_SigPro <- list("Signatures" = sig_sigPro, "Loadings" = Theta_sigPro)


#------------------------------------------------------------ # Figure 5 Panel (A)
CosmicSig <- CompressiveNMF::COSMIC_v3.4_ID83_GRCh37
df_res <- rbind(data.frame(Method = "CompNMF", match_signatures(out_CompNMF$Signatures, CosmicSig), loadings = rowMeans(out_CompNMF$Weights)),
                data.frame(Method = "CompNMF+cosmic", match_signatures(out_CompNMF_cosmic$Signatures, CosmicSig), loadings = rowMeans(out_CompNMF_cosmic$Weights)),
                data.frame(Method = "SigProfiler", match_signatures(out_SigPro$Signatures, CosmicSig), loadings = rowMeans(out_SigPro$Loadings)),
                data.frame(Method = "SignatureAnalyzer", match_signatures(out_ARD$Signature.norm, CosmicSig), loadings = rowMeans(out_ARD$Exposure)))
df_res$best_cosmic <- factor(df_res$best_cosmic, levels = rev(paste0("ID", 1:23)))
df_res$Method <- factor(df_res$Method, levels = c("CompNMF", "CompNMF+cosmic", "SigProfiler", "SignatureAnalyzer"))


pcosine <-ggplot(df_res) +
  geom_point(aes(x = Method, y = best_cosmic, fill = cosine_sim, size = loadings), 
             shape = 21, stroke = 0.5) +
  theme_bw()+
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = -45, hjust = 1), 
        legend.key.height = unit(1, "null")) +
  scale_x_discrete(position = "top") +
  ylab("Closest cosmic ID signature") +
  scale_size_continuous("Avg.\nloadings", breaks = c(1, 100, 200, 300))+
  scale_fill_gradient(name = "Cosine\nsimilarity", low = "white", 
                      high = "darkblue", na.value = NA) 

pSign <- CompressiveNMF:::plot.ID.signature(out_CompNMF_cosmic$Signatures, 
                                            lowCI = out_CompNMF_cosmic$CI$lowCI_sigs,
                                            highCI = out_CompNMF_cosmic$CI$highCI_sigs) + 
  theme(strip.background = element_rect(fill = c("white")))

p_both <- ggpubr::ggarrange(pcosine, pSign, nrow = 1, widths = c(1,2))
ggsave(plot = p_both ,
       filename = "figures/Indel_pancAdeno_application.pdf", width = 9.64, height = 3.66)

#------------------------------------------------------------ Effective sample sizes
rbind(data.frame(method = "CompNMF + cosmic", t(out_CompNMF_cosmic$ESS[c(1,3,5)])),
      data.frame(method = "CompNMF", t(out_CompNMF$ESS[c(1,3,5)])))

#------------------------------------------------------------ RMSEs
rmse <- c("CompNMF" = sqrt(mean((Xid - out_CompNMF$Signatures %*% out_CompNMF$Weights)^2)),
  "CompNMFcos" = sqrt(mean((Xid - out_CompNMF_cosmic$Signatures %*% out_CompNMF_cosmic$Weights)^2)),
  "SigAnalyzer" = sqrt(mean((Xid - out_ARD$Signature.norm %*% out_ARD$Exposure)^2)),
  "SigProfiler" = sqrt(mean((Xid - out_SigPro$Signatures %*% out_SigPro$Loadings)^2)))
round(rmse, 2)

################################################## Supplements plots

#------ CompNMF + cosmic
Sig_CompNMF_cosmic <- out_CompNMF_cosmic$Signatures[rownames(CosmicSig), ]
Theta_CompNMF_cosmic <- apply(out_CompNMF_cosmic$Weights, 2, function(x) x/sum(x))
# re-arrange the signatures and re-name the
Sig_CompNMF_cosmic <- Sig_CompNMF_cosmic[, sort(colnames(Sig_CompNMF_cosmic))]
Theta_CompNMF_cosmic <- Theta_CompNMF_cosmic[sort(rownames(Theta_CompNMF_cosmic)), ]
colnames(Sig_CompNMF_cosmic)[grepl("_new", colnames(Sig_CompNMF_cosmic))] <- paste0("IDN", 1:sum(grepl("_new", colnames(Sig_CompNMF_cosmic))))
rownames(Theta_CompNMF_cosmic)[grepl("_new", rownames(Theta_CompNMF_cosmic))] <- paste0("IDN", 1:sum(grepl("_new", rownames(Theta_CompNMF_cosmic))))

# Plot the signatures
pSig_CompNMF_cos <- CompressiveNMF:::plot.ID.signature(Sig_CompNMF_cosmic)

# Keep the same clustering order to ensure visibility
Wmat <- t(apply(Theta_CompNMF_cosmic, 2, function(x) x/sum(x)))
clust <- cutree(hclust(dist(Wmat)), k = 6)
order_loadings <- order(clust)

# Plot the loadings
pW_CompNMF_cos <- plot_sorted_loadings(Wmat = Wmat, order_loadings = order_loadings) +
  scale_fill_manual(values = loadings_palette)

ggarrange(pSig_CompNMF_cos, pW_CompNMF_cos)
ggsave("figures/PancAdenoCA_CompNMFcos.pdf", width = 12.56, height = 3.44)

#------ CompNMF
Sig_CompNMF <- out_CompNMF$Signatures[rownames(CosmicSig), ]
Theta_CompNMF <- out_CompNMF$Weights
match <- match_MutSign(Sig_CompNMF_cosmic, R_hat = Sig_CompNMF)$match

Sig_CompNMF <- Sig_CompNMF[, match]; colnames(Sig_CompNMF) <- paste0("IDN", letters[1:ncol(Sig_CompNMF)])
Theta_CompNMF <- Theta_CompNMF[match, ]; rownames(Theta_CompNMF) <- paste0("IDN", letters[1:nrow(Theta_CompNMF)])

pSig_CompNMF <- CompressiveNMF:::plot.ID.signature(Sig_CompNMF)
Wmat <- t(apply(Theta_CompNMF, 2, function(x) x/sum(x)))
pW_CompNMF <- plot_sorted_loadings(Wmat = Wmat, order_loadings = order_loadings) +
  scale_fill_manual(values = loadings_palette)

ggarrange(pSig_CompNMF, pW_CompNMF)
ggsave("figures/PancAdenoCA_CompNMF.pdf", width = 12.56, height = 3.44)


#------ SigProfiler
Sig_SigPro <- out_SigPro$Signatures[rownames(CosmicSig), ]
Theta_SigPro <- out_SigPro$Loadings

# Re-order rows and columns
match <- match_signatures_Hungarian(Sig_CompNMF_cosmic, Sig_SigPro)
Sig_SigPro <- Sig_SigPro[, match]; colnames(Sig_SigPro) <- paste0("IDN", letters[1:ncol(Sig_SigPro)])
Theta_SigPro <- Theta_SigPro[match, ]; rownames(Theta_SigPro) <- paste0("IDN", letters[1:nrow(Theta_SigPro)])

# Plots
pSig_SigPro <- CompressiveNMF:::plot.ID.signature(Sig_SigPro)
Wmat <- t(apply(Theta_SigPro, 2, function(x) x/sum(x)))
pW_SigPro <- plot_sorted_loadings(Wmat = Wmat, order_loadings = order_loadings) +
  scale_fill_manual(values = loadings_palette[-1])

ggarrange(pSig_SigPro, pW_SigPro)
ggsave("figures/PancAdenoCA_SigPro.pdf", width = 12.56, height = 3.44)

#------ SignatureAnalyzer
Sig_ARD <- out_ARD$Signature.norm[rownames(CosmicSig), ]
Theta_ARD <- out_ARD$Exposure

match <- match_signatures_Hungarian(Sig_CompNMF_cosmic, Sig_ARD)
match <- c(match, setdiff(1:ncol(Sig_ARD), match))

Sig_ARD <- Sig_ARD[, match]; colnames(Sig_ARD) <- paste0("ID", letters[1:ncol(Sig_ARD)])
Theta_ARD <- Theta_ARD[match, ]; rownames(Theta_ARD) <- paste0("ID", letters[1:nrow(Theta_ARD)])

# Plots
pSig_ARD <- CompressiveNMF:::plot.ID.signature(Sig_ARD)
Wmat <- t(apply(Theta_ARD, 2, function(x) x/sum(x)))
pW_ARD <- plot_sorted_loadings(Wmat = Wmat, order_loadings = order_loadings) +
  scale_fill_manual(values = loadings_palette)

ggarrange(pSig_ARD, pW_ARD)
ggsave("figures/PancAdenoCA_ARD.pdf", width = 12.56, height = 3.44)


################################################################################
# Reprocess the ICGC data to extract the mutations. 
################################################################################
reprocessICGC_data <- FALSE
if(reprocessICGC_data){
  col_names_data <- colnames(read_tsv("../final_consensus_passonly.snv_mnv_indel.icgc.controlled.maf", n_max = 0, show_col_types = FALSE))
  
  nskip <- 1
  n_max <- 2e6
  parse_maf <- TRUE
  Mut_matrix_list <- list()
  j <- 1
  while(parse_maf) {
    print(j)
    # Read the MAF file
    pcawg_maf <- read_tsv("../final_consensus_passonly.snv_mnv_indel.icgc.controlled.maf", 
                          n_max = n_max, skip = nskip, col_names = col_names_data, 
                          show_col_types = FALSE)
    pcawg_maf$Tumor_Sample_Barcode <- paste0(pcawg_maf$Tumor_Sample_Barcode, "::", 
                                             pcawg_maf$Project_Code, "::",
                                             pcawg_maf$Donor_ID)
    maf_nrows <- nrow(pcawg_maf)
    if(maf_nrows == 0){
      break
    }
    # Make it into maf format
    maf <- maftools::read.maf(pcawg_maf)
    # Calculate the mutation matrix for all mutation types
    mutMatrix <- sig_tally(maf, mode = "ALL", cores = 20)
    Mut_matrix_list[[j]] <- mutMatrix
    
    # Update
    j <- j + 1
    nskip <- nskip + n_max
  }
  
  saveRDS(Mut_matrix_list, file = "Mut_matrix_listPCAWG.rds.gzip", compress="gzip")


  Mut_matrix_list <- readRDS("Mut_matrix_listPCAWG.rds.gzip")
  
  Mut_matrix_list[[1]]$SBS_96
  lapply(Mut_matrix_list, function(x) sum(x$SBS_96))
  
  
  Mut_matrix_list <- readRDS("Mut_matrix_listPCAWG.rds.gzip")
  Xsbs <- Mut_matrix_list[[1]]$SBS_96
  for(i in 2:length(Mut_matrix_list)){
    patients <- rownames(Xsbs)
    Xsbs_temp <- Mut_matrix_list[[i]]$SBS_96
    patients_temp <- rownames(Xsbs_temp)
    p_shared <- which(patients %in% patients_temp)
    Xsbs[patients[p_shared], ] <- Xsbs[patients[p_shared], ] + Xsbs_temp[patients[p_shared], ]
    Xsbs <- rbind(Xsbs, Xsbs_temp[patients_temp[-p_shared], ])
  }
  
  Xid <- Mut_matrix_list[[1]]$ID_83
  for(i in 2:length(Mut_matrix_list)){
    patients <- rownames(Xid)
    Xid_temp <- Mut_matrix_list[[i]]$ID_83
    patients_temp <- rownames(Xid_temp)
    p_shared <- which(patients %in% patients_temp)
    Xid[patients[p_shared], ] <- Xid[patients[p_shared], ] + Xid_temp[patients[p_shared], ]
    Xid <- rbind(Xid, Xid_temp[patients_temp[-p_shared], ])
  }
  
  
  names_all <- str_split(rownames(Xid), "::", simplify = TRUE)
  table(names_all[, 2])
  Xid_panc <- Xid[names_all[, 2] == "Panc-AdenoCA", ]
  df_id_panc <- t(Xid_panc)
  df_id_panc <- data.frame(cbind("Mutation" = rownames(df_id_panc), df_id_panc))
  write_tsv(df_id_panc, file = "data/Indels_PancAdenoCA.tsv")

}


