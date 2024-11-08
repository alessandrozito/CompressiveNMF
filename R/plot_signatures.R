################################################################################
# Plotting functions
################################################################################ 
# Useful file to plot mutational signatures
plot_cosmic_signature <- function(sign = "SBS5",
                                  cosmic_data = NULL,
                                  add_aetiology = TRUE,
                                  sign_inferred = NULL,
                                  palette =c("#084B9D", "#2473D0", "#A4C1ED", "#F2DAAC", "#F29D52", "#E82C36")){
  if(is.null(cosmic_data)){
    load("~/CompressiveNMF/data/Cosmic_data.rdata")
  }
  
  # Prepare data
  df_plot <- cosmic_data %>% 
    dplyr::select(Triplet, Mutation, all_of(sign))
  colnames(df_plot)[3] <- "prob"
  if(!is.null(sign_inferred)){
    sign_inferred <- sign_inferred/sum(sign_inferred) # Renormalize for comparison
    df_plot$prob_inferred <- sign_inferred
  }
  
  # Make the plot
  p <- ggplot(df_plot, aes(x = Triplet, y = prob, fill = Mutation))+
    geom_bar(stat = "identity") +
    facet_wrap(~Mutation, scales = "free_x", nrow = 1)+
    ylab("Mutation probability")+
    xlab("Mutational channel")+
    theme_minimal()+
    scale_fill_manual(values = palette)+
    theme(aspect.ratio = 1,
          legend.position = "none",
          axis.text.x = element_text(angle = 90, color = "gray35",
                                     vjust = .5, size = 8, margin = margin(t = -5)), 
          panel.grid.major.x = element_blank())+
    guides(fill = guide_legend(show = FALSE))
  title <- paste0("Signature ", sign)
  # Add aetiology
  if(add_aetiology){
    load("~/CompressiveNMF/data/COSMICsignatures_aetiology.rdata")
    df_aet <- aetiology_data %>% filter(signature == sign)
    caption <- paste0("Proposed aetiology: ", c(df_aet$aetiology_short))
    p <- p + labs(title = title, caption = caption)
  }  else {
    p <- p + labs(title = title)
  }
  
  if(!is.null(sign_inferred)){
    # sign_inferred <- apply(as.matrix(sign_inferred), 2, function(x) x/sum(x))
    # df_inferred <- cbind(df_plot[, 1:2], sign_inferred) %>%
    #   gather(key = "Algorithm", value = "prob", -Triplet, -Mutation) %>%
    #   group_by(Algorithm) %>%
    #   mutate(cosine_sim = round(cosine(df_plot$prob, prob), 4), 
    #          Algorithm2 = paste0(Algorithm, " - ", cosine_sim))
    p <- p + #geom_point(aes(x = Triplet, y =prob_inferred))+
      geom_errorbar(aes(x = Triplet, ymin =prob_inferred, ymax =prob_inferred))
    sim <- cosine(df_plot$prob, df_plot$prob_inferred)
    p <- p + labs(subtitle = paste0("Cosine similarity with inferred signature: ", round(sim, 4)))
  }
  
  return(p)
}

plot_signature <- function(signature,
                           sign_inferred = NULL,
                           CI = NULL,
                           palette =c("#084B9D", "#2473D0", "#A4C1ED", "#F2DAAC", "#F29D52", "#E82C36")) {
  load("~/CompressiveNMF/data/Cosmic_data.rdata")
  df_plot <- cosmic_data %>% 
    dplyr::select(Triplet, Mutation) %>%
    dplyr::mutate(prob = signature)
  if(!is.null(sign_inferred)){
    sign_inferred <- sign_inferred/sum(sign_inferred) # Renormalize for comparison
    df_plot$prob_inferred <- sign_inferred
  }
  if(!is.null(CI)){
    #sign_inferred <- sign_inferred/sum(sign_inferred) # Renormalize for comparison
    df_plot$lowCI <- CI[, 1]
    df_plot$highCI <- CI[, 2]
  }
  
  p <- ggplot(df_plot, aes(x = Triplet, y = prob, fill = Mutation))+
    geom_bar(stat = "identity") +
    facet_wrap(~Mutation, scales = "free_x", nrow = 1)+
    ylab("Mutation probability")+
    xlab("Mutational channel")+
    theme_minimal()+
    scale_fill_manual(values = palette)+
    theme(aspect.ratio = 1,
          legend.position = "none",
          axis.text.x = element_text(angle = 90, color = "gray35",
                                     vjust = .5, size = 8, margin = margin(t = -5)), 
          panel.grid.major.x = element_blank())+
    labs(title = "Signature")
  if(!is.null(sign_inferred)){
    # sign_inferred <- apply(as.matrix(sign_inferred), 2, function(x) x/sum(x))
    # df_inferred <- cbind(df_plot[, 1:2], sign_inferred) %>%
    #   gather(key = "Algorithm", value = "prob", -Triplet, -Mutation) %>%
    #   group_by(Algorithm) %>%
    #   mutate(cosine_sim = round(cosine(df_plot$prob, prob), 4), 
    #          Algorithm2 = paste0(Algorithm, " - ", cosine_sim))
    p <- p + #geom_point(aes(x = Triplet, y =prob_inferred))+
      geom_errorbar(aes(x = Triplet, ymin =prob_inferred, ymax =prob_inferred))
    sim <- cosine(df_plot$prob, df_plot$prob_inferred)
    p <- p + labs(subtitle = paste0("Cosine similarity with inferred signature: ", round(sim, 4)))
  }
  
  if(!is.null(CI)){
    p <- p + 
      geom_errorbar(aes(x = Triplet, ymin =lowCI, ymax = highCI), color = "forestgreen", width = 0.5)
  }
  
  return(p)
  
}

plot_matrix_signature <- function(sigMat = NULL,
                                  out_CompNMF = NULL,
                                  add_cosine_cosmic = TRUE,
                                  aspect.ratio = 1, 
                                  palette =c("#084B9D", "#2473D0", "#A4C1ED", "#F2DAAC", "#F29D52", "#E82C36"), 
                                  palette_errors =c("#A4C1ED", "darkblue", "#2473D0", "#F29D52", "#E82C36", "darkred")) {
  load("~/CompressiveNMF/data/Cosmic_data.rdata")
  if(!is.null(sigMat)){
    df_plot <- cosmic_data %>% 
      dplyr::select(Triplet, Mutation) %>%
      cbind(sigMat) %>%
      gather(key = "signature", value = "prob", -Triplet, -Mutation)
  } else if (!is.null(out_CompNMF)) {
    sigMat <- out_CompNMF$Signatures
    # Get the best chain
    nchains <- length(out_CompNMF$mcmc_out)
    post <-lapply(out_CompNMF$mcmc_out, function(x) postprocess_mcmc_out(x, 0.05))
    logposterior <- sapply(1:nchains, function(i) post[[i]]$logpost)
    id_best <- which.max(logposterior)
    
    # Get the posterior samples for all signatures
    chain <- out_CompNMF$mcmc_out[[id_best]]
    nonzero_sign <- which(colMeans(chain$Mu) > 0.05)
    R_chain <- chain$Signatures[, , nonzero_sign]
    
    # Calculate lowCI and highCI
    lowCI <- apply(R_chain, c(2, 3), function(x) quantile(x, 0.05))
    colnames(lowCI) <- names(nonzero_sign)
    highCI <- apply(R_chain, c(2, 3), function(x) quantile(x, 0.95))
    colnames(highCI) <- names(nonzero_sign)
    
    # Create the dataset
    df_plot <- cosmic_data %>% 
      dplyr::select(Triplet, Mutation) %>%
      cbind(sigMat) %>%
      gather(key = "signature", value = "prob", -Triplet, -Mutation) %>%
      left_join(cosmic_data %>% 
                  dplyr::select(Triplet, Mutation) %>%
                  cbind(lowCI) %>%
                  gather(key = "signature", value = "lowCI", -Triplet, -Mutation),
                by = c("Triplet", "Mutation", "signature")) %>%
      left_join(cosmic_data %>% 
                  dplyr::select(Triplet, Mutation) %>%
                  cbind(highCI) %>%
                  gather(key = "signature", value = "highCI", -Triplet, -Mutation),
                by = c("Triplet", "Mutation", "signature"))
  }
  
  if(add_cosine_cosmic){
    # Match to comsic
    data_text <- match_to_cosmic(sigMat)
    data_text$label <- paste0(data_text$best_cosmic, " - ", round(data_text$cosine_sim, 2))
    df_plot <- df_plot %>% left_join(data.frame(data_text), by = "signature")
    df_plot$label[df_plot$Mutation != "T>G"] <- NA
  }
  
  # Make the plot
  p <- ggplot(df_plot, aes(x = Triplet, y = prob, fill = Mutation))+
    geom_bar(stat = "identity", width = 0.75) +
    facet_grid(signature ~ Mutation, scales = "free")+
    ylab("Mutation probability")+
    xlab("Mutational channel")+
    theme_light()+
    scale_fill_manual(values = palette)+
    theme(
          legend.position = "none",
          axis.text.x = element_text(angle = 90, color = "gray35",
                                     vjust = .5, size = 5, margin = margin(t = -5)), 
          panel.spacing.x=unit(0, "lines"),
          panel.spacing.y=unit(0.1, "lines"), 
          panel.grid = element_blank(), 
          strip.background = element_rect(color="grey50", fill="white", linetype="solid"), 
          strip.text = element_text(colour = "grey25"))
  if(!is.null(out_CompNMF)){
    p <- p + 
      #geom_bar(aes(x = Triplet, y =lowCI), stat = "identity", width = 0.75) +
      #geom_bar(aes(x = Triplet, y =highCI), stat = "identity", width = 0.75, alpha = 0.8) +
      #geom_errorbar(aes(x = Triplet, ymin =lowCI, ymax = lowCI, color = Mutation))+
      geom_linerange(aes(x = Triplet, ymin =lowCI, ymax = highCI, color = Mutation), lty="11")+
      #geom_errorbar(aes(x = Triplet, ymin =highCI, ymax = highCI, color = Mutation))+
      scale_color_manual(values = palette_errors)
  }
  
  
  if(add_cosine_cosmic){
    # Match to comsic
    p <- p + geom_text(aes(x = Inf, y = Inf, label = label), 
                       size =1.6,
                  hjust = 1,
                  vjust = 2)
  }
  suppressWarnings(print(p))
}



plot_matrix_signature_v2 <- function(sigMat,
                                     sigChain = NULL,
                                     add_cosine_cosmic = TRUE,
                                     aspect.ratio = NULL,
                                     dims = c(2, 3),
                                     n_signs_to_add = 0,
                                     palette =c("#084B9D", "#2473D0", "#A4C1ED", "#F2DAAC", "#F29D52", "#E82C36"), 
                                     palette_errors =c("#A4C1ED", "darkblue", "#2473D0", "#F29D52", "#E82C36", "darkred")) {
  load("~/CompressiveNMF/data/Cosmic_data.rdata")
if(is.null(sigChain)){
  if(n_signs_to_add > 0){
    sigMat <- cbind(sigMat, matrix(.5, nrow = nrow(sigMat), ncol = n_signs_to_add))
  }
    df_plot <- cosmic_data %>% 
      dplyr::select(Triplet, Mutation) %>%
      cbind(sigMat) %>%
      gather(key = "signature", value = "prob", -Triplet, -Mutation)
  } else {
    # Calculate lowCI and highCI
    lowCI <- apply(sigChain, dims, function(x) quantile(x, 0.05))
    colnames(lowCI) <- colnames(sigMat)
    highCI <- apply(sigChain, dims, function(x) quantile(x, 0.95))
    colnames(highCI) <- colnames(sigMat)
    
    if(n_signs_to_add > 0){
      sigMat <- cbind(sigMat, matrix(0.5, nrow = nrow(sigMat), ncol = n_signs_to_add))
      lowCI <- cbind(lowCI, matrix(0.5, nrow = nrow(sigMat), ncol = n_signs_to_add))
      highCI <- cbind(highCI, matrix(0.5, nrow = nrow(sigMat), ncol = n_signs_to_add))
    }
    # Create the dataset
    df_plot <- cosmic_data %>% 
      dplyr::select(Triplet, Mutation) %>%
      cbind(sigMat) %>%
      gather(key = "signature", value = "prob", -Triplet, -Mutation) %>%
      left_join(cosmic_data %>% 
                  dplyr::select(Triplet, Mutation) %>%
                  cbind(lowCI) %>%
                  gather(key = "signature", value = "lowCI", -Triplet, -Mutation),
                by = c("Triplet", "Mutation", "signature")) %>%
      left_join(cosmic_data %>% 
                  dplyr::select(Triplet, Mutation) %>%
                  cbind(highCI) %>%
                  gather(key = "signature", value = "highCI", -Triplet, -Mutation),
                by = c("Triplet", "Mutation", "signature"))
  }
  
  if(add_cosine_cosmic){
    # Match to comsic
    if(n_signs_to_add>0){
      colnames(sigMat)[colnames(sigMat) == ""] <- paste0("V", 1:n_signs_to_add)  
    }
    data_text <- match_to_cosmic(sigMat)
    data_text$label <- paste0(data_text$best_cosmic, " - ", round(data_text$cosine_sim, 2))
    df_plot <- df_plot %>% left_join(data.frame(data_text), by = "signature")
    df_plot$label[df_plot$Mutation != "T>G"] <- NA
  }
  
  # Make the plot
  p <- ggplot(df_plot, aes(x = Triplet, y = prob, fill = Mutation))+
    geom_bar(stat = "identity", width = 0.75) +
    facet_grid(signature ~ Mutation, scales = "free")+
    ylab("Mutation probability")+
    xlab("Mutational channel")+
    theme_light()+
    scale_fill_manual(values = palette)+
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 90, color = "black",
                                 vjust = 0.5, size = 6, margin = margin(t = 0)), 
      axis.text.y = element_blank(), 
      panel.spacing.x=unit(0, "lines"),
      panel.spacing.y=unit(0.1, "lines"), 
      panel.grid = element_blank(), 
      strip.background = element_rect(color="grey50", fill="white", linetype="solid"), 
      strip.text = element_text(colour = "black", size = 10))
  if(!is.null(sigChain)){
    p <- p + 
      geom_linerange(aes(x = Triplet, ymin =lowCI, ymax = highCI), color = "grey60")#+
  }
  
  
  if(add_cosine_cosmic){
    # Match to comsic
    p <- p + geom_text(aes(x = Inf, y = Inf, label = label), 
                       size =3,
                       hjust = 1,
                       vjust = 2)
  }
  suppressWarnings(print(p))
}


plot_weights <- function(Wmat, col_palette = c("darkblue", "#2473D0", "#A4C1ED", "#F2DAAC" ,"#E82C36", "#F29D52",  "red4", "grey45")) {
  pW <- Wmat %>%
    as.data.frame() %>%
    mutate(patient = rownames(Wmat))%>%
    gather(key = "Signature", value = "weight", -patient) %>%
    mutate(Signature = as.factor(Signature),
           patient = as.factor(patient))%>%
    mutate(patient = ordered(patient, unique(patient))) %>%
    ggplot() +
    theme_bw() + 
    geom_bar(aes(x=patient, y = weight, fill = Signature), stat = "identity") +
    scale_x_discrete(position = "top") +
    theme(axis.text.x = element_text(angle = -40, color = "gray35",
                                     vjust = .5, size = 8.5, hjust = 0.95), axis.title.x = element_blank(), 
          legend.position = "bottom", 
          )+
    scale_fill_manual("", values = col_palette)+
    ylab("Activity %") +
    guides(fill = guide_legend(ncol = 4))
  
  return(pW)
} 


# load("~/CompressiveNMF/data/Cosmic_data.rdata")
# cosmic_data$SBS7a
# library(LaplacesDemon)
# c(rdirichlet(1, 10 * cosmic_data$SBS5))
# plot_cosmic_signature("SBS5", sign_inferred = c(rdirichlet(1, 10 * cosmic_data$SBS5)))



#plot_cosmic_signature("SBS5")
#plot_cosmic_signature("SBS7a")
#plot_cosmic_signature("SBS7b")
#plot_cosmic_signature("SBS7c")
#plot_cosmic_signature("SBS7d")
#plot_cosmic_signature("SBS25")
#plot_cosmic_signature("SBS17a")



# #-------------------------------------------------------------------------------
# out_CompNMF_cosmic2$Signatures[, 1] %*% out_CompNMF_cosmic2$Weights[1, ]
# 
# # # Function to plot the weights
# Wmat <- t(apply(out_CompNMF_cosmic2$Weights, 2, function(x) x/sum(x)))
# #Wmat <-t(out_ARD_pcawg$Exposure.norm)
# rownames(Wmat) <- colnames(X)
# 
# Wmat %>%
#   as.data.frame()%>%
#   mutate(patient = rownames(Wmat))%>%
#   gather(key = "Signature", value = "weight", -patient) %>%
#   mutate(Signature = as.factor(Signature),
#          patient = as.factor(patient))%>%
#   mutate(patient = ordered(patient, unique(patient))) %>%
#   ggplot(aes(x = patient, y = Signature, fill = weight)) +
#   geom_raster()+
#   coord_fixed(ratio = 1.5)+
#   scale_x_discrete(breaks = round(seq(1, 100, length.out = 15)))+
#   theme_bw()+
#   #scale_fill_continuous("Weight")+
#   geom_point(aes(col = weight), shape = 15, size = 2)+
#   scale_fill_gradient(low = "white", high = "white",
#                       space = "Lab",na.value = "grey50",guide="none")+
#   #scale_color_gradient2("Weight",low = "white", high = "blue",mid = "red", midpoint = 800,
#   #                  space = "Lab",na.value = "grey50",guide = "colourbar")
#   scale_color_gradient("Weight",low = "white", high = "blue",
#                         space = "Lab",na.value = "grey50",guide = "colourbar")

# sigAD <- c(read_csv("ADsig-approx.csv")$sig)
# plot_signature(sigAD) + xlab("Mutation type (single trinucleotide variant with trinucleotide context)")
# plot_signature
# 
# 
# load("~/CompressiveNMF/data/Cosmic_data_known_aetiology.rdata")
# cbind(cosmic_data$Channel, sigAD)
# 
# 
# palette <- c("#084B9D", "#2473D0", "#A4C1ED", "#F2DAAC", "#F29D52", "#E82C36")
# cosmic_data %>% 
#   dplyr::select(Triplet, Mutation) %>%
#   arrange(Mutation, Triplet)%>%
#   dplyr::mutate(prob = sigAD) %>%
#   ggplot(aes(x = Triplet, y = prob, fill = Mutation))+
#   geom_bar(stat = "identity", width = 0.7) +
#   facet_wrap(~ Mutation, scales = "free_x", nrow = 1)+
#   ylab("Mutation probability")+
#   xlab("Mutation type (single nucleotide variant with trinucleotide context)")+
#   theme_minimal()+
#   labs(title = "Mutational signature associated with Alzheimer's disease")+
#   scale_fill_manual(values = palette)+
#   theme(
#         legend.position = "none",
#         axis.text.x = element_text(angle = 90, color = "gray35",
#                                    vjust = .5, size = 5.5, margin = margin(t = -5)), 
#         panel.grid.major.x = element_blank(), 
#         plot.title = element_text(hjust = 0.5), 
#         strip.text = element_text(colour = "gray30", face = "bold"))
# ggsave("ADsig.pdf", height = 2.39, width = 7.25)
#   
# 
#   
# dfAD <- read_csv("41586_2022_4640_MOESM6_ESM.csv")
# ggplot(dfAD %>%
#           filter(`Neuron included after signature-based filtering (as shown in Ext. Data Fig. 1g)` == "Included"), 
#        aes(x = Age, 
#            y = `Estimated SNVs (per autosomal genome, post-filtering)`, 
#            color = Diagnosis, 
#            shape = Diagnosis))+
#   #geom_point(size = 2.5)+
#   geom_point(size = 1.5, stroke = .55)+
#   geom_smooth(method=lm, se=FALSE, show.legend = FALSE)+
#   theme_minimal()+
#   theme(
#         legend.position = "top")+
#   ylab("Somatic SNVs / neuron")+
#   scale_shape_manual(labels = c("AD", "Controls"), values = c(2, 1))+
#   scale_color_manual(labels = c("AD", "Controls"), values = c("#E82C36", "#084B9D"))
# ggsave("AD_vs_control.pdf", height = 2.39, width = 3.1)
# 
# 
# 
# # #-------------------------------------------------------------------------------
