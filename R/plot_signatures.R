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


