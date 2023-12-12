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
    load("data/Cosmic_data.rdata")
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
    load("data/COSMICsignatures_aetiology.rdata")
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
                           palette =c("#084B9D", "#2473D0", "#A4C1ED", "#F2DAAC", "#F29D52", "#E82C36")) {
  load("data/Cosmic_data.rdata")
  df_plot <- cosmic_data %>% 
    dplyr::select(Triplet, Mutation) %>%
    dplyr::mutate(prob = signature)
  if(!is.null(sign_inferred)){
    sign_inferred <- sign_inferred/sum(sign_inferred) # Renormalize for comparison
    df_plot$prob_inferred <- sign_inferred
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
  return(p)
  
}


#plot_cosmic_signature("SBS5")
#plot_cosmic_signature("SBS7a")
#plot_cosmic_signature("SBS7b")
#plot_cosmic_signature("SBS7c")
#plot_cosmic_signature("SBS7d")
#plot_cosmic_signature("SBS25")
#plot_cosmic_signature("SBS17a")
