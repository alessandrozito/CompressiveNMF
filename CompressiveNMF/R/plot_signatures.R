## This function plots the signatures from CompressiveNMF
# Old palette
# c("#084B9D", "#2473D0", "#A4C1ED", "#F2DAAC", "#F29D52", "#E82C36")
# Cosmic palette
# c("#40BDEE", "#020202", "#E52925", "#CCC9CA", "#A3CF62", "#ECC5C5")
plot.SBS.signature <- function(signatures,
                               lowCI = NULL,
                               highCI = NULL,
                               palette = c("#40BDEE", "#020202", "#E52925", "#CCC9CA", "#A3CF62", "#ECC5C5")) {
  #load("~/CompressiveNMF/data/Cosmic_data.rdata")
  signatures <- as.matrix(signatures)
  names_sig <- rownames(signatures)
  df_plot <- data.frame(signatures) %>%
    dplyr::mutate(Channel = names_sig,
           Triplet = apply(stringr::str_split(names_sig, "", simplify = TRUE), 1,
                           function(x) paste0(x[c(1,3,7)], collapse = "")),
           Mutation = apply(stringr::str_split(names_sig, "", simplify = TRUE), 1,
                            function(x) paste0(x[c(3,4,5)], collapse = "")),
           Mutation = as.factor(Mutation)) %>%
    gather(key = "Sig", value = "Prob", -Channel, -Triplet, -Mutation)

  if(!is.null(lowCI) & !is.null(highCI)){
    df_plot <- df_plot %>%
      dplyr::left_join(data.frame(lowCI) %>%
                  dplyr::mutate(Channel = names_sig,
                                Triplet = apply(stringr::str_split(names_sig, "", simplify = TRUE), 1,
                                                function(x) paste0(x[c(1,3,7)], collapse = "")),
                                Mutation = apply(stringr::str_split(names_sig, "", simplify = TRUE), 1,
                                                 function(x) paste0(x[c(3,4,5)], collapse = "")),
                                Mutation = as.factor(Mutation)) %>%
                  gather(key = "Sig", value = "lowCI", -Channel, -Triplet, -Mutation),
              by = c("Channel", "Triplet", "Mutation", "Sig")) %>%
      dplyr::left_join(data.frame(highCI) %>%
                  dplyr::mutate(Channel = names_sig,
                                Triplet = apply(stringr::str_split(names_sig, "", simplify = TRUE), 1,
                                                function(x) paste0(x[c(1,3,7)], collapse = "")),
                                Mutation = apply(stringr::str_split(names_sig, "", simplify = TRUE), 1,
                                                 function(x) paste0(x[c(3,4,5)], collapse = "")),
                                Mutation = as.factor(Mutation)) %>%
                  gather(key = "Sig", value = "highCI", -Channel, -Triplet, -Mutation),
                by = c("Channel", "Triplet", "Mutation", "Sig"))
  }

  p <- ggplot(df_plot, aes(x = Triplet, y = Prob, fill = Mutation))+
    geom_bar(stat = "identity", width = 0.7) +
    facet_grid(Sig~Mutation, scales = "free")+
    theme_minimal()+
    scale_fill_manual(values = palette)+
    theme(
      legend.position = "none",
      axis.title = element_blank(),
      axis.text.y = element_blank(),
      axis.text.x = element_text(angle = 90, color = "gray35",
                   vjust = .5, size = 6.5, margin = margin(t = -4)),
      panel.grid = element_blank(),
      panel.spacing.x=unit(0, "lines"),
      panel.spacing.y=unit(0,"lines"))

  if(!is.null(lowCI) & !is.null(highCI)){
    p <- p +
      geom_linerange(aes(x = Triplet, ymin =lowCI, ymax = highCI), color = "grey65")
  }

  return(p)

}

plot.DBS.signature <- function(signatures,
                               lowCI = NULL,
                               highCI = NULL,
                               palette =c("#03BDEE", "#0366CD", "#A3CE62", "#076501", "#FD9798",
                                          "#E42925", "#FEB064", "#FD8004", "#CB99FC", "#4B0C9B")) {
  #load("~/CompressiveNMF/data/Cosmic_data.rdata")
  signatures <- as.matrix(signatures)
  if(is.null(colnames(signatures))){
    colnames(signatures) <- paste0("Sig", 1:ncol(signatures))
  }
  names_sig <- rownames(signatures)
  df_plot <- data.frame(signatures) %>%
    dplyr::mutate(Channel = names_sig,
                  Pair = sub(".*>", "", names_sig),
                  Mutation = as.factor(paste0(sub(">(.*)", "", names_sig), ">NN"))) %>%
    gather(key = "Sig", value = "Prob", -Channel, -Pair, -Mutation)

  if(!is.null(lowCI) & !is.null(highCI)){
    if(is.null(colnames(lowCI))){
      colnames(lowCI) <- colnames(signatures)
    }
    if(is.null(colnames(highCI))){
      colnames(highCI) <- colnames(signatures)
    }

    df_plot <- df_plot %>%
      dplyr::left_join(data.frame(lowCI) %>%
                         dplyr::mutate(Channel = names_sig,
                                       Pair = sub(".*>", "", names_sig),
                                       Mutation = as.factor(paste0(sub(">(.*)", "", names_sig), ">NN"))) %>%
                         gather(key = "Sig", value = "lowCI", -Channel, -Pair, -Mutation),
                       by = c("Channel", "Pair", "Mutation", "Sig")) %>%
      dplyr::left_join(data.frame(highCI) %>%
                         dplyr::mutate(Channel = names_sig,
                                       Pair = sub(".*>", "", names_sig),
                                       Mutation = as.factor(paste0(sub(">(.*)", "", names_sig), ">NN"))) %>%
                         gather(key = "Sig", value = "highCI", -Channel, -Pair, -Mutation),
                       by = c("Channel", "Pair", "Mutation", "Sig"))
  }

  p <- ggplot(df_plot, aes(x = Pair, y = Prob, fill = Mutation))+
    geom_bar(stat = "identity", width = 0.7) +
    facet_grid(Sig~Mutation, scales = "free", space='free_x')+
    theme_minimal()+
    scale_fill_manual(values = palette)+
    theme(
      legend.position = "none",
      axis.title = element_blank(),
      axis.text.y = element_blank(),
      axis.text.x = element_text(angle = 90, color = "gray35",
                   vjust = .5, size = 6, margin = margin(t = -5)),
      panel.grid = element_blank(),
      panel.spacing.x=unit(0.1, "lines"),
      panel.spacing.y=unit(0.2,"lines"))

  if(!is.null(lowCI) & !is.null(highCI)){
    p <- p +
      geom_linerange(aes(x = Pair, ymin =lowCI, ymax = highCI), color = "grey65")
  }

  return(p)

}

plot.ID.signature <- function(signatures,
                               lowCI = NULL,
                               highCI = NULL,
                               palette =c("#FEBE6E", "#FD8004", "#A3CE62", "#076501", "#FBCAB4",
                                          "#FC886F", "#EE4635", "#C11717", "#D1E1F2", "#94C2E0",
                                          "#4C95C8", "#1862AB", "#E4E3E8", "#B5B7D8", "#8584BD", "#614199")){
  signatures <- as.matrix(signatures)
  if(is.null(colnames(signatures))){
    colnames(signatures) <- paste0("Sig", 1:ncol(signatures))
  }
  names_sig <- rownames(signatures)

  df_plot <- data.frame(signatures) %>%
    dplyr::mutate(Channel = names_sig) %>%
    left_join(CompressiveNMF::df_Indel_names, by = 'Channel') %>%
    gather(key = "Sig", value = "Prob", -Channel, -Type, -Mutation)

  if(!is.null(lowCI) & !is.null(highCI)){
    if(is.null(colnames(lowCI))){
      colnames(lowCI) <- colnames(signatures)
    }
    if(is.null(colnames(highCI))){
      colnames(highCI) <- colnames(signatures)
    }

    df_plot <- df_plot %>%
      dplyr::left_join(data.frame(lowCI) %>%
                         dplyr::mutate(Channel = names_sig) %>%
                         left_join(CompressiveNMF::df_Indel_names, by = 'Channel') %>%
                         gather(key = "Sig", value = "lowCI", -Channel, -Type, -Mutation),
                       by = c("Channel", "Type", "Mutation", "Sig")) %>%
      dplyr::left_join(data.frame(highCI) %>%
                        dplyr::mutate(Channel = names_sig) %>%
                        left_join(CompressiveNMF::df_Indel_names, by = 'Channel') %>%
                        gather(key = "Sig", value = "highCI", -Channel, -Type, -Mutation),
                      by = c("Channel", "Type", "Mutation", "Sig"))
  }

  p <- ggplot(df_plot, aes(x = Mutation, y = Prob, fill = Type))+
    geom_bar(stat = "identity", width = 0.5) +
    facet_grid(Sig~Type, scales = "free", space='free_x')+
    theme_minimal()+
    scale_fill_manual(values = palette)+
    theme(
      legend.position = "none",
      axis.title = element_blank(),
      axis.text.y = element_blank(),
      axis.text.x = element_text(angle = 90, color = "gray35",
                                 vjust = .5, size = 6, margin = margin(t = -5)),
      panel.grid = element_blank(),
      panel.spacing.x=unit(0.1, "lines"),
      panel.spacing.y=unit(0.2,"lines"))

  if(!is.null(lowCI) & !is.null(highCI)){
    p <- p +
      geom_linerange(aes(x = Mutation, ymin =lowCI, ymax = highCI), color = "grey65")
  }
  return(p)
}


plot_weights <- function(Wmat) {
  clust <- cutree(hclust(dist(Wmat)), k = 10)
  Wmat <- Wmat[order(clust), ]
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

#' Plotting function for CompressiveNMF
#' @param object An object of class \code{CompressiveNMF} description
#' @export
plot.CompressiveNMF <- function(object, type = "signatures", ...){
  if(type == "signatures") {
    I <- nrow(as.matrix(object$Signatures))
    # Calculate the CI
    sigChain <- object$mcmc_out[[object$selected_chain]]$Signatures
    sigChain <- sigChain[, , names(object$RelWeights)]
    if(length(object$RelWeights) == 1){
      lowCI <- as.matrix(apply(sigChain, 2, function(x) quantile(x, 0.05)))
      highCI <- as.matrix(apply(sigChain, 2, function(x) quantile(x, 0.95)))
    } else {
      lowCI <- apply(sigChain, c(2,3), function(x) quantile(x, 0.05))
      highCI <- apply(sigChain, c(2,3), function(x) quantile(x, 0.95))
    }
    if(I == 96) {
      # Plot the SBS signatures mutational signatures
      plot.SBS.signature(signatures = object$Signatures, lowCI = lowCI, highCI = highCI)
    } else if (I == 78) {
      plot.DBS.signature(signatures = object$Signatures, lowCI = lowCI, highCI = highCI)
    } else if (I == 83) {
      plot.ID.signature(signatures = object$Signatures, lowCI = lowCI, highCI = highCI)
    }
  } else if (type == "loadings") {
    Wnorm <- t(apply(object$Weights, 2, function(x) x/sum(x)))
    plot_weights(Wnorm)
  }
}




