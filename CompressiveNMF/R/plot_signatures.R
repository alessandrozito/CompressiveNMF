## This function plots the signatures from CompressiveNMF
# Old palette
# c("#084B9D", "#2473D0", "#A4C1ED", "#F2DAAC", "#F29D52", "#E82C36")
# Cosmic palette
# c("#40BDEE", "#020202", "#E52925", "#CCC9CA", "#A3CF62", "#ECC5C5")
#' Plot SBS signatures
#' @param signatures Matrix of SBS signatures
#' @param lowCI Matrix of lower credible interval
#' @param highCI Matrix of lower credible interval
#' @export
plot_SBS_signature <- function(signatures,
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
    tidyr::gather(key = "Sig", value = "Prob", -Channel, -Triplet, -Mutation)

  if(!is.null(lowCI) & !is.null(highCI)){
    df_plot <- df_plot %>%
      dplyr::left_join(data.frame(lowCI) %>%
                  dplyr::mutate(Channel = names_sig,
                                Triplet = apply(stringr::str_split(names_sig, "", simplify = TRUE), 1,
                                                function(x) paste0(x[c(1,3,7)], collapse = "")),
                                Mutation = apply(stringr::str_split(names_sig, "", simplify = TRUE), 1,
                                                 function(x) paste0(x[c(3,4,5)], collapse = "")),
                                Mutation = as.factor(Mutation)) %>%
                    tidyr::gather(key = "Sig", value = "lowCI", -Channel, -Triplet, -Mutation),
              by = c("Channel", "Triplet", "Mutation", "Sig")) %>%
      dplyr::left_join(data.frame(highCI) %>%
                  dplyr::mutate(Channel = names_sig,
                                Triplet = apply(stringr::str_split(names_sig, "", simplify = TRUE), 1,
                                                function(x) paste0(x[c(1,3,7)], collapse = "")),
                                Mutation = apply(stringr::str_split(names_sig, "", simplify = TRUE), 1,
                                                 function(x) paste0(x[c(3,4,5)], collapse = "")),
                                Mutation = as.factor(Mutation)) %>%
                    tidyr::gather(key = "Sig", value = "highCI", -Channel, -Triplet, -Mutation),
                by = c("Channel", "Triplet", "Mutation", "Sig"))
  }

  p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = Triplet, y = Prob, fill = Mutation))+
    ggplot2::geom_bar(stat = "identity", width = 0.7) +
    ggplot2::facet_grid(Sig~Mutation, scales = "free")+
    ggplot2::theme_minimal()+
    ggplot2::scale_fill_manual(values = palette)+
    ggplot2::theme(
      legend.position = "none",
      axis.title = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 90, color = "gray35",
                   vjust = .5, size = 6.5, margin = ggplot2::margin(t = -4)),
      panel.grid = ggplot2::element_blank(),
      panel.spacing.x = ggplot2::unit(0, "lines"),
      panel.spacing.y = ggplot2::unit(0,"lines"))

  if(!is.null(lowCI) & !is.null(highCI)){
    p <- p +
      ggplot2::geom_linerange(ggplot2::aes(x = Triplet, ymin =lowCI, ymax = highCI), color = "grey65")
  }

  return(p)

}

#' Plot Double-based-substitution signatures
#' @param signatures Matrix of SBS signatures
#' @param lowCI Matrix of lower credible interval
#' @param highCI Matrix of lower credible interval
#' @export
plot_DBS_signature <- function(signatures,
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
    tidyr::gather(key = "Sig", value = "Prob", -Channel, -Pair, -Mutation)

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
                         tidyr::gather(key = "Sig", value = "lowCI", -Channel, -Pair, -Mutation),
                       by = c("Channel", "Pair", "Mutation", "Sig")) %>%
      dplyr::left_join(data.frame(highCI) %>%
                         dplyr::mutate(Channel = names_sig,
                                       Pair = sub(".*>", "", names_sig),
                                       Mutation = as.factor(paste0(sub(">(.*)", "", names_sig), ">NN"))) %>%
                         tidyr::gather(key = "Sig", value = "highCI", -Channel, -Pair, -Mutation),
                       by = c("Channel", "Pair", "Mutation", "Sig"))
  }

  p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = Pair, y = Prob, fill = Mutation))+
    ggplot2::geom_bar(stat = "identity", width = 0.7) +
    ggplot2::facet_grid(Sig~Mutation, scales = "free", space='free_x')+
    ggplot2::theme_minimal()+
    ggplot2::scale_fill_manual(values = palette)+
    ggplot2::theme(
      legend.position = "none",
      axis.title = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 90, color = "gray35",
                   vjust = .5, size = 6, margin = ggplot2::margin(t = -5)),
      panel.grid = ggplot2::element_blank(),
      panel.spacing.x = ggplot2::unit(0.1, "lines"),
      panel.spacing.y = ggplot2::unit(0.2,"lines"))

  if(!is.null(lowCI) & !is.null(highCI)){
    p <- p +
      ggplot2::geom_linerange(ggplot2::aes(x = Pair, ymin =lowCI, ymax = highCI), color = "grey65")
  }

  return(p)

}

#' Plot Indel signatures
#' @param signatures Matrix of SBS signatures
#' @param lowCI Matrix of lower credible interval
#' @param highCI Matrix of lower credible interval
#' @export
plot_ID_signature <- function(signatures,
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
    dplyr::left_join(CompressiveNMF::df_Indel_names, by = 'Channel') %>%
    tidyr::gather(key = "Sig", value = "Prob", -Channel, -Type, -Mutation)

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
                         dplyr::left_join(CompressiveNMF::df_Indel_names, by = 'Channel') %>%
                         tidyr::gather(key = "Sig", value = "lowCI", -Channel, -Type, -Mutation),
                       by = c("Channel", "Type", "Mutation", "Sig")) %>%
      dplyr::left_join(data.frame(highCI) %>%
                        dplyr::mutate(Channel = names_sig) %>%
                        dplyr::left_join(CompressiveNMF::df_Indel_names, by = 'Channel') %>%
                        tidyr::gather(key = "Sig", value = "highCI", -Channel, -Type, -Mutation),
                      by = c("Channel", "Type", "Mutation", "Sig"))
  }

  p <- ggplot(df_plot, ggplot2::aes(x = Mutation, y = Prob, fill = Type))+
    ggplot2::geom_bar(stat = "identity", width = 0.5) +
    ggplot2::facet_grid(Sig~Type, scales = "free", space='free_x')+
    ggplot2::theme_minimal()+
    ggplot2::scale_fill_manual(values = palette)+
    ggplot2::theme(
      legend.position = "none",
      axis.title = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 90, color = "gray35",
                                 vjust = .5, size = 6, margin = margin(t = -5)),
      panel.grid = ggplot2::element_blank(),
      panel.spacing.x= ggplot2::unit(0.1, "lines"),
      panel.spacing.y= ggplot2::unit(0.2,"lines"))

  if(!is.null(lowCI) & !is.null(highCI)){
    p <- p +
      ggplot2::geom_linerange(ggplot2::aes(x = Mutation, ymin =lowCI, ymax = highCI), color = "grey65")
  }
  return(p)
}

#' Plot the signature loadings (exposures or weights)
#' @param signatures Matrix of SBS signatures
#' @param lowCI Matrix of lower credible interval
#' @param highCI Matrix of lower credible interval
#' @export
plot_weights <- function(Wmat, clust = NULL, nclust = 10) {
  Wmat <- t(Wmat)
  if(is.null(clust)){
    clust <- stats::cutree(stats::hclust(dist(Wmat)), k = nclust)
  } else {
    if(length(clust) != nrow(Wmat)){
      stop("Length of clust must be equal to ncol(Wmat)")
    }
  }
  Wmat <- Wmat[order(clust), ]
  pW <- Wmat %>%
    as.data.frame() %>%
    dplyr::mutate(patient = rownames(Wmat))%>%
    tidyr::gather(key = "Signature", value = "weight", -patient) %>%
    dplyr::mutate(Signature = as.factor(Signature),
           patient = as.factor(patient))%>%
    dplyr::mutate(patient = ordered(patient, unique(patient))) %>%
    ggplot() +
    ggplot2::theme_minimal() +
    ggplot2::geom_bar(ggplot2::aes(x=patient, y = weight, fill = Signature), stat = "identity") +
    ggplot2::scale_x_discrete(position = "top") +
    ggplot2::theme(axis.text.x =  element_blank(),
          axis.title.x = element_blank(),
          legend.position = "right",
          panel.grid = element_blank())+
    ggplot2::ylab("Loadings")

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
      plot_SBS_signature(signatures = object$Signatures, lowCI = lowCI, highCI = highCI)
    } else if (I == 78) {
      plot_DBS_signature(signatures = object$Signatures, lowCI = lowCI, highCI = highCI)
    } else if (I == 83) {
      plot_ID_signature(signatures = object$Signatures, lowCI = lowCI, highCI = highCI)
    }
  } else if (type == "loadings") {
    Wnorm <- apply(object$Weights, 2, function(x) x/sum(x))
    plot_weights(Wnorm)
  }
}




