# This file samples form the inverse Kummer distribution. It uses the ratio of 
# uniform method.
library(rust)
library(bridgesampling)
library(tidyverse)

rInvKummer <- function(nsamples = 1, alpha, beta, gamma, delta) {
  # Function with the log_pdf of the kummer distribution
  log_pdf_Kummer <- function(x, alpha, beta, gamma, delta){
    (alpha - 1) * log(x) - beta * x - gamma * log(1 + delta * x)
  }
  
  # Function to sample from the Kummer distribution
  rKummer_ru <- function(nsamples = 1, alpha, beta,gamma, delta){
    x2 <- suppressWarnings(rust::ru(logf = log_pdf_Kummer, d = 1, n = nsamples, 
                   lower = 1e-12, init = 1,
                   alpha = alpha, beta = beta, gamma = gamma, delta = delta))
    return(c(x2$sim_vals))
  }
  
  nsamples <- floor(nsamples)
  samples <- 1/rKummer_ru(nsamples, alpha, beta, gamma, delta)
  return(samples)
}


dInvKummer <- function(x, alpha, beta, gamma, delta, nsamples_birdge = 1e5){
  # Pdf of the Inverse Kummer distribution
  log_pdf_InvKummer <- function(x, alpha, beta, gamma, delta){
    (- (alpha - gamma) - 1) * log(x) - gamma * log(1 + x/delta) - beta/x
  }
  # Sample from the distribution
  samples <- rInvKummer(nsamples_birdge, alpha, beta, gamma, delta)
  # Calculate the normalizing constant via bridgesampling
  mat_samples <- as.matrix(samples)
  colnames(mat_samples) <- "x"
  normconst <- bridgesampling::bridge_sampler(mat_samples,
                                              log_posterior = function(pars, data)
                                                log_pdf_InvKummer(x = pars,
                                                                  alpha = data$alpha,
                                                                  beta = data$beta,
                                                                  gamma = data$gamma,
                                                                  delta = data$delta),
                                              data = list("alpha" = alpha, "beta" = beta, "gamma" = gamma, "delta" = delta),
                                              lb = c("x" = 0),
                                              ub = c("x" = Inf), silent = T)$logml
  
  # Return the value of the density
  logpdf <- log_pdf_InvKummer(x, alpha, beta, gamma, delta) - normconst
  return(exp(logpdf))
}


############################################### Make figure for the shrinkage
make_plot <- FALSE
if(make_plot){

  # right panel, shrinkage induced by the weights
  set.seed(10)
  J <- c(20, 50, 100, 500)
  mean_y <- seq(0, 5, length.out = 200)
  df_plot2 <- data.frame()
  for(j in 1:length(J)){
    print(j)
    res <- sapply(mean_y, function(y) {
      a <- 1
      alpha <- 2 * a * J[j] + 1
      beta <- epsilon * a * J[j]
      delta <- a
      gamma <- 2 * J[j] * (y/2) + a * J[j]
      dd <- rInvKummer(10000, alpha, beta, gamma, delta)
      c(mean(dd), quantile(dd, 0.1), quantile(dd, 0.9))
    })
    df_temp <- data.frame(t(res), "J" = J[j], "Y" = mean_y)
    df_plot2 <- rbind(df_plot2,df_temp)
  }
  
  color_values <- c("darkblue", "lightblue", "#FCB92B", "#DB5902", "#C30027") 
  df_plot$Mutations <- as.factor(df_plot$Mutations)
  df_plot$J <- as.factor(df_plot$J)
  ggplot(df_plot) +
    #geom_segment(x = epsilon, xend = epsilon, y = 0, yend = Inf, linetype = "dashed", color = "gray")+
    geom_line(aes(x = x, y = density, color = Mutations, linetype = J), linewidth = 0.75)+
    theme_bw()+
    theme(aspect.ratio = .6, legend.position = "top")+
    ylab("Density")+
    facet_grid(~"Density of the Inverse Kummer")+
    scale_color_manual(expression(bar(Y)[k]), values = color_values)+
    xlab(expression(mu[k]))
  ggsave("figures/InvKummer.pdf", height = 3.58, width=4.53)
  
  ggplot(df_plot2)+
    geom_line(aes(x = Y, y = V1, col = as.factor(J)))+
    geom_line(aes(x = Y, y = X10., col = as.factor(J)), linetype = "dashed")+
    geom_line(aes(x = Y, y = X90., col = as.factor(J)), linetype = "dashed")+
    geom_segment(x = 0, xend = 5, y = 0, yend = 5, linewidth = 0.1)+
    scale_color_manual("J", values = color_values)+
    theme_bw() +theme(aspect.ratio = .6, legend.position = "top")+
    facet_wrap(~"Posterior relevance weights") + 
    xlab(expression(bar(Y)[k]))+
    ylab(expression(mu[k]))
  ggsave("figures/shrinkage.pdf", height = 3.58, width=4.53)
  

}


















