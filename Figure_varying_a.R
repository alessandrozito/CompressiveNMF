# In this file, we plot the compressive "elbow" for varing values of a. We also
# include the non-compressive case. 
library(foreach)
library(doParallel)
library(tidyverse)

# Source function to sample from the Inverse Kummer
source("R/rInvKummer.R")

color_values <- c("darkblue", "lightblue", "#FCB92B", "#DB5902", "#C30027") 
# Hyperparameters for the simulation
a_range <- c(0.1, 0.5, 1, 2)
compressive <- c(TRUE, FALSE)
J <- c(5, 20, 50, 100, 500)
params <- as.data.frame(expand_grid(a_range, compressive, J))
N_fixed <- 10
nsamples <- 5000
#nsamples <- 1000
epsilon <- 0.001
#epsilon <- 0.5
# Range of values of y for the plot
mean_y <- seq(0, 5, length.out = 100)

# Set the numer of cores to run it in parallel
ncores <- 23
registerDoParallel(ncores)

df_plot <- foreach(j = 1:nrow(params), .combine = "rbind")  %dopar% {
  # Find simulation parameters
  a <- params$a_range[j]
  compressive <- params$compressive[j] 
  J <- params$J[j] 
  # Sample from the Inverse Kummer
  if(compressive){
    a0 <- a * J + 1
    b0 <- a * J * epsilon
  } else {
    a0 <- a * N_fixed + 1
    b0 <- a * N_fixed * epsilon
  }
  # Sample from the Inverse Kummer
  res <- sapply(mean_y, function(y) {
    a <- 1
    alpha <-  a0 + a * J
    beta <- b0
    delta <- a
    gamma <- J * y + a * J
    dd <- rInvKummer(nsamples, alpha, beta, gamma, delta)
    c(mean(dd), quantile(dd, 0.1), quantile(dd, 0.9))
  })
  df_temp <- data.frame(t(res), 
                        "J" = J, 
                        "Y" = mean_y, 
                        "compressive" = compressive, 
                        "a" = a)
}


p_plot <- ggplot(df_plot %>%
         mutate(a_verb = paste0("a = ", a),
                comp = case_when(compressive ~ "Compressive",
                                 TRUE ~ "Fixed")))+
  geom_segment(aes(x = a, xend = a, y = 0, yend = Inf), color = "pink", 
               linetype = "dotted", linewidth = 0.1)+
  geom_line(aes(x = Y, y = V1, col = as.factor(J), linetype = as.factor(J)), linewidth = 0.6)+
  #geom_line(aes(x = Y, y = X10., col = as.factor(a)), linetype = "dashed")+
  #geom_line(aes(x = Y, y = X90., col = as.factor(a)), linetype = "dashed")+
  geom_segment(x = 0, xend = 5, y = 0, yend = 5, color = "red", linewidth = 0.1, linetype = "solid")+
  
  #scale_color_manual("J", values = color_values)+
  scale_color_manual("J", values = hcl.colors(6, palette = "Blues 3")) +
  scale_linetype_manual("J", values = rev(c("solid", "longdash", "dashed", "dotdash", "dotted"))) +
  theme_test() +
  theme(aspect.ratio = 1, 
        legend.position = "right", 
        panel.spacing = unit(0, "lines"))+
  facet_wrap(~"Posterior relevance weights") + 
  xlab(expression(bar(Y)[k]))+
  ylab(expression(mu[k]))+
  facet_grid(comp~a_verb)

ggsave(plot = p_plot, filename = "figures/compressive_vs_fixed_varying_a.pdf", 
       width = 7.78, height = 4.27)









