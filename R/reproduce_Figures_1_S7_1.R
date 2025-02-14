# In this file, we plot the compressive "elbow" for varing values of a. We also
# include the non-compressive case. 
library(foreach)
library(doParallel)
library(tidyverse)

# Source function to sample from the Inverse Kummer
source("R/rInvKummer.R")


# Figure 1 Panel (A)
#--------------------------------------------------------------- Panel (A) - bottom 

x <- seq(0.001, 8, 0.01)
y_list <- c(3, 5, 7, 10)
J_list <- c(10, 20)
df_plot <- data.frame()
a <- 1
epsilon <- 0.001
for(J in J_list){
  for(y in y_list){
    alpha <- 2 * a * J + 1
    beta <- epsilon * a * J
    delta <- a
    gamma <- 2 * J * (y/2) + a * J
    dens <- dInvKummer(x, alpha, beta, gamma, delta, nsamples_birdge = 10000)
    df_temp <- data.frame(x = x, dens = dens, y = y, J = J)
    df_plot <- rbind(df_plot, df_temp)
  }
}

#--------------------------------------------------------------- Panel (A) - top
x <- seq(0.0001, 0.002, 0.00001)
y_list <- c(0)
J_list <- c(10, 20, 100)
df_plot2 <- data.frame()
a <- 1
epsilon <- 0.001
for(J in J_list){
  for(y in y_list){
    alpha <- 2 * a * J + 1
    beta <- epsilon * a * J
    delta <- a
    gamma <- 2 * J * (y/2) + a * J
    dens <- dInvKummer(x, alpha, beta, gamma, delta, nsamples_birdge = 10000)
    df_temp <- data.frame(x = x, dens = dens, y = y, J = J)
    df_plot2 <- rbind(df_plot2, df_temp)
  }
}

#--------------------------------------------------------------- Merge

color_values <- c("darkblue", "lightblue", "#FCB92B", "#DB5902", "#C30027") 

df_plot2$values <- "Density when y = 0"
df_plot$values <- "Density when y > a"

df_all <- rbind(df_plot, df_plot2)
df_all$y = as.factor(df_all$y)
df_all$J = as.factor(df_all$J)
ggplot(df_all) +
  geom_segment(x = epsilon, xend = epsilon, y = 0, yend = Inf,linetype ="dotted", color = "gray50")+
  geom_line(aes(x = x, y = dens, linetype = J, color = y), linewidth = 0.75)+
  theme_bw()+
  theme(legend.position = "right")+
  ylab("Density")+
  facet_wrap(.~values, scales = "free", ncol = 1)+
  scale_color_manual(expression(bar(Y)[k]), values = color_values)+
  xlab(expression(mu[k]))+ylab("Posterior density")
ggsave("figures/densities.pdf", height = 2.90, width=4.10)



#----------------------------- Figure S7.1
a_range <- c(0.1, 0.5, 1, 2)
compressive <- c(TRUE, FALSE)
J <- c(5, 20, 50, 100, 200)
params <- as.data.frame(expand_grid(a_range, compressive, J))
N_fixed <- 10
nsamples <- 5000
epsilon <- 0.001
# Range of values of y for the plot
mean_y <- seq(0, 5, length.out = 100)

# Set the numer of cores to run it in parallel
ncores <- 23
registerDoParallel(ncores)

df_plot_a <- foreach(j = 1:nrow(params), .combine = "rbind")  %dopar% {
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


p_plot <- ggplot(df_plot_a %>%
                   mutate(a_verb = paste0("a = ", a),
                          comp = case_when(compressive ~ "Compressive",
                                           TRUE ~ "Fixed")))+
  geom_segment(aes(x = a, xend = a, y = 0, yend = Inf), color = "grey40", 
               linetype = "dotdash", linewidth = 0.25)+
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


#----------------------------------- Figure 1 Panel (B)
p_elbow <- ggplot(df_plot_a %>%
         filter(compressive == TRUE, 
                a >= 0.5, J %in% c(20, 50, 200)) %>%
         mutate(a_verb = paste0("a = ", a), 
                J = as.factor(J))) +
  geom_segment(aes(x = a, xend = a, y = 0, yend = Inf), color = "grey70", 
               linetype = "dotdash", linewidth = 0.25)+
  geom_line(aes(x = Y, y = V1), color = "red")+
  geom_line(aes(x = Y, y = X10., col = J, linetype = J))+
  geom_line(aes(x = Y, y = X90., col = J, linetype = J))+
  geom_segment(x = 0, xend = 5, y = 0, yend = 5, color = "black", 
               linetype = "solid", linewidth = 0.1)+
  facet_grid(~a_verb) +
  #xlim(c(0, 4)) +
  theme_test() +
  #theme(aspect.ratio = 1)+
  scale_color_manual("J", values = hcl.colors(6, palette = "Blues 3")[c(4,3,1)]) +
  #scale_linetype_manual("J", values = c("dashed", "dashed","dashed"))+
  xlab(expression(bar(Y)[k]))+
  ylab(expression(mu[k]))
ggsave(plot = p_elbow, "figures/shrinkage.pdf", height = 2.29, width=5.71)

