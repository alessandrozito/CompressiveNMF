
library(tidyverse)
source("R/rInvKummer.R")

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
ggsave("figures/densities.pdf", height = 3.58, width=4.53)

#--------------------------------------------------------------- Panel (B)
set.seed(10)
J <- c(20, 50, 100, 500)
mean_y <- seq(0, 5, length.out = 100)
epsilon <- 0.001
df_plot3 <- data.frame()
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
  df_plot3 <- rbind(df_plot3,df_temp)
}

ggplot(df_plot3)+
  geom_line(aes(x = Y, y = V1, col = as.factor(J)))+
  geom_line(aes(x = Y, y = X10., col = as.factor(J)), linetype = "dashed")+
  geom_line(aes(x = Y, y = X90., col = as.factor(J)), linetype = "dashed")+
  geom_segment(x = 0, xend = 5, y = 0, yend = 5, linewidth = 0.1)+
  scale_color_manual("J", values = color_values)+
  theme_bw() +theme(aspect.ratio = .6, legend.position = "right")+
  facet_wrap(~"Posterior relevance weights") + 
  xlab(expression(bar(Y)[k]))+
  ylab(expression(mu[k]))
ggsave("figures/shrinkage.pdf", height = 3.58, width=4.53)


