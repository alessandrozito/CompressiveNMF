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
  
  set.seed(10)
  J <- c(20, 50, 100, 500)
  mean_y <- seq(0, 5, length.out = 100)
  epsilon <- 0.01
  df_plot3 <- data.frame()
  for(j in 1:length(J)){
    print(j)
    res <- sapply(mean_y, function(y) {
      a <- 0.25
      alpha <- 2 * a * J[j] + 1
      beta <- epsilon * a * J[j]
      delta <- a
      gamma <- 2 * J[j] * (y/2) + a * J[j]
      dd <- rInvKummer(1000, alpha, beta, gamma, delta)
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
  





# Plot the cdf
set.seed(10)
epsilon <- 0.01
seq_mu <- seq(0.001, 500, by = 0.001)
y_list <- c(1, 2, 3)
J_list <- c(10, 50)
df_cdf <- df_pdf <- data.frame()
for(J in J_list){
  for(y in y_list){
    a <- 1
    alpha <- 2 * a * J + 1
    beta <- epsilon * a * J
    delta <- a
    gamma <- 2 * J * (y/2) + a * J
    dd <- rInvKummer(2e5, alpha, beta, gamma, delta)
    qq2 <- ecdf(dd)
    dens <- density(dd)
    df_pdf <- rbind(df_pdf, data.frame("x" = dens$x, "pdf" =  dens$y, "y" = y, "J" = J))
    df_cdf <- rbind(df_cdf, data.frame("mu" = seq_mu, "cdf" =  qq2(seq_mu), "y" = y, "J" = J))
  }
}

library(scales)
options(scipen = 999)
color_values <- c("darkblue", "lightblue", "#FCB92B", "#DB5902", "#C30027") 
df_cdf$y <- as.factor(df_cdf$y)
df_cdf$J <- as.factor(df_cdf$J)
ggplot(df_cdf) +
  #geom_segment(x = epsilon, xend = epsilon, y = 0, yend = Inf, linetype = "dashed", color = "gray")+
  geom_line(aes(x = mu, y = cdf, color = y, linetype = J), linewidth = 0.75)+
  theme_bw()+
  theme(aspect.ratio = .6, legend.position = "right")+
  ylab("")+
  facet_grid(~"Cdf of the Inverse Kummer")+
  scale_color_manual(expression(bar(Y)[k]), values = color_values)+
  #xlab(expression(mu[k])) + 
  xlab("")+
  scale_x_continuous(trans=log10_trans())
  


# Densities for a few values y < a
df_pdf2 <- data.frame


# Densities for varying J (ie. concentration)




df_pdf$y <- as.factor(df_pdf$y)
df_pdf$J <- as.factor(df_pdf$J)
ggplot(df_pdf) +
  #geom_segment(x = epsilon, xend = epsilon, y = 0, yend = Inf, linetype = "dashed", color = "gray")+
  geom_line(aes(x = x, y = pdf, color = y, linetype = J), linewidth = 0.75)+
  theme_bw()+
  theme(aspect.ratio = .6, legend.position = "right")+
  ylab("")+
  facet_wrap(~'Density')+
  scale_color_manual(expression(bar(Y)[k]), values = color_values)+
  xlab(expression(mu[k]))


qq <- sapply(seq_mu, function(x) mean(dd<x))
plot(seq_mu, qq)

qq2
plot(qq2)


set.seed(10)
J <- 50
a <- 1
mean_y <- seq(0, a, by = 0.01)
epsilon <- 0.01
pp <- sapply(mean_y, function(y){
  alpha <- 2 * a * J + 1
  beta <- epsilon * a * J
  delta <- a
  gamma <- 2 * J * (y/2) + a * J
  dd <- rInvKummer(1000, alpha, beta, gamma, delta)
  mean(dd)
})

sol <- function(a, y, epsilon=0.01){
  b <- y + epsilon - a
  2 * a * epsilon / (sqrt(b^2 + 8 * a * epsilon) - b)
}

approx <- function(a, y, epsilon=0.01) {
  epsilon * a * (a-y)/( (a-y)^2 + epsilon * (a+y))
}
  
plot(mean_y, pp, type = "l")
lines(mean_y, sol(a, mean_y), col = "red")
lines(mean_y, approx(a, mean_y), col = "green")


df_plot2 <- data.frame()
for(j in 1:length(J)){
  print(j)
  res <- sapply(mean_y, function(y) {
    a <- 0.25
    alpha <- 2 * a * J[j] + 1
    beta <- epsilon * a * J[j]
    delta <- a
    gamma <- 2 * J[j] * (y/2) + a * J[j]
    dd <- rInvKummer(1000, alpha, beta, gamma, delta)
    c(mean(dd), quantile(dd, 0.1), quantile(dd, 0.9))
  })
  df_temp <- data.frame(t(res), "J" = J[j], "Y" = mean_y)
  df_plot2 <- rbind(df_plot2,df_temp)
}

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



set.seed(10)
y_true <- 0.01
J_seq <- seq(10, 100000, by = 100)
y_seq <- sapply(J_seq, function(J) mean(rpois(J, y_true)))
a <- 1
epsilon <- 0.01
pp <- sapply(1:length(J_seq), function(i){
  J <- J_seq[i]; y <- y_seq[i]
  alpha <- 2 * a * J + 1
  beta <- epsilon * a * J
  delta <- a
  gamma <- 2 * J * (y/2) + a * J
  dd <- rInvKummer(1000, alpha, beta, gamma, delta)
  mean(dd)
})
plot(J_seq, pp, type = "l")
abline(h = sol(a, y_true, epsilon), col = "red")



# Panel (A)
x <- seq(0.001, 8, 0.01)
y_list <- c(3, 5, 7, 10)
J_list <- c(10, 20)
df_plot <- data.frame()
for(J in J_list){
  for(y in y_list){
    a <- 1
    epsilon <- 0.01
    alpha <- 2 * a * J + 1
    beta <- epsilon * a * J
    delta <- a
    gamma <- 2 * J * (y/2) + a * J
    dens <- dInvKummer(x, alpha, beta, gamma, delta, nsamples_birdge = 10000)
    df_temp <- data.frame(x = x, dens = dens, y = y, J = J)
    df_plot <- rbind(df_plot, df_temp)
  }
}

color_values <- c("darkblue", "lightblue", "#FCB92B", "#DB5902", "#C30027") 
# df_plot$y <- as.factor(df_plot$y)
# df_plot$J <- as.factor(df_plot$J)
# p1 <- ggplot(df_plot) +
#   #geom_segment(x = epsilon, xend = epsilon, y = 0, yend = Inf, linetype = "dashed", color = "gray")+
#   geom_line(aes(x = x, y = dens, color = y, linetype = J), linewidth = 0.75)+
#   theme_bw()+
#   theme(, axis.title = element_blank(), legend.position = "right")+
#   ylab("Density")+
#   facet_grid(~"Density when y > a")+
#   scale_color_manual(expression(bar(Y)[k]), values = color_values)+
#   xlab(expression(mu[k]))



x <- seq(0.001, 0.02, 0.0001)
y_list <- c(0)
J_list <- c(10, 20, 100)
df_plot2 <- data.frame()
for(J in J_list){
  for(y in y_list){
    a <- 1
    epsilon <- 0.01
    alpha <- 2 * a * J + 1
    beta <- epsilon * a * J
    delta <- a
    gamma <- 2 * J * (y/2) + a * J
    dens <- dInvKummer(x, alpha, beta, gamma, delta, nsamples_birdge = 10000)
    df_temp <- data.frame(x = x, dens = dens, y = y, J = J)
    df_plot2 <- rbind(df_plot2, df_temp)
  }
}

color_values <- c("darkblue", "lightblue", "#FCB92B", "#DB5902", "#C30027") 
# df_plot2$y <- as.factor(df_plot2$y)
# df_plot2$J <- as.factor(df_plot2$J)
# p2 <- ggplot(df_plot2) +
#   geom_segment(x = epsilon, xend = epsilon, y = 0, yend = Inf, color = "gray")+
#   geom_line(aes(x = x, y = dens, linetype = J), linewidth = 0.75)+
#   theme_bw()+
#   theme(legend.position = "right", axis.title = element_blank())+
#   ylab("Density")+
#   facet_grid(~"Density when y = 0")+
#   scale_color_manual(expression(bar(Y)[k]), values = color_values)
# p2

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

# Panel (B)
set.seed(10)
J <- c(20, 50, 100, 500)
mean_y <- seq(0, 5, length.out = 100)
epsilon <- 0.01
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








}
