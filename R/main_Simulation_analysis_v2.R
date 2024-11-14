# This file replicates Figure 2 in the main document and Figures S1-S4 in the supplement
library(tidyverse)

# Source the functions
source("R/SignatureAnalyzer.R")
source("R/SigProfilerExtractor.R")
source("R/signeR.R")
source("R/PoissonCUSP.R")
source("R/CompressiveNMF.R")
source("R/plot_signatures.R")
source("R/Postprocess_functions.R")
source("R/plot_signatures.R")

color_values <- c("darkblue", "royalblue", "#02DB9E", "#FCB92B", "#DB5902", "#C30027", "darkred")
labels <- c("CompNMF + cosmic", "CompNMF", "signeR", "SigProfiler", "SignatureAnalyzer", "PoissonCUSP", "BNMF")

# Load the output
df_F1 <- read_csv(file = "output/main_simulation/df_F1.csv")
df_all <- read_csv(file = "output/main_simulation/simulation_output.csv")

#-------- Figure 2 - Panel (A) - Plot for number of signatures K
df_all %>%
  filter(theta == 100) %>%
  mutate(J = as.factor(J),
         overd = paste0("Overd.=", as.factor(overd)), 
         K_tot = 4 + K_new,
         J = factor(paste0("J = ", J), levels=c("J = 50", "J = 100", "J = 200")),
         K_new2 = as.factor(paste0("K = ", K_new + 4)),
         K_new2 = fct_reorder(K_new2, as.integer(K_new))) %>%
  dplyr::select(Method, J, overd, K_new2, K,K_tot) %>%
  ggplot() +
  geom_hline(aes(yintercept = K_tot), color = "grey55")+
  geom_boxplot(aes(y = K, color = Method, fill = Method), alpha = 0.4, outlier.shape = NA)+
  theme_bw() +
  scale_color_manual(values = color_values, labels = labels)+
  scale_fill_manual(values = color_values, labels = labels)+
  facet_grid(K_new2~overd+J, scales = "free")+
  scale_y_continuous(breaks = seq(from = 0, to = 20, by = 2))+
  theme(legend.position ="top",
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-5,-5,-5,-5),
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(), 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        panel.spacing = unit(0, "lines"))+
  guides(color = guide_legend(title.position = "left", title.hjust = 1, nrow = 1))+
  ylab("Estimated number of signatures")
#ggsave(p_Ktot, filename = "figures/K_tot.pdf", width = 8.55, height = 3.15)


#-------- Figure 2 - Panel (B) - Plot for sensitivity and precision
x <- seq(0.001, 0.99999, length.out = 100)
y <- seq(0.001, 0.99999, length.out = 100)
grid <- expand.grid(x = x, y = y)
grid$f <- 2 * grid$x * grid$y / (grid$x + grid$y)

p_sens_prec <- df_all %>%
  dplyr::select(Method, J, overd, K_new, Precision, Sensitivity) %>%
  mutate(J = as.factor(J),
         overd = paste0("Overdispersion = ", as.factor(overd)), 
         K_new2 = paste0("K = ", K_new + 4),
         K_new2 = fct_reorder(K_new2, as.integer(K_new)),
         Method = case_when(Method == "4.SigPro" ~ "3.SigPro",
                            Method == "3.signeR" ~ "4.signeR",
                            Method == "6.CUSP" ~ "1.CUSP",
                            Method == "5.ARD" ~ "2.ARD",
                            Method == "1.CompNMFcos" ~ "6.CompNMFcos",
                            Method == "2.CompNMF" ~ "5.CompNMF", 
                            Method == "7.BayesNMF" ~ "7.BayesNMF")) %>%
  group_by(Method, K_new2, overd, J) %>% 
  summarise_all(mean) %>%
  ggplot() +
  scale_color_manual(values = rev(color_values), labels = rev(labels))+
  geom_contour(data = grid, aes(x = x, y = y, z = f), linetype = "dashed", alpha = 0.3, color = "red") +
  metR::geom_text_contour(data = grid, aes(x = x, y = y, z = f), stroke = 0.2, alpha = 1, color = "pink")+
  geom_point(aes(x = Precision, y = Sensitivity, color = Method, shape = J), 
             size = 2, stroke = .85) +
  theme_bw() +
  facet_grid(~overd+K_new2)+
  scale_shape_manual(values = c(1, 0, 2))+
  xlim(c(0.4, 1))+ ylim(c(0.75, 1))+
  xlab("Precision")+ylab("Sensitivity")+
  theme(legend.position ="top",
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-5,-5,-5,-5))+
  guides(color = guide_legend(reverse=TRUE))
ggsave(p_sens_prec, filename = "figures/Prec_sensitivity_plot.pdf", width = 8.51, height = 3.15)


df_all %>%
  dplyr::select(Method, J, overd, K_new, Precision, Sensitivity) %>%
  mutate(J = as.factor(J),
         overd = paste0("Overdispersion = ", as.factor(overd)), 
         K_new2 = paste0("K = ", K_new + 4),
         K_new2 = fct_reorder(K_new2, as.integer(K_new)),
         Method = case_when(Method == "4.SigPro" ~ "4.SigPro",
                            Method == "3.signeR" ~ "5.signeR",
                            Method == "6.CUSP" ~ "2.CUSP",
                            Method == "5.ARD" ~ "3.ARD",
                            Method == "1.CompNMFcos" ~ "7.CompNMFcos",
                            Method == "2.CompNMF" ~ "6.CompNMF", 
                            Method == "7.BayesNMF" ~ "1.BayesNMF")) %>%
  group_by(Method, K_new2, overd, J) %>% 
  summarise_all(mean) %>%
  ggplot() +
  scale_color_manual(values = rev(color_values), labels = rev(labels))+
  geom_contour(data = grid, aes(x = x, y = y, z = f), linetype = "dashed", alpha = 0.5, color = "black") +
  metR::geom_text_contour(data = grid, aes(x = x, y = y, z = f), stroke = 0.2, alpha = 1, color = "grey55")+
  geom_point(aes(x = Precision, y = Sensitivity, color = Method, shape = J), 
             size = 2, stroke = 1) +
  theme_bw() +
  facet_grid(~overd+K_new2)+
  scale_shape_manual(values = c(1, 0, 2))+
  xlim(c(0.2, 1))+ ylim(c(0.2, 1))+
  xlab("Precision")+ylab("Sensitivity")+
  theme(legend.position ="top",
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-5,-5,-5,-5))+
  guides(color = guide_legend(reverse=TRUE))

#-------- Figure 2 - Panel (C) - Plot for F1 score when overdispersion = 0.15
pF1 <- df_F1 %>%
  mutate(J = as.factor(J),
         overd = paste0("Overdispersion=", as.factor(overd)), 
         J = factor(paste0("J = ", J), levels=c("J = 50", "J = 100", "J = 200")),
         K_new2 = paste0("K = ", K_new + 4),
         K_new2 = fct_reorder(K_new2, as.integer(K_new)),) %>%
  ggplot()+
  geom_line(aes(x = cutoff, y = mean, color = Method), linewidth = .6)+
  facet_grid(overd ~ K_new2 + J)+
  theme_bw()+
  scale_color_manual(values = color_values, labels = labels) +
  xlab("Cutoff")+
  ylab("Average F1 score")+
  theme(legend.position ="top",
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-5,-5,-5,-5), 
        axis.text = element_text(size = 5.5))+
  xlim(c(0.8, 1))+
  guides(color = guide_legend(title.position = "left", title.hjust = 0.5, nrow = 1))
ggsave(plot = pF1, "figures/F1_score.pdf", width = 8.51, height = 2.50)



#-------- Figure 3 - Execution time
ptime <- df_all %>%
  group_by(Method, J) %>%
  dplyr::select(time) %>%
  summarise(time = mean(time/60)) %>% 
  as.data.frame() %>%
  ggplot()+
  theme_bw()+
  facet_grid(.~"Computational time")+
  geom_point(aes(x = J, y = time, color = Method))+
  geom_line(aes(x = J, y = time, color = Method))+
  ylab("Average time (m)")+
  xlab("Sample size J")+
  scale_color_manual(values = color_values, labels = labels)
ggsave(plot = ptime, "figures/Simulation_time.pdf", width = 5.25, height = 3.17)


# Effective sample size
pESS <- df_all %>%
  dplyr::select(Method, J, overd, ESS_Sig_mean, ESS_Theta_mean, ESS_relweight_mean) %>%
  gather(key = "quantity", value = "ESS", -Method, -J, -overd) %>%
  drop_na() %>%
  mutate(quantity = case_when(quantity == "ESS_Sig_mean" ~ "Signatures", 
                              quantity == "ESS_Theta_mean" ~ "Loadings", 
                              quantity == "ESS_relweight_mean" ~ "Relevance weights"), 
         quantity = factor(quantity, levels = c("Signatures", "Loadings","Relevance weights")), 
         Method = case_when(Method == "4.SigPro" ~ "4.SigPro",
                            Method == "3.signeR" ~ "signeR",
                            Method == "6.CUSP" ~ "2.CUSP",
                            Method == "5.ARD" ~ "3.ARD",
                            Method == "1.CompNMFcos" ~ "CompNMF + cosmic",
                            Method == "2.CompNMF" ~ "CompNMF", 
                            Method == "7.BayesNMF" ~ "BayesNMF"), 
         overd = paste0("Overd.=", as.factor(overd)),
         J = as.factor(J)) %>%
  ggplot() +
  geom_boxplot(aes(x = Method, y = ESS, color = J)) +
  facet_grid(overd~quantity) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  ylab("Average effective sample size")+
  xlab("")
ggsave(plot = pESS, "figures/Simulation_effectiveSize.pdf", width = 8.88, height = 4.40)

################################################################################
# Figures in the supplementary material
################################################################################

p_overd0 <- df_all %>%
  filter(theta == 100, overd == 0) %>%
  dplyr::select(Method, J, K_new, overd, rmse_Signatures, rmse_Weights) %>%
  mutate(J = factor(paste0("J = ", J), levels=c("J = 50", "J = 100", "J = 200")),
         Signatures = rmse_Signatures,
         Loadings = rmse_Weights,
         overd = paste0("Overdispersion = ", as.factor(overd)), 
         K_new2 = paste0("K = ", K_new + 4),
         K_new2 = fct_reorder(K_new2, as.integer(K_new))) %>%
  dplyr::select(-rmse_Signatures, -rmse_Weights, -K_new) %>%
  gather(key = "Metric", value = "value", -Method, -J, -overd, -K_new2) %>%
  ggplot(aes(x = J, y = value, fill = Method, color = Method))+
  geom_boxplot(alpha = 0.4)+
  theme_bw()+
  scale_color_manual( values = color_values, labels = labels)+
  scale_fill_manual(values = color_values, labels = labels)+
  facet_grid(Metric ~ K_new2, scales = "free_y") +
  xlab("Number of samples J")+
  ylab("RMSE")+
  theme(
    legend.position ="right",
    legend.margin=margin(0,0,0,0),
    legend.box.margin=margin(-5,-5,-5,-5),
    panel.grid.major.x = element_blank(), 
    panel.grid.minor.x = element_blank(), 
  )+
  labs(title = "Overdispersion = 0")
ggsave(plot = p_overd0, filename = "figures/RMSE_plot_overd0.pdf", width = 8.51, height = 3.61)

p_overd15 <- df_all %>%
  filter(theta == 100, overd == 0.15) %>%
  dplyr::select(Method, J, K_new, overd, rmse_Signatures, rmse_Weights) %>%
  mutate(J = factor(paste0("J = ", J), levels=c("J = 50", "J = 100", "J = 200")),
         Signatures = rmse_Signatures,
         Loadings = rmse_Weights,
         overd = paste0("tau = ", as.factor(overd)), 
         K_new2 = paste0("K = ", K_new + 4),
         K_new2 = fct_reorder(K_new2, as.integer(K_new))) %>%
  dplyr::select(-rmse_Signatures, -rmse_Weights, -K_new) %>%
  gather(key = "Metric", value = "value", -Method, -J, -overd, -K_new2) %>%
  ggplot(aes(x = J, y = value, fill = Method, color = Method))+
  geom_boxplot(alpha = 0.4)+
  theme_bw()+
  scale_color_manual( values = color_values, labels = labels)+
  scale_fill_manual(values = color_values, labels = labels)+
  facet_grid(Metric ~ K_new2 , scales = "free_y") +
  #theme(axis.title.y = element_blank())+
  xlab("Number of samples J")+
  ylab("RMSE")+
  theme(#aspect.ratio = 1,
    legend.position ="right",
    legend.margin=margin(0,0,0,0),
    legend.box.margin=margin(-5,-5,-5,-5),
    panel.grid.major.x = element_blank(), 
    panel.grid.minor.x = element_blank(), 
  )+
  labs(title = "Overdispersion = 0.15")
ggsave(plot = p_overd15, filename = "figures/RMSE_plot_overd15.pdf", width = 8.51, height = 3.61)

#---- Counts and Lambda
p_counts0 <- df_all %>%
  filter(theta == 100, overd == 0) %>%
  dplyr::select(Method, J, K_new, overd, rmse_Lambda, rmse_Counts) %>%
  mutate(J = as.factor(J),
         Lambda = rmse_Lambda,
         Counts = rmse_Counts,
         overd = paste0("Overdispersion = ", as.factor(overd)), 
         K_new2 = paste0("K = ", K_new + 4),
         K_new2 = fct_reorder(K_new2, as.integer(K_new))) %>%
  dplyr::select(-rmse_Counts, -rmse_Lambda, -K_new) %>%
  gather(key = "Metric", value = "value", -Method, -J, -overd, -K_new2) %>%
  ggplot(aes(x = J, y = value, fill = Method, color = Method))+
  geom_boxplot(alpha = 0.4)+
  theme_bw()+
  scale_color_manual( values = color_values, labels = labels)+
  scale_fill_manual(values = color_values, labels = labels)+
  facet_grid(Metric ~ K_new2, scales = "free_y") +
  theme(#aspect.ratio = 1,
    legend.position ="right",
    legend.margin=margin(0,0,0,0),
    legend.box.margin=margin(-5,-5,-5,-5),
    panel.grid.major.x = element_blank(), 
    panel.grid.minor.x = element_blank(), 
  )+
  xlab("Number of samples J")+
  ylab("RMSE")+
  labs(title = "Overdispersion = 0")
ggsave(plot = p_counts0, filename = "figures/RMSE_plot_counts0.pdf", width = 8.51, height = 3.61)

p_counts15 <- df_all %>%
  filter(theta == 100, overd == 0.15) %>%
  dplyr::select(Method, J, K_new, overd, rmse_Lambda, rmse_Counts) %>%
  mutate(J = as.factor(J),
         Lambda = rmse_Lambda,
         Counts = rmse_Counts,
         overd = paste0("Overdispersion = ", as.factor(overd)), 
         K_new2 = paste0("K = ", K_new + 4),
         K_new2 = fct_reorder(K_new2, as.integer(K_new))) %>%
  dplyr::select(-rmse_Counts, -rmse_Lambda, -K_new) %>%
  gather(key = "Metric", value = "value", -Method, -J, -overd, -K_new2) %>%
  ggplot(aes(x = J, y = value, fill = Method, color = Method))+
  geom_boxplot(alpha = 0.4)+
  theme_bw()+
  scale_color_manual( values = color_values, labels = labels)+
  scale_fill_manual(values = color_values, labels = labels)+
  facet_grid(Metric ~ K_new2, scales = "free_y") +
  theme(#aspect.ratio = 1,
    legend.position ="right",
    legend.margin=margin(0,0,0,0),
    legend.box.margin=margin(-5,-5,-5,-5),
    panel.grid.major.x = element_blank(), 
    panel.grid.minor.x = element_blank(), 
  )+
  xlab("Number of samples J")+
  ylab("RMSE")+
  labs(title = "Overdispersion = 0.15")
ggsave(plot = p_counts15, filename = "figures/RMSE_plot_counts15.pdf", width = 8.51, height = 3.61)




