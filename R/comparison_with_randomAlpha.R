library(CompressiveNMF)
library(tidyverse)
library(coda)
library(foreach)
library(doParallel)
devtools::document()
data <- readRDS("~/CompressiveNMF/output/main_simulation/Scenario_50_overdisp_0_Knew_2_theta_100/data.rds.gzip")
data_mispec <- readRDS("~/CompressiveNMF/output/main_simulation/Scenario_50_overdisp_0.15_Knew_2_theta_100/data.rds.gzip")

burnin <- 4000
nsamples <- 1000

registerDoParallel(20)
set.seed(42, kind = "L'Ecuyer-CMRG")

# Correctly specified case
res <- foreach(i = 1:20) %dopar% {
  X <- data[[i]]$X

  tryCatch(
    {
      CompressiveNMF_randa(
        X,
        random_a = TRUE,
        nsamples = nsamples,
        burnin  = burnin
      )
    },
    error = function(e) {
      message(sprintf("Error in iteration %d: %s", i, e$message))
      return(e$message)  # or some placeholder object
    }
  )
}
out <- lapply(1:20, function(i) Postprocess_Compressive(resComp = res[[i]], data = data[[i]]))

# Mispecified case
res_mispec <- foreach(i = 1:20) %dopar% {
  X <- data_mispec[[i]]$X
  
  tryCatch(
    {
      CompressiveNMF_randa(
        X,
        random_a = TRUE,
        nsamples = nsamples,
        burnin  = burnin
      )
    },
    error = function(e) {
      message(sprintf("Error in iteration %d: %s", i, e$message))
      return(e$message)  # or some placeholder object
    }
  )
}
out_mispec <- lapply(1:20, function(i) Postprocess_Compressive(resComp = res_mispec[[i]], 
                                                               data = data_mispec[[i]]))


df1 <- as.data.frame(t(sapply(1:20, function(i) out[[i]]$results, simplify = TRUE)))
df1$J <- 50
df1$overd <- 0
df2 <- as.data.frame(t(sapply(1:20, function(i) out_mispec[[i]]$results, simplify = TRUE)))
df2$J <- 50
df2$overd <- 0.15
df1$a <- df2$a <- "Random a"

df_all <- read_csv(file = "~/CompressiveNMF/output/main_simulation/simulation_output.csv")
df_fixed <- df_all %>%
  filter(Method == "2.CompNMF", J == 50, K_new == 2)
df_fixed$a <- "Fixed a"

df_join <- bind_rows(df_fixed, df1, df2)


pK <- ggplot(df_join, aes(x = as.factor(overd), y = K, color = a)) +
  geom_hline(yintercept = 6, linetype = "dashed") +
  geom_boxplot(outlier.shape = NA) +
  theme_bw() +
  geom_point(position = position_jitterdodge(jitter.width = 0.1, 
                                             jitter.height = 0, 
                                             dodge.width = 0.75))+
  ylim(c(4, 18))+
  xlab("Overdispersion")+
  theme(legend.position = "top")+
  facet_grid(~"Estimated K")
  

a_list <- lapply(1:20, function(i) {
  data.frame(
    run = as.factor(i),                        # Make 'i' a factor for grouping/coloring
    a_value = res[[i]]$mcmc_out[[1]]$a         # Extract the 'a' vector
  )
})

a_list_mispec <- lapply(1:20, function(i) {
  data.frame(
    run = as.factor(i),                        # Make 'i' a factor for grouping/coloring
    a_value = res_mispec[[i]]$mcmc_out[[1]]$a         # Extract the 'a' vector
  )
})


# Combine them all into one long data frame
df_a <- bind_rows(a_list) %>%
  left_join(df1 %>% mutate(run = as.factor(1:20)), by = "run")

pa1 <- ggplot(df_a, aes(x = a_value, color = as.factor(K), fill = run)) +
  geom_density(alpha = 0.5) +
  theme_bw() +
  scale_fill_manual(values = rep(NA, 20)) +
  guides(fill = "none") +
  scale_color_viridis_d(name = "Estimated\nK", option = "mako", begin = 0.2, end = 0.8)+
  labs(x = "a", y = "Density") +
  facet_grid(~"Overdispersion = 0")

df_a2 <- bind_rows(a_list_mispec)%>%
  left_join(df2 %>% mutate(run = as.factor(1:20)), by = "run")

pa2 <- ggplot(df_a2, aes(x = a_value, color = as.factor(K), fill = run)) +
  geom_density(alpha = 0.5) +
  theme_bw() +
  scale_fill_manual(values = rep(NA, 20)) +
  scale_color_viridis_d(name = "Estimated\nK", option = "mako", begin = 0.2, end = 0.8)+
  guides(fill = "none") +
  labs(x = "a", y = "Density") +
  facet_grid(~"Overdispersion = 0.15")

pK + pa1 + pa2


