# Add effective sample size to the results

# Load libraries
library(tidyverse)
library(CompressiveNMF)
library(coda)
library(lsa)

source("R/Postprocess_functions.R")

# Useful function
create_directory <- function(dir) {
  if (!dir.exists(dir)) {
    dir.create(dir)
  }
}

simulation_dir <- "~/CompressiveNMF/output/main_simulation/"

overdispersion_list <- c(0, 0.15)
J_list <- c(50, 100, 200)
K_new_list <- c(2, 6)
theta <- 100

# overdispersion <- 0
# J <- 50
# K_new <- 2

for(K_new in K_new_list){
  for(J in J_list) {
    for(overdispersion in overdispersion_list){
      out_dir <- paste0(simulation_dir, "Scenario_", J, "_overdisp_", 
                        overdispersion, "_Knew_", K_new, "_theta_", theta)
      data_all <- readRDS(paste0(out_dir, "/data.rds.gzip"))
      print(c(K_new, overdispersion, J))
      #------------------------------------------------- CompressiveNMF 
      print("CompNMF")
      resComp <- lapply(1:20, function(i){
        # Read the result
        nm <- paste0(out_dir, "/CompressiveNMF/results_", i, ".rds.gzip")
        out <- readRDS(nm)
        Postprocess_Compressive(out, data_all[[i]])
      })
      # Save the merged output
      saveRDS(resComp, file = paste0(out_dir, "/CompressiveNMF.rds.gzip"))
      
      #-------------------------------------------------- CompressiveNMF + cosmic
      print("CompNMFcosmic")
      resComp_cos <- lapply(1:20, function(i){
        # Read the result
        nm <- paste0(out_dir, "/CompressiveNMF_cosmic/results_", i, ".rds.gzip")
        out <- readRDS(nm)
        Postprocess_Compressive(out, data_all[[i]])
      })
      # Save the merged output
      saveRDS(resComp_cos, file = paste0(out_dir, "/CompressiveNMF_cosmic.rds.gzip"))
      
      #-------------------------------------------------- signeR
      print("signeR")
      resSigneR <- lapply(1:20, function(i){
        # Read the result
        nm <- paste0(out_dir, "/signeR/results_", i, ".rds.gzip")
        out <- readRDS(nm)
        Postprocess_signeR(out, data_all[[i]])
      })
      # Save the merged output
      saveRDS(resSigneR, file = paste0(out_dir, "/signeR.rds.gzip"))
      
      # #-------------------------------------------------- PoissonCUSP
      # print("PoissonCUSP")
      # resCUSP <- lapply(1:20, function(i){
      #   # Read the result
      #   nm <- paste0(out_dir, "/PoissonCUSP/results_", i, ".rds.gzip")
      #   out <- readRDS(nm)
      #   Postprocess_PoissonCUSP(out, data_all[[i]])
      # })
      # # Save the merged output
      # saveRDS(resCUSP, file = paste0(out_dir, "/PoissonCUSP.rds.gzip"))
      
    }
  }
}






