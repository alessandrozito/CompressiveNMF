# This function runs sigprofiler extractor using the package called
# sigminer
library(sigminer)
library(reticulate)

# Note: the matrix has to be in the format of J x 96 format, columns are signatures
sigprofiler <- function(X, out_dir_sigpro = "sigprofiler_out/", Kmin = 2, Kmax = 8, cores = 1, is_exome = FALSE, ...){
  # Record computational time
  t_start = Sys.time()
  # Configurate python to run SigProfiler extractor
  use_python(reticulate::conda_list()[2, 2])
  py_config()
  
  # Run sigprofiler extractor
  try(sigminer::sigprofiler_extract(nmf_matrix = t(X),
                                output = out_dir_sigpro, 
                                range = Kmin:Kmax, init_method = "random", 
                                py_path = reticulate::conda_list()[2, 2], 
                                refit_plot = FALSE, refit = FALSE, nrun = 5L,
                                sigprofiler_version = "1.1.23",is_exome = is_exome,
                                cores = cores), 
      outFile = paste0(out_dir,"try_out.txt"))
  
  # Take now the suggested solution and save it
  
  # Signatures
  solution_sign <- paste0(out_dir_sigpro, "/SBS96/Suggested_Solution/SBS96_De-Novo_Solution/Signatures/SBS96_De-Novo_Signatures.txt")
  Signatures <- suppressMessages(readr::read_table(file = solution_sign))
  channels <- Signatures$MutationsType
  Signatures <- as.matrix(Signatures[, -1])
  rownames(Signatures) <- channels
  
  # Weights
  solution_weights <- paste0(out_dir_sigpro, "/SBS96/Suggested_Solution/SBS96_De-Novo_Solution/Activities/SBS96_De-Novo_Activities.txt")  
  Weights <- suppressMessages(readr::read_table(file = solution_weights))  
  Weights <- t(as.matrix(Weights[, -1]))
  
  time <- difftime(Sys.time(), t_start, units = "secs")[[1]]
  
  
  return(list("Signatures" = Signatures,
              "Weights" = Weights,
              "K" = ncol(Signatures),
              "time" = time))
  
}

