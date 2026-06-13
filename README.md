# CompressiveNMF

Code to reproduce figures, simulation results, and real data analysis results from the paper 

> Zito, A. and Miller, J. W. (2024) - Compressive Bayesian non-negative matrix factorization for
mutational signatures analysis, arxiv: https://arxiv.org/abs/2404.10974

## As a first step, install the CompressiveNMF package

```
install.packages("CompressiveNMF_1.0.0.tar.gz")
library(CompressiveNMF)
```

## Codes to repreduce the figures:

### Interpretation
* Figure 1 and Figure S12 - Density of the Inverse Kummer and concentration behavior
  - `R/reproduce_Figures_1_S12.R`
  
### Simulation
* Figure 2, and Figure S4 to S9 in the Supplementary material - Main simulation
  - `R/main_Simulation_analysis.R`
* Figure S1 and S2 - Realignement and label switching
  - `R/test_time_CompNMF.R`
* Figure S10 - Computational time for varying M and N
  - `R/Simulation_sparsity_indels_suppl.R`
* Figure S11 - Sparse indel simulation
  - `R/Simulation_sparsity_indels_suppl.R`
* Figure S12 - CompNMF vs Infinite factorization models
  - `R/Simulation_Compressive_vs_InfiniteFactors.R`
* Figure S14, S15 and S16 - Fixed vs Compressive hyperprior
  - `R/Simulation_fixed_vs_compressive.R`
* Figure S17 - Sensitivity to epsilon and K
  - `R/Simulation_sensitivity_epsilon_K.R`
* Figure S18 - Sensitivity to a and alpha
  - `R/Simulation_sensitivity_alpha_a_suppl.R`
* Figure S19 - Effect of hyperprior on alpha
  - `R/comparison_with_randomAlpha.R`
  
### Application

* Figure 3, 4, and Figure S20 to S26, Table S1 in the Supplementary material - 21 breast cancer application
  - `R/Application_21Breast.R`
* Figure S27 and S28 - Sensitivity to hyperparameters in 21 breast cancer
  - `R/Application_21Breast_sensitivity_analysis.R`
* Figure S29 and S30 - Indel Panc-AdenoCA application
  - `R/Application_PancAdeno_indels.R`

## Codes for reproducing simulation and applications

* Simulation: `R/main_Simulation.R`. This file produces the files `output/main_simulation/df_F1.csv` and `output/main_simulation/simulation_output.csv`, which are needed to reproduce the figures. The original simulation was run on 20 cores in parallel on a cluster. To reproduce the plots, run  `R/main_Simulation_analysis.R`

* Application on 21 breast cancer data: `R/Application_21Breast.R`. This file produces the files is `output/Application_21brca` and the figures. To reproduce the figure only, set `rerun = FALSE`.

* Application on indels Panc-AdenoCA: `R/Application_PancAdeno_indels.R`. This file produces the files is `output/Application_ID_panc` and the figures. To reproduce the figure only, set `rerun = FALSE`.

## R code description, file by file

### Methods

- `R/CompressiveNMF.R` - Implementation of the `CompressiveNMF` method. It runs both the unsupervised version, and the semi-supervised one with the informative prior based on the COSMIC signatures. 

- `R/SignatureAnalyzer.R` - Wrapper for the function `sig_auto_extract` from the `R` package `sigminer` (Wang et al., 2021) 

- `R/signeR.R` - Wrapper for the function `signeR` from the `R` package `signeR` (Rosales et al., 2016) 

- `R/SigProfilerExtractor.R` - Wrapper for the function `sigprofiler_extract` from the `R` package `sigminer` (Wang et al., 2021) 

- `R/Poisson_InfiniteFactors.R` - Our implementation of the shrinkage process proposed by Legramanti et al. (2020), Kowal and Canale (2023) and Batthacharya and Dunson (2011), adapted to Poisson factorization 

- Folder `R/run_BayesNMF_python` - All the python scripts to run the BayesNMD ARD (Brouwer et al. (2017))

- Folder `R/run_SigProfiler_python` - All the python scripts to run SigProfiler extractor for indels.

### Additional functions to postprocess the output and to plot the signatures

- `R/Postprocess_functions.R` - Functions to post-process the output from every method. 

- `R/rInvKummer.R` - Function to evaluate the density and generate random samples from the inverse Kummer distribution

- `R/plot_signatures.R` - Functions to plot the cosmic signature. For example,`plot_cosmic_signature("SBS5")` plots signature SBS5.  

- `R/tune_Cosmic_priors.R` - Functions to tune the hyperparamers `betah` in the `CompressiveNMF` (informative priors) 

### Miscellanea

- `R/get_COSMIC_data.R` - Function to scrape the COSMIC v3.4 signature data and the proposed aetiologies from https://cancer.sanger.ac.uk/signatures/sbs/ 

- `src/tune_betah.cpp`  - Function to tune the hyperparameter `betah` in the `CompressiveNMF` method when specifying informative priors. See the file `data/optimal_betah.rdata`, and refer to the Supplementary material.


#### References

> Brouwer, T., Frellsen, J., and Lio, P. (2017). Comparative Study of Inference Methods for Bayesian Nonnegative Matrix Factorisation. In Proceedings of the European Conference on Machine Learning and Principles and Practice of Knowledge Discovery in Databases
https://link.springer.com/chapter/10.1007/978-3-319-71249-9_31

> Wang S., Li H., Song M., Tao Z., Wu T., He Z., et al. (2021) Copy number signature analysis tool and its application in prostate cancer reveals distinct mutational processes and clinical outcomes. PLoS Genetics 17(5): e1009557.
https://doi.org/10.1371/journal.pgen.1009557

> Rosales R. A., Drummond R. D., Valieris R., Dias-Neto E., da Silva I. T. (2017) signeR: an empirical Bayesian approach to mutational signature discovery, Bioinformatics, Volume 33, Issue 1, 8–16. https://doi.org/10.1093/bioinformatics/btw572

> Legramanti S., Durante D., Dunson D. B. (2020) Bayesian cumulative shrinkage for infinite factorizations, Biometrika, Volume 107, Issue 3, 745–752. https://doi.org/10.1093/biomet/asaa008

> Kowal, D. and Canale, A. (2023) Semiparametric Functional Factor Models with Bayesian Rank Selection, Bayesian Analysis, 10.1214/23-BA1410

> Bhattacharya, A. and Dunson D. B. (2011) Sparse Bayesian infinite factor models, Biometrika, Volume 98, Issue 2, Pages 291–306, https://doi.org/10.1093/biomet/asr013
