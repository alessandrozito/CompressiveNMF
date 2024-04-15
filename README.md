# CompressiveNMF

Code to reproduce figures, simulation results, and real data analysis results from the paper 

> Zito, A. and Miller, J. W. (2024) - Compressive Bayesian non-negative matrix factorization for
mutational signatures analysis, arxiv:

## Codes to repreduce the figures:

* Figure 1 - Density of the Inverse Kummer and concentration behavior
  - `R/reproduce_Figure1.R`
* Figure 2, 3 and Figure S1 - S4 in the Supplementary material
  - `R/main_Simulation_analysis.R`
* Figure 4, 5 and Figure S5 - S11 in the Supplementary material
  - `R/Application_21Breast.R`

## Codes for reproducing simulation and applications

* Simulation: `R/main_Simulation.R`. This file produces the files `output/main_simulation/df_F1.csv` and `output/main_simulation/simulation_output.csv`, which are needed to reproduce the figures. The original simulation was run on 20 cores in parallel on a cluster. To reproduce the plots, run  `R/main_Simulation_analysis.R`

* Application on 21 breast cancer data: `R/Application_21Breast.R`. This file produces the files is `output/Application_21brca` and the figures. To reproduce the figure only, simply set `rerun = FALSE`.

## R code description, file by file

### Methods

- `R/CompressiveNMF.R`
  - Implementation of the `CompressiveNMF` method. It runs both the unsupervised version, and the semi-supervised one with the informative prior based on the COSMIC signatures. 

- `R/SignatureAnalyzer.R`
  - Wrapper for the function `sig_auto_extract` from the `R` package `sigminer` (Wang et al., 2021) 

- `R/signeR.R`
  - Wrapper for the function `signeR` from the `R` package `signeR` (Rosales et al., 2016) 

- `R/SigProfilerExtractor.R`
  - Wrapper for the function `sigprofiler_extract` from the `R` package `sigminer` (Wang et al., 2021) 

- `R/PoissonCUSP.R`
  - Our implementation of the shrinkage process proposed by Legramanti et al. (2020), adapted to the Poisson factorizaton. 


### Additional functions to postprocess the output and to plot the signatures

- `R/Postprocess_functions.R`
  - Functions to post-process the output from every method. 

- `R/rInvKummer.R`
  - Function to evaluate the density and generate random samples from the inverse Kummer distribution

- `R/plot_signatures.R`
  - Functions to plot the cosmic signature. For example,`plot_cosmic_signature("SBS5")` plots signature SBS5.  


### Miscellanea

- `R/get_COSMIC_data.R`
  - Function to scrape the COSMIC v3.4 signature data and the proposed aetiologies from https://cancer.sanger.ac.uk/signatures/sbs/ 

- `src/tune_betah.cpp` 
  - Function to tune the hyperparameter `betah` in the `CompressiveNMF` method when specifying informative priors. See the file `data/optimal_betah.rdata`, and refer to the Supplementary material.

#### References
> Wang S., Li H., Song M., Tao Z., Wu T., He Z., et al. (2021) Copy number signature analysis tool and its application in prostate cancer reveals distinct mutational processes and clinical outcomes. PLoS Genetics 17(5): e1009557.
https://doi.org/10.1371/journal.pgen.1009557

> Rosales R. A., Drummond R. D., Valieris R., Dias-Neto E., da Silva I. T. (2017) signeR: an empirical Bayesian approach to mutational signature discovery, Bioinformatics, Volume 33, Issue 1, 8–16. https://doi.org/10.1093/bioinformatics/btw572

> Legramanti S., Durante D., Dunson D. B. (2020) Bayesian cumulative shrinkage for infinite factorizations, Biometrika, Volume 107, Issue 3, 745–752. https://doi.org/10.1093/biomet/asaa008

