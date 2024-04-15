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

## Codes for reproducing results

* Simulation: `R/main_Simulation.R`

* Application: `R/Application_21Breast.R`


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
  - Our implementation of the shrinkage process proposed by Legramanti et al. (2020), adapted to the Poisson factorizaton

### Files to run the Simulation

- `R/main_Simulation.R`:

- `R/main_Simulation_analysis.R`:


### Files to run the Application

- `R/Application_21Breast.R`

### Additional functions to postprocess the output and to plot the signatures

- `R/Postprocess_functions.R`

- `R/rInvKummer.R`

- `R/plot_signatures.R`


### Miscellanea

- `R/get_COSMIC_data.R`

- `src/tune_betah.cpp`

### Citations
> Wang S, Li H, Song M, Tao Z, Wu T, He Z, et al. (2021) Copy number signature analysis tool and its application in prostate cancer reveals distinct mutational processes and clinical outcomes. PLoS Genetics 17(5): e1009557.
https://doi.org/10.1371/journal.pgen.1009557
> Rosales R A, Drummond R D, Valieris R, Dias-Neto E, da Silva I T. (2017) signeR: an empirical Bayesian approach to mutational signature discovery, Bioinformatics, Volume 33, Issue 1, 8–16. https://doi.org/10.1093/bioinformatics/btw572


