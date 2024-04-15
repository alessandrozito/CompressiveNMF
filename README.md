# CompressiveNMF

Code to reproduce figures, simulation results, and real data analysis results from the paper 

> Zito, A. and Miller, J. W. (2024) - Compressive Bayesian non-negative matrix factorization for
mutational signatures analysis, arxiv:

## Codes to repreduce the figures:

* Figure 1 - Density of the Inverse Kummer and concentration behavior
  - `R/rInvKummer.R`
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
  - The file contains our

- `R/SignatureAnalyzer.R`

- `R/signeR.R`

- `R/SigProfilerExtractor.R`

- `R/PoissonCUSP.R`

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



