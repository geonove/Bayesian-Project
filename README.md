# Bayesian Inference for Dynamic Graphical Models

This repository contains the following files to run a Bayesian inference algorithm for Dynamic Graphical Models:

- `HMM_DAGS.R`, which contains the code necessary to run the final model on a simulated dataset
- `real_data.R`, which contains the code necessary to run the final model on a real dataset
- `dataset.xlsx`, which is the real dataset used in `real_data.R`, taken from the file `FinancialMarketData.xlsx`
- `GS_DAGS.R`, which is a simplified version of `HMM_DAGS.R` and considers a mixture model instead of the more general hidden Markov model (HMM)

Additionally, this repository contains three folders:

- `dags` folder, which contains the necessary functions to deal with dags
- `gibbs` and `hmm` folders, which contain code that was used in the development of the final codes in `GS_DAGS.R` and `HMM_DAGS.R`
