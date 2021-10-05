# Microbial Community Coalescence
## Introduction
This repository contains the code to perform simualtions and analysis for Lechón et al. The role of competition versus cooperation in microbial community coalescence (2021). 

The repository contains three folders: `code` which contains the code to generate communities, run simualtions and plot figures, `data` which is used to store simulation results and `sandbox` which contains...

## Running simulations

In order to run the simulations in the manuscript you'll need to run two scripts. First run `code/assembly_clean.py` which generates parent communities, simulates till equilibrium and saves the output in the `data` folder. To run the coalescence simualtions run `code/coalescence_clean.py` which will simulate community coalescence and save the output to `data`.

## Plotting results

To plot the simualtion results as seen in the manuscript run any of the respective `R` scripts in the `code` folder:
```
├── figure_2.R
├── figure_3.R
├── figure_4.R
├── figure_S3.R
├── figure_S4.R
├── figure_S6.R
```
which will save the figures in the `sandbox` folder. 

