# Repository for the paper titled "Reviving collapsed plant-pollinator networks from a single species."
Authors: Gaurav Baruah, Meike Wittmann

### This repository contains the following Rscripts, and data used for producing the figures in the main-text as well as the for the supplementary data (in folder appendix_data for reproducing the figures of the supplementary).

The codes and scripts are prepared by Gaurav Baruah.

#### Data: 
1. The folder `datasets_1` contains data for all plant-pollinator networks used in the paper.

2. `figure_4_network_data.RData`, `figure_4_species_level_data.RData`,are used to reproduce the figure 4 used in the main-text. These data has been generated through eco-evolutionary dynamic simulations of 115 plant-pollinator network as described in the main-text. `Mean_trait_data.RData` is the mean trait values of species of all the 115 networks that are used for network ressurection analysis, or section 4.3 (Perturbation regime) of the main-text (used in figure 2, figure 3, figure 4). `Upper_bound_abundance_response_time.RData` is a miscellaneous data that was used for miscellaneous results (not used in main-text).
3. `hysteresis.RData` is the data used to produce figure 1 and figure S1, contains the data used in producing the figure 1. `hysteresis_species_data.RData` is the species level data that was used to produce figure 1D.
4. `Empirical_data.RData` is the data extracted from Kaiser-Bunbury et al 2017 paper which can be obtained from the open website [repository] (http://www.ecologia.ib.usp.br/iwdb/html/
kaiser-bunbury_et_al_2017.html) and was used to produce figure 1E-1F.

#### R scripts:

1. `01_ODE_Function.R` is the script that details all functions required for the eco-evolutionary dynamical simulation and analysis and figure plotting.
2. `01_hyst_functions.R` is the R script that has all the functions required to produce figure 1 .
3. `01_Figure_1_hysteresis.R` R script to reproduce figure 1.
4. `figure_2.R` R script that does the eco-evolutionary simulation to reproduce figure 2.
5. `figure_3_a.R` R script that does the eco-evolutionary simulation to reproduce part of figure 3 (i.e., Figure 3 A-H)
6. `figure_3_b.R`  R script that does the eco-evolutionary simulation to reproduce second part of figure 3 (i.e., Figure 3 I-N)
7. `figure_4.R` R script to reproduce figure 4 and figure 5. 
