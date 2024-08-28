# Repository for the paper titled "Reviving collapsed plant-pollinator networks from a single species."
Authors: Gaurav Baruah, Meike Wittmann
For Journal: PLOS Biology
### This repository contains the following Rscripts, and data used for producing the figures in the main-text as well as the for the supplementary data (in folder appendix_data for reproducing the figures of the supplementary).

The codes and scripts are prepared by Gaurav Baruah.

#### Simulated R Data and Empirical data from web-of-life database and Kaiser-Bunbury et al 2017 paper: 
1. The folder `datasets_1` contains incidence data for the plant-pollinator networks used in the manuscript.

2. `figure_4_network_data.RData`, `figure_4_species_level_data.RData`,are used to reproduce the figure 4 used in the main-text. These data has been generated through eco-evolutionary dynamic simulations of 115 plant-pollinator network as described in the main-text. `Mean_trait_data.RData` is the mean trait values of species of all the 115 networks that are used for network ressurection analysis, or section 4.3 (Perturbation regime) of the main-text (used in figure 2, figure 3, figure 4). `Upper_bound_abundance_response_time.RData` is a miscellaneous data that was used for miscellaneous results (not used in main-text).
3. `hysteresis_net_data.RData` and `hysteresis_species_data.RData` are the data used to produce figure 1. `hysteresis_species_data.RData` is the species level data that was used to produce figure 1D.
4. `Empirical_data.RData` is the data extracted from Kaiser-Bunbury et al 2017 paper which can be obtained from the open website [repository] (http://www.ecologia.ib.usp.br/iwdb/html/
kaiser-bunbury_et_al_2017.html) and was used to produce figure 1E-1F.
5. `03_figure_data.RData` simulated summary data used to produce the figure 3 and used in R script `fig3a.R` and `fig3b.R`.
6. `Mean_trait_data.RData` is the quasi- eco-evolutionray equilibria mean trait data that are subsequently used to simulate whether networks could be revived or not. Used in `fig2.R` and for producing figure 3.

#### R scripts:

1. `01_ODE_Function.R` is the script that details all functions required for the eco-evolutionary dynamical simulation and analysis and figure plotting.
2. `01_hyst_functions.R` is the R script that has all the functions required to produce figure 1 .
3. `fig1.R` R script to reproduce figure 1.
4. `fig2.R` R script that does the eco-evolutionary simulation to reproduce figure 2.
5. `fig3a.R` R script that does the eco-evolutionary simulation to reproduce part of figure 3 (i.e., Figure 3 A-H)
6. `fig3b.R`  R script that does the eco-evolutionary simulation to reproduce second part of figure 3 (i.e., Figure 3 I-N)
7. `fig4_5.R` R script to reproduce figure 4 and figure 5. 
