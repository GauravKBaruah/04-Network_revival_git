# Repository for the paper titled "Reviving collapsed networks from a single species: the importance of trait variation and network architecture."
Authors: Gaurav Baruah, Meike Wittmann

### This repository contains the following Rscripts, and data used for producing the figures

The codes and scripts are prepared by Gaurav Baruah.

#### Data: 
The folder `datasets_1` contains data for all plant-pollinator networks used in the paper.
`figure_4_network_data.RData`, `figure_4_species_level_data.RData`,are used to reproduce the figure 4 used in the main-text. These data has been generated through eco-evolutionary dynamic simulations of 115 plant-pollinator network as described in the main-text. `Mean_trait_data.RData` is the mean trait values of species of all the 115 networks that are used for network ressurection analysis, or section 4.3 (Perturbation regime) of the main-text (used in figure 2, figure 3, figure 4). `Upper_bound_abundance_response_time.RData` is a miscellaneous data that was used for miscellaneous results (not used in main-text). `hysteresis.RData` is the data used to produce figure 1 and figure S1, contains the data used in producing the figure 1.

#### R scripts:

`01_ODE_Function.R` is the script that details all functions required for the eco-evolutionary dynamical simulation and analysis and figure plotting.
`01_hyst_functions.R` is the R script that has all the functions required to produce figure 1 .
`01_Figure_1_hysteresis.R` R script to reproduce figure 1.
`figure_2.R` R script that does the eco-evolutionary simulation to reproduce figure 2.
`figure_3_a.R` R script that does the eco-evolutionary simulation to reproduce part of figure 3 (i.e., Figure 3 A-H)
`figure_3_b.R`  R script that does the eco-evolutionary simulation to reproduce second part of figure 3 (i.e., Figure 3 I-N)
`figure_4.R` R script to reproduce figure 4.
