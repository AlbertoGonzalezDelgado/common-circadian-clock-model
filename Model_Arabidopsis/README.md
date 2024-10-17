# Arabidopsis FT simulation in LD and SD photoperiods
This directory contains the files and scripts needed to model FT expression in Arabidopsis using ordinary differential equations (ODE) that receive as input the temporal expression of circadian clock genes ([circadian_clock_info](data/Circadian_clock.csv) and [circadian_clock_expression](data/circadian_clock_fit.tsv) files). 


* [Arabidopsis_model.m](Arabidopsis_model.m) script applies simulated annealing and grid search algorithms to obtain the best parameters to apply in the equation to fit the FT model to LD photoperiod according to a cost function. Usage: Arabidopsis_model.m <circadian_clock_expression_file> <circadian_clock_info_file>

* [print_results.m](print_results.m) script generates a plot of FT simulations using best parameters for LD and SD photoperiod. Usage: print_result.m <circadian_clock_expression_file> <circadian_clock_info_file>
   
