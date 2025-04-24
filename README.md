# Resilience Animal Networks

This repository includes all data and code necessary to reproduce the paper "Resilience of animal societies to anthropogenic removals is higher when key network positions are protected: evidence from three long-term data sets".

To reproduce the paper run the scripts in the following order: 

- Generate the networks with `generate_networks.R`
- Calculate differences in efficiency with `efficiency_after_removal.R`
- Then run the simulations with `run_all_analyses.R` (This requires a lot of computational power and might take quite a long time).
- `plot_script.R` includes the script for all the plots. 
