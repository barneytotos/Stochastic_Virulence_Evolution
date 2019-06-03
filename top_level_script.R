###########################################################################################################
## Parasite evolution and host defense evolution, code based on BMB's Nearly Neutral Parasite Evolution, ##
## expanded to deal with multiple axes, a tradeoff curve, and host evolution                             ##
###########################################################################################################

### required packages
source("needed_packages.R")

### ggplot theme
source("ggplot_theme.R")

### required functions to run the model
source("funs_SIR.R")

### Load previous results?
load_previous <- FALSE

if (load_previous == FALSE) {
### set up parameters
source("parameter_setup.R")

### run the sim
source("res_tol_sim.R")
  
} else {
### load previous results
res_1000_all <- readRDS("res_out/res_1000_all2b.Rds")
}

### If plotting from previously saved results, need the following scripts:
 ## 1) res_tol_determ_manual_setup.R
 ## 2) determ_tidy.R
 ## 3) stochastic_tidy.R
 ## 4) res_tol_plots_final.R

### plot some results
source("res_tol_plots.R")

#####
## Manual exploration of runs across a large parameter set
#####

### Set up parameter values
source("res_tol_stochas_manual_setup.R")
source("res_tol_determ_manual_setup.R")
source("res_tol_ad_manual_setup.R")

### Examine the results
source("stochastic_tidy.R")
source("determ_tidy.R")
source("res_tol_plots_final.R")

