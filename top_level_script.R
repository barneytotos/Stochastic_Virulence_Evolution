###########################################################################################################
## Parasite evolution and host defense evolution, code based on BMB's Nearly Neutral Parasite Evolution, ##
## expanded to deal with multiple axes, a tradeoff curve, and host evolution                             ##
###########################################################################################################

## Note to any how may have found their way here by way of EEID or stumbled across this repo on their own:
 ## The code here is a work in progress and not super clean / fully prepared for sharing. You should be able to
 ## work your way through by following this script, but errors may be encountered (especially in plotting)

## required packages
source("needed_packages.R")

## ggplot theme
source("ggplot_theme.R")

## required functions to run the model. **Note: This script is poorly named. The model is an SIS model.
source("funs_SIR.R")

#####
## Options to run the model and obtain results include
#####
# (1) parameter_setup.R _then_ res_tol_sim.R. Best to use for single runs to explore parameter values
# (2) XXXX_manual_setup.R _then_ XXXX_tidy.R. Large parameter space exploration and initial cleaning for plotting

#####
## Option 1
#####

## set up parameters
source("parameter_setup.R")

## run the sim
source("res_tol_sim.R")
  
#####
## Option 2
#####

## Set up parameter values
source("res_tol_stochas_manual_setup.R")
source("res_tol_determ_manual_setup.R")
source("res_tol_ad_manual_setup.R")

## Examine the results
source("stochastic_tidy.R")
source("determ_tidy.R")
source("res_tol_plots_final.R")

### Can also just load previous results
# res_1000_all <- readRDS("res_out/res_1000_all.Rds")


## If plotting from previously saved results, need the following scripts for res_tol_plots_final.R to work because:
# (1) res_tol_determ_manual_setup.R: need the parameter values from this script to set up tradeoff curve and adaptive landscape
# (2) determ_tidy.R: load and clean saved results from 
# (3) stochastic_tidy.R
# (4) res_tol_plots_final.R

#####
## Scripts that aren't really needed
#####
