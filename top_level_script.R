###########################################################################################################
## Parasite evolution and host defense evolution, code based on BMB's Nearly Neutral Parasite Evolution, ##
## expanded to deal with multiple axes, a tradeoff curve, and the very beginnings of host evolution      ##  
##                                                                                                       ##
## Optimally this code would be fully integrated with BMB's code (model options) so that the two models  ##
## could be toggled between with a knob; however, I haven't gotten there yet. Only the simplest piece    ##
## of that model can be run from these scripts (beta evolution in discrete time only)                    ##
##                                                                                                       ##
## Most recently semi-cleaned as of the week of June 17, 2019 in prep for transitioning to new work      ##
## All notes contained within all scripts are up to date as of the week of June 17, 2019                 ##
###########################################################################################################

## Notes: 
 ## 1) This repo was created originally for EEID 2019. The poster I presented at that conference is
   ## contained within this repo.
 ## 2) If you have made your way here by way of EEID or stumbled across this repo on your own:
   ## The code here is a work in progress and not super clean / fully prepared for sharing. You should be able to
   ## work your way through by following this script, but errors may be encountered (especially in plotting)

## See NOTE: ** for a couple of key notes in a few scripts

## BMB: what *exactly* does this mean?
host_mort <- FALSE

## load (and install if necessary!) required packages
source("needed_packages.R")

## ggplot theme
source("ggplot_theme.R")

## Load equired functions to run the model.
## The SIR script is poorly named: it is really "SIS with host mortality". 
## funs_SIS.R is SIS with no host mortality
if (host_mort) {
    source("funs_SIR.R")
} else {
    source("funs_SIS.R")   
}

#####
## Options to run the model and store results:
## 1. parameter_setup.R *then* res_tol_sim.R. Best to use for single runs to explore parameter values
## 2. XXXX_manual_setup.R *then* XXXX_tidy.R. Large parameter space exploration and initial cleaning for plotting
## (Can recover what I fit for thesis and poster)

#####
## Option 1
#####

## set up parameters
if (host_mort) {
    source("parameter_setup.R")
} else {
    source("parameter_setup_nm.R")   
}

## run the sim. For manual plotting exploration the output from this script to use is:
 ## res_1000_stochas   : all
 ## res_1000_stochas_s : runs summarized into quantiles
if (host_mort) {
    source("res_tol_sim.R")
} else {
    ## single run
    source("vir_evo_sim.R")
}
  
#####
## Option 2
#####

if (host_mort) {
    ## Set up parameter values
    source("res_tol_stochas_manual_setup.R")
    source("res_tol_determ_manual_setup.R")
    source("res_tol_ad_manual_setup.R")
} else {
    ## This script now has lots of exploration and debugging stuff. Don't run through source() but it
    ## does remain a useful script to look through
    ## source("vir_evo_stochas_manual_setup.R")  
    ## streamlined script for tradeoff only
    source("vir_evo_stochas_to_bulk.R")
    ## streamlined script for efficiency model
    source("vir_evo_stochas_eff_bulk.R")   
}

## Examine the results
if (host_mort) {
    source("stochastic_tidy.R")
    source("determ_tidy.R")
}

## Can also just load previous results... Email me for thesis or poster "raw data" (sim results), too large for github 

## If plotting from previously saved results, need the following scripts for res_tol_plots_final.R to work because:
# (1) res_tol_determ_manual_setup.R: need the parameter values from this script to set up tradeoff curve and adaptive landscape 
# (2) determ_tidy.R: load and clean saved results from 
# (3) stochastic_tidy.R
# (4) res_tol_plots_final.R

## Plot
if (host_mort) {
    source("res_tol_plots_final.R")
} else {
    ## Plot one model at a time
    tradeoff_only <- TRUE
    source("vir_evo_hypercube_plotting.R")
}

 ## NOTE: ** res_tol_plots_final.R is pretty incomplete. Can probably step through this script, but I warn you (*whoever
  ## is finding this note*) that that plotting script is quite convoluted. Either work your way through it or do your
   ## own plotting from the results you fit.

