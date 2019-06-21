##########################################
## Load and clean up deterministic runs ## 
##########################################

## Script should be run in conjunction with stochastic_tidy.R

## First grab the AD solutions, some place you previously saved the file, or in the workspace...
#grad_ascent_tot      <- readRDS("res_out/res_out_AD/AD_grad_ascent.Rds")
#grad_ascent_sing_tot <- readRDS("res_out/res_out_AD/AD_grad_ascent_sing.Rds")

## If you are saving a large number of determnistic runs (RD), because of how the output from the
 ## deterministic runs are saved, the Rds files need to be loaded after res_tol_determ_manual_setup.R
  ## is run with the parameters used to run the sim that generated that Rds file (and other is fine, 
   ## but the parameters need to be contained)
  
## To match stochastic results to deterministic results, search for parameter values for a given
 ## stochastic run in the deterministic parater value set, then load the appropriate saved Rds

## start with a given stochastic run, find a deterministic run with comperable parameters
unique(stochas.res_subset$mut_sd)

####################
## For plotting a single parameter combination on the adaptive landscape need to subset
## down to a single parameter combination. For the other plots, use the median
## values for all parameters except all values for a single parameters of focus
####################
params %>% 
  filter(
  eff_scale  == 30
, Imat_seed1 == 72
, Imat_seed2 == 28) %>%
  dplyr::select(param_num)
  
## and load that row
  # determ.res_subset <- readRDS("res_out/res_out_RD/res_1000_all_determ_6.Rds")
## Temp R0 column for plotting on the surface. May actually want to plot calculated R0 later
determ.res_subset  <- determ.res_subset[[1]]
## NOTE: ** Eeek! The deterministic data frame has the names flipped!!!
 ## need to change by_row in mat_out in res_tol_determ_manual_setup.R I believe
names(determ.res_subset)[c(3, 4, 5)] <- c("R0", "tune", "agg")

## What is recorded in this data frame is actual proportion of I from the deterministic run (sum = 1 - S)
 ## Need to scale density to calculate mean and sd
determ.res_subset <- determ.res_subset %>% 
  group_by(Time) %>%
  mutate(R0 = R0 / sum(R0))

## also subset the adaptive dynamics results
grad_ascent.res_subset <- grad_ascent_tot %>%
  filter(
    eff_scale  == 30
  , Imat_seed1 == 87
  , Imat_seed2 == 13)

grad_ascent_sing.res_subset <- grad_ascent_sing_tot %>%
  filter(
    eff_scale  == 30
  , Imat_seed1 == 87
  , Imat_seed2 == 13)
