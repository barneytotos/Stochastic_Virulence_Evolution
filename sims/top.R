### Batch runs for paper

## FALSE runs SIS, TRUE is not supported right now
host_mort <- FALSE

## load (and install if necessary!) required packages
source("needed_packages.R")

## ggplot theme
source("ggplot_theme.R")

## funs
source("hpevosim/R/funs_SIS.R")  

## model options:
 ## "nt"  = no tradeoff
 ## "to"  = tradeoff only
 ## "eff" = efficiency
 ## [later] "tune" = tuning
model.choice <- "eff"

## determinisitc or stochastic
deterministic <- TRUE

## Note:  
 ## for hypercube set: num_points to some large value and num_runs to 1
 ## for stochasticity at one combo set: num_points to some small value and num_runs to a larger value
num_runs     <- 1

## Number of parameter combinations: 1000 for other models 2000 here because of TRUE, FALSE above
num_points <-  
 if (model.choice == "nt") {
    2000
 } else {
    1000
 }

## placing in to top for ease of debugging...
num_points <- 10

## load params. ifelse for deterministic inside 
source("sims/params.R")

## run sim. Requires different setup so separate deterministic outside
if (deterministic) {
 source("sims/run_d.R") 
} else {
 source("sims/run_s.R")
}

