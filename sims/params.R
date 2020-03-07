if (model.choice == "nt") {
   
no_tradeoff      <- TRUE
parasite_tuning  <- FALSE
tradeoff_only    <- FALSE
agg_eff_adjust   <- FALSE
     
} else if (model.choice == "to") {

no_tradeoff      <- FALSE
parasite_tuning  <- FALSE
tradeoff_only    <- TRUE
agg_eff_adjust   <- FALSE
   
} else if (model.choice == "eff") {
   
no_tradeoff      <- FALSE
parasite_tuning  <- FALSE
## tradeoff_only == FALSE and parasite_tuning == FALSE puts model to "efficiency"
tradeoff_only    <- FALSE
## Could have efficiency decrease with an increase in recovery rate, but for now ignore this
 ## have draws for these two traits be independent 
agg_eff_adjust   <- FALSE
   
}

## Generate way more than could possibly be needed
if (!file.exists("lhs_samps_eff.csv")) {
lhs     <- randomLHS(100000, 10) 
write.csv(lhs, "lhs_samps_eff.csv")
} else {
lhs     <- read.csv("lhs_samps_eff.csv")
lhs     <- lhs[, -1]
}


## Lots of parameter values need to be different for deterministic vs stochastic so repeat this whole 
 ## section twice instead of layers of ifelse inside of assignments
if (!deterministic) {

params.all  <- data.frame(
## Parameters to establish what model is run
   no_tradeoff           = no_tradeoff
 , parasite_tuning       = parasite_tuning
 , tradeoff_only         = tradeoff_only
 , agg_eff_adjust        = agg_eff_adjust

## If the traits evolving don't start in different places can have a different script for each of these
 ## Completely ignored if model isn't "nt"
 , nt_mut_var_pos_trait  = ifelse(round(qunif(lhs[, 1], min = 1, max = 2)) == 1, T, F)
   
## Unclear about starting with a range for the trait that is evolving
 , pos_trait0            =
   if (model.choice == "nt") {
      qunif(lhs[, 2], min = 0.01, max = 0.50)
 } else if (model.choice == "to") {
      ## ignored in 'to', transmission is just a function of gamma (neg_trait)
      qunif(lhs[, 2], min = 0.01, max = 0.50)                            
 } else if (model.choice == "eff") {
      ## relative beta in 'eff'. Give full range, but I think there will be relatively few parameter combos that don't lead to extinction
       ## see parameter below
      qunif(lhs[, 2], min = 0.01, max = 0.99)
 }
 , neg_trait0            =
   if (model.choice == "nt") {
      qunif(lhs[, 3], min = 0.01, max = 0.50)
 } else if (model.choice == "to") {
      ## full range allowed here, though may lead to a lot of extinctions...
      qunif(lhs[, 3], min = 0.01, max = 0.99)                            
 } else if (model.choice == "eff") {
      ## full range allowed here, though may lead to a lot of extinctions...
      qunif(lhs[, 3], min = 0.01, max = 0.99)        
 }
   
## Changes by model
 , mut_mean              = 
   if (model.choice == "nt") {
    qunif(lhs[, 4], min = -0.50, max = 0.0)
 } else if (model.choice == "to") {
   ## for now assume no bias in tradeoff only or efficiency model
    0                          
 } else if (model.choice == "eff") {
    0      
 }
   
## Parameters that are the same for all stochastic models
 , mu                    = qunif(lhs[, 5], min = 0.001, max = 0.10)
 , mut_sd                = qunif(lhs[, 6], min = 0.01, max = 0.50) 
 , N                     = round(qunif(lhs[, 7], min = 100, max = 2500))
 , gamma0                = qunif(lhs[, 8], min = 0.01, max = 0.4) 
   
## used only for the tradeoff models but include them. Still a bit unclear of the best parameter values here
 , power_c             = qunif(lhs[, 9], min = 0.0005, max = 0.001)
 , power_exp           = qunif(lhs[, 10], min = 1.5, max = 5.5) 
 
## Extra required parameters that are the same for all models
 , R0_init               = 2
 , deterministic         = deterministic
 , nt                    = 2e5
## number of time points state is recorded
 , nrpt                  = 400
)

} else {
   
params.all  <- data.frame(
## Parameters to establish what model is run
   no_tradeoff           = no_tradeoff
 , parasite_tuning       = parasite_tuning
 , tradeoff_only         = tradeoff_only
 , agg_eff_adjust        = agg_eff_adjust

## If the traits evolving don't start in different places can have a different script for each of these
 ## Completely ignored if model isn't "nt"
 , nt_mut_var_pos_trait  = ifelse(round(qunif(lhs[, 1], min = 1, max = 2)) == 1, T, F)
   
## Unclear about starting with a range for the trait that is evolving
 , pos_trait0            =
   if (model.choice == "nt") {
      qunif(lhs[, 2], min = 0.01, max = 0.50)
 } else if (model.choice == "to") {
      ## ignored in to, transmission is just a function of gamma (neg_trait)
      qunif(lhs[, 2], min = 0.01, max = 0.50)                            
 } else if (model.choice == "eff") {
      ## relative beta in eff
      qunif(lhs[, 2], min = 0.01, max = 0.50)
 }
 , neg_trait0            =
   if (model.choice == "nt") {
      qunif(lhs[, 3], min = 0.01, max = 0.50)
 } else if (model.choice == "to") {
      ## full range allowed here, though may lead to a lot of extinctions...
      qunif(lhs[, 3], min = 0.01, max = 0.99)                            
 } else if (model.choice == "eff") {
      ## full range allowed here, though may lead to a lot of extinctions...
      qunif(lhs[, 3], min = 0.01, max = 0.99)        
 }
   
## Changes by model
 , mut_mean              = 
   if (model.choice == "nt") {
    qunif(lhs[, 4], min = -0.20, max = 0.0)
 } else if (model.choice == "to") {
   ## for now assume no bias in tradeoff only or efficiency model
    0                          
 } else if (model.choice == "eff") {
    0      
 }
   
## Parameters that are the same for all deterministic models
 , mu                    = qunif(lhs[, 5], min = 0.05, max = 0.40)
## mut_sd doesn't mean anything in the deterministic model but leave it around for now
 , mut_sd                = qunif(lhs[, 6], min = 0.01, max = 0.50) 

 , N                     = round(qunif(lhs[, 7], min = 100, max = 2500))
 , gamma0                = qunif(lhs[, 8], min = 0.01, max = 0.4) 
 
## used only for the tradeoff models but include them. Still a bit unclear of the best parameter values here
 ## need a much higher power_c for deterministic formulation however
 , power_c             = qunif(lhs[, 9], min = 0.0005, max = 0.001) * 1000
 , power_exp           = qunif(lhs[, 10], min = 1.5, max = 5.5) 
   
## Extra required parameters that are the same for all models
 , R0_init               = 2
 , deterministic         = deterministic
   
## Not used in deterministic runs
 , nt                    = 2e5
   
## it could potentially be important to change this to be dependent on mu as in the stochastic model
 ## but leave it for now; could be wasting a lot of time though....
 , determ_length         = 400     
 , determ_timestep       = 5       
 , lsoda_hini            = 0.5      
) 

## mu + mut_mean/2 has to be positive for how bias and mutation works
## Not the most elegant, but throwing away parameter values in a random way should be ok...
 ## Code set up to loop over many clusters of parameters until X number are reached, so shouldn't be
  ## much of an issue if many are thrown out
params.all <- params.all %>% filter(mu + mut_mean/2 > 0)
   
}

params.all <- params.all %>% mutate(param_num = seq(1, nrow(params.all)))

