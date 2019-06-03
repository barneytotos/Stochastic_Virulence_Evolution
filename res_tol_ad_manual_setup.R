#########################################################################
## Setup AD runs with similar parameter values to the stochastic model ##
#########################################################################

nt            <- 5e5
num_points    <- 1500
rptfreq       <- max(nt / num_points, 1) 
nrpt          <- nt %/% rptfreq
num_runs      <- 250
deterministic <- FALSE

mut_link_h <- make.link("log")
mut_link_p <- make.link("logit")
#alpha0    <- c(seq(0.01, 0.99, by = 0.01), 0.999)
#tuning    <- c(seq(0.01, 0.99, by = 0.01), 0.999)
alpha0     <- c(seq(
        mut_link_p$linkfun(0.01), mut_link_p$linkfun(0.99)
        , length = 100))
tuning     <- c(seq(
        mut_link_p$linkfun(0.01), mut_link_p$linkfun(0.99)
        , length = 100))

## For the AD model all we really care about is the shpae of the surface, parameters such as 
 ## mu, mut_sd, N don't mean anything

params <- data.frame(
   nt                  = nt
 , rptfreq             = rptfreq
 , nrpt                = nrpt
 , mut_var             = "beta"
 , d                   = 0.01
 , mu                  = 1
 , mut_mean            = 1
 , mut_sd              = 1
 , tol0                = 1
 , res0                = 1
 , mut_host_mean_shift = 1      
 , mut_host_sd_shift   = 1
 , mut_host_mu_shift   = 100000000000
 , mut_host_res_bias   = 0
 , host_dyn_only       = FALSE
 , power_c             = 2
 , power_exp           = 2
 , b_decay             = 2.3
 , b                   = 0.5
 , N                   = 1
 , balance_birth       = FALSE
 , stochastic_birth    = TRUE
 , fill_birth          = TRUE
 , agg_eff_adjust      = TRUE
 , parasite_tuning     = TRUE
 , eff_scale           = 3 # rep(c(10, 30, 50), 3)
 , R0_init             = 2
 , determ_length       = 1000
 , determ_timestep     = 2
 , lsoda_hini          = 1
 , Imat_seed1          = 13 # rep(c(13, 87, 13), each = 3)
 , Imat_seed2          = 87 # rep(c(87, 13, 13), each = 3)
 , numbins             = 1000
 , deterministic       = deterministic)

######
## Run the sims
######

for (i in 1:nrow(params)) {

  print(i / nrow(params))
  
## Gradient Ascent
## Need power c and power exp?
grad_ascent <- with(params
, par_evo_AD(
  c          = power_c[i]
, curv       = power_exp[i]
, eff_scale  = eff_scale[i]
, mut_link   = mut_link_p
, numbins    = numbins[i]
   ## same parameter as in the RD model, but used a bit differently
, Iseed      = c(Imat_seed1[i], Imat_seed2[i]) * (numbins[i] / length(alpha0))
, simul_mut  = TRUE
, max_range  = FALSE
, debug1     = FALSE
, debug1_val = 759
  ))

grad_ascent <- transform(
  grad_ascent
  , param_num  = i
  , eff_scale  = params$eff_scale[i]
  , numbins    = params$numbins[i]
  , Imat_seed1 = params$Imat_seed1[i]
  , Imat_seed2 = params$Imat_seed2[i])

## Maximum range in which a postive R0 can be obtained
grad_ascent_sing <- with(params
, par_evo_AD(
  c          = power_c[i]
, curv       = power_exp[i]
, eff_scale  = eff_scale[i]
, mut_link   = mut_link_p
, numbins    = numbins[i]
   ## same parameter as in the RD model, but used a bit differently
, Iseed      = c(Imat_seed1[i], Imat_seed2[i]) * (numbins[i] / length(alpha0))
, simul_mut  = FALSE
, max_range  = FALSE
, debug1     = FALSE
, debug1_val = 360
  ))

grad_ascent_sing <- transform(
  grad_ascent_sing
  , param_num  = i
  , eff_scale  = params$eff_scale[i]
  , numbins    = params$numbins[i] 
  , Imat_seed1 = params$Imat_seed1[i]
  , Imat_seed2 = params$Imat_seed2[i])

if (i == 1) {
  grad_ascent_tot      <- grad_ascent
  grad_ascent_sing_tot <- grad_ascent_sing
} else {
  grad_ascent_tot      <- rbind(grad_ascent_tot, grad_ascent)
  grad_ascent_sing_tot <- rbind(grad_ascent_sing_tot, grad_ascent_sing)
    
}

}

saveRDS(grad_ascent_tot, "res_out/res_out_AD/AD_grad_ascent.Rds")
saveRDS(grad_ascent_sing_tot, "res_out/res_out_AD/AD_grad_ascent_sing.Rds")


######
## Stochastic AD version for thesis defense
######

i = 1 ## just for single parameter set
for (j in 1:10) {
  
  print(j / nrow(params))
  
## Gradient Ascent
## Need power c and power exp?
grad_ascent3 <- with(params
, par_evo_AD_rand(
  c          = power_c[i]
, curv       = power_exp[i]
, eff_scale  = eff_scale[i]
, mut_link   = mut_link_p
, numbins    = numbins[i]
   ## same parameter as in the RD model, but used a bit differently
, Iseed      = c(Imat_seed1[i], Imat_seed2[i]) * (numbins[i] / length(alpha0))
, simul_mut  = TRUE
, max_range  = FALSE
, debug1     = FALSE
, debug1_val = 759
  ))

grad_ascent3 <- transform(
  grad_ascent3
  , param_num  = i
  , run_num    = j
  , eff_scale  = params$eff_scale[i]
  , numbins    = params$numbins[i]
  , Imat_seed1 = params$Imat_seed1[i]
  , Imat_seed2 = params$Imat_seed2[i])

if (j == 1) {
  grad_ascent_tot3      <- grad_ascent3
} else {
  grad_ascent_tot3      <- rbind(grad_ascent_tot3, grad_ascent3)
    
}
  
}

saveRDS(grad_ascent_tot, "res_out/res_out_AD/AD_grad_ascent_stochas.Rds")
