#########################################################################
## Run run_sim many times across either parameter values or, first,    ##
## a single set of parameter values to capture magnitude of the effect ##
## of stochasticity as a debug like tool                               ##
#########################################################################

i = 1
for (i in 1:nrow(params)) {
  
  print(i / nrow(params))
  
  time_check <- print(system.time(res_1000 <- with(params
        , run_sim(
   nt                  = nt[i]
 , rptfreq             = rptfreq[i]
 , mut_var             = mut_var[i]
 , seed                = seed[i]
 , d                   = d[i]
 , mu                  = mu[i]
 , alpha0              = alpha0[i]
 , tune0               = tune0[i]
 , tol0                = tol0[i]
 , res0                = res0[i]
 , mut_host_mean_shift = mut_host_mean_shift[i]
 , mut_host_sd_shift   = mut_host_sd_shift[i]
 , mut_host_mu_shift   = mut_host_mu_shift[i]
 , mut_host_res_bias   = mut_host_res_bias[i]
 , host_dyn_only       = host_dyn_only[i]
 , mut_mean            = mut_mean[i]
 , mut_sd              = mut_sd[i]
 , power_c             = power_c[i]
 , power_exp           = power_exp[i]
 , b_decay             = b_decay[i]
 , b                   = b[i]
 , N                   = N[i]
 , balance_birth       = balance_birth[i]
 , stochastic_birth    = stochastic_birth[i]
 , fill_birth          = fill_birth[i]
 , agg_eff_adjust      = agg_eff_adjust[i]
 , parasite_tuning     = parasite_tuning[i]
 , tradeoff_only       = tradeoff_only[i]
 , eff_scale           = eff_scale[i]
 , debug               = FALSE
 , debug2              = FALSE
 , debug3              = FALSE
 , host_evo_delay      = FALSE
 , host_evo_delay_start= 30
 , host_evo_delay_stop = 140
 , progress            = TRUE
 , R0_init             = R0_init[i]
 , deterministic       = deterministic[i]
## Some defaults here for deterministic run. Ignored if deterministic = FALSE
 , determ_length       = 400
 , determ_timestep     = 1
 , lsoda_hini          = 0.50
## Choose these from the parameter values data frame so that the deterministic run starts from the 
 ## same place as the stochastic run
 , Imat_seed           = c(
#   which(c(round(seq(0.01, 0.99, by = 0.01), 2), 0.999) == params$alpha0[i])
# , which(c(round(seq(0.01, 0.99, by = 0.01), 2), 0.999) == params$tune0[i]))
   15        ## tune
 , 85        ## alpha
          )))))
 
  ## clean up run i (add parameters and remove NA if the system went extinct)   
res_1000 <- res_1000 %>% mutate(param_num = i, elapsed_time = time_check[3])  
res_1000 <- left_join(res_1000, params, by = "param_num")
res_1000 <- res_1000[complete.cases(res_1000), ]

  if (i == 1) {
res_1000_all <- res_1000
  } else {
res_1000_all <- rbind(res_1000_all, res_1000)
  }

if (i %% 20 == 0) {
  saveRDS(res_1000_all, "res_1000_all2.Rds")
}

}

res_1000_all <- transform(res_1000_all, mut_mean = as.factor(mut_mean))
#res_1000_all <- transform(res_1000_all, run = as.numeric(run))
res_1000_all <- transform(res_1000_all, tol0 = as.factor(tol0))
res_1000_all <- transform(res_1000_all, res0 = as.factor(res0))

## Little bit of code for debugging/plotting single runs
res_1000_all     <- transform(res_1000_all, R0 = 0)
res_1000_stochas <- res_1000_all

res_1000_stochas <- rbind(res_1000_all, res_1000_stochas)


######
## Summarize multiple runs 
######

res_1000_stochas_s <- res_1000_stochas %>%
  group_by(time) %>%
  summarize(
    q05_mean_plalpha = quantile(mean_plalpha, 0.05)
  , q25_mean_plalpha = quantile(mean_plalpha, 0.25)
  , q50_mean_plalpha = quantile(mean_plalpha, 0.50)
  , q75_mean_plalpha = quantile(mean_plalpha, 0.75)
  , q95_mean_plalpha = quantile(mean_plalpha, 0.95)
  , q05_mean_plbeta  = quantile(mean_plbeta, 0.05)
  , q25_mean_plbeta  = quantile(mean_plbeta, 0.25)
  , q50_mean_plbeta  = quantile(mean_plbeta, 0.50)
  , q75_mean_plbeta  = quantile(mean_plbeta, 0.75)
  , q95_mean_plbeta  = quantile(mean_plbeta, 0.95)) %>%
  mutate(R0 = 0)


############
## AD run ##
############

## Just the deterministic version run here, can also run the stochastic AD with:
 ## par_evo_AD_rand function

## Gradient Ascent
## Need power c and power exp?
grad_ascent <- par_evo_AD(
  c          = 4
, curv       = 2
, eff_scale  = 30
, mut_link   = mut_link_p
, numbins    = 1000
, Iseed      = c(950, 50)  ## tune, agg
, simul_mut  = TRUE
, max_range  = FALSE
, debug1     = FALSE
, debug1_val = 360)

## Maximum range in which a postive R0 can be obtained
grad_max_range <- par_evo_AD(
  c          = 4
, curv       = 2
, eff_scale  = 50
, mut_link   = mut_link_p
, numbins    = 1000 
, Iseed      = c(50, 950)
, simul_mut  = TRUE
, max_range  = TRUE)

