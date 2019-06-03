##################################################
## Run across all needed parameter combinations ##
##################################################

##########
### Deterministic runs
##########

nt            <- 2e5
num_points    <- 600
rptfreq       <- max(nt / num_points, 1) 
nrpt          <- nt %/% rptfreq
num_runs      <- 1
deterministic <- FALSE ## run the model via reaction-diffusion deterministic?

## Set up parameters in a manner to correspond to output data frame
params <- expand.grid(
   nt                  = nt
 , rptfreq             = rptfreq
 , nrpt                = nrpt
 , mut_var             = "beta"
 , d                   = 0.01
 , mu                  = c(0.005, 0.010, 0.050, 0.10)
 , mut_mean            = 0
 , mut_sd              = c(0.01) 
 , alpha0              = c(0.05, 0.95) 
 , tune0               = c(0.05, 0.95)
 , tol0                = 1
 , res0                = 1
 , run                 = seq(1, num_runs, by = 1)
 , mut_host_mean_shift = 1      
 , mut_host_sd_shift   = 1
 , mut_host_mu_shift   = 100000000000
 , mut_host_res_bias   = 0
 , host_dyn_only       = FALSE
 , power_c             = 4
 , power_exp           = 2
 , b_decay             = 2.3
 , b                   = 0.5
 , N                   = 1000
 , balance_birth       = FALSE
 , stochastic_birth    = TRUE
 , fill_birth          = TRUE
 , agg_eff_adjust      = TRUE
 , parasite_tuning     = TRUE
 , eff_scale           = c(10, 30, 50, 70, 90)
 , R0_init             = 2
 , deterministic       = deterministic
  )

params <- transform(params
  , param_num = seq(1, nrow(params))
  , seed      = sample(1:1e5, nrow(params), replace = FALSE))

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
 , eff_scale           = eff_scale[i]
 , debug2              = FALSE
 , debug3              = FALSE
 , host_evo_delay      = FALSE
 , host_evo_delay_start= 30
 , host_evo_delay_stop = 140
 , progress            = TRUE
 , R0_init             = R0_init[i]
 , deterministic       = deterministic[i]
## Some defaults here for deterministic run. Ignored if deterministic = FALSE
 , determ_length       = 600
 , determ_timestep     = 5
 , lsoda_hini          = 0.2
## Choose these from the parameter values data frame so that the deterministic run starts from the 
 ## same place as the stochastic run
 , Imat_seed           = c(
   which(c(round(seq(0.01, 0.99, by = 0.01), 2), 0.999) == params$alpha0[i])
 , which(c(round(seq(0.01, 0.99, by = 0.01), 2), 0.999) == params$tune0[i]))
          ))))
 
  ## clean up run i (add parameters and remove NA if the system went extinct)   
res_1000 <- res_1000 %>% mutate(param_num = i, elapsed_time = time_check[3])  
res_1000 <- left_join(res_1000, params, by = "param_num")
res_1000 <- res_1000[complete.cases(res_1000), ]

  if (i == 1) {
res_1000_all <- res_1000
  } else {
res_1000_all <- rbind(res_1000_all, res_1000)
  }

if ((i/20 %% 1) == 0) {
  saveRDS(res_1000_all, "res_1000_all_determ.Rds")
}

}

saveRDS(res_1000_all, "res_1000_all_determ.Rds")

##########
### Stochastic runs
##########

nt            <- 2e5
num_points    <- 600
rptfreq       <- max(nt / num_points, 1) 
nrpt          <- nt %/% rptfreq
num_runs      <- 1
deterministic <- FALSE ## run the model via reaction-diffusion deterministic?

## Set up parameters in a manner to correspond to output data frame
params <- expand.grid(
   nt                  = nt
 , rptfreq             = rptfreq
 , nrpt                = nrpt
 , mut_var             = "beta"
 , d                   = 0.01
 , mu                  = c(0.001, 0.005, 0.01)
 , mut_mean            = 0
 , mut_sd              = c(0.01, 0.05, 0.20) 
 , alpha0              = c(0.03) # c(0.03, 0.97) 
 , tune0               = c(0.97) # c(0.03, 0.97)
 , tol0                = 1
 , res0                = 1
 , run                 = seq(1, num_runs, by = 1)
 , mut_host_mean_shift = 1      
 , mut_host_sd_shift   = 1
 , mut_host_mu_shift   = 100000000000
 , mut_host_res_bias   = 0
 , host_dyn_only       = FALSE
 , power_c             = 0.75
 , power_exp           = 2
 , b_decay             = 2.3
 , b                   = 0.5
 , N                   = 500 # c(100, 500, 2000)
 , balance_birth       = FALSE
 , stochastic_birth    = TRUE
 , fill_birth          = TRUE
 , agg_eff_adjust      = TRUE
 , parasite_tuning     = TRUE
 , eff_scale           = 30 # c(10, 30, 50)
 , R0_init             = 2
 , deterministic       = deterministic)

params <- params %>% filter(!(alpha0 == 0.97 & tune0 == 0.97))

params <- transform(params
  , param_num = seq(1, nrow(params))
  , seed      = sample(1:1e5, nrow(params), replace = FALSE))

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
 , eff_scale           = eff_scale[i]
 , debug2              = FALSE
 , debug3              = FALSE
 , host_evo_delay      = FALSE
 , host_evo_delay_start= 30
 , host_evo_delay_stop = 140
 , progress            = TRUE
 , R0_init             = R0_init[i]
 , deterministic       = deterministic[i]
## Some defaults here for deterministic run. Ignored if deterministic = FALSE
 , determ_length       = 600
 , determ_timestep     = 5
 , lsoda_hini          = 0.2
## Choose these from the parameter values data frame so that the deterministic run starts from the 
 ## same place as the stochastic run
 , Imat_seed           = c(
   which(c(round(seq(0.01, 0.99, by = 0.01), 2), 0.999) == params$alpha0[i])
 , which(c(round(seq(0.01, 0.99, by = 0.01), 2), 0.999) == params$tune0[i]))
          ))))
 
  ## clean up run i (add parameters and remove NA if the system went extinct)   
res_1000 <- res_1000 %>% mutate(param_num = i, elapsed_time = time_check[3])  
res_1000 <- left_join(res_1000, params, by = "param_num")
res_1000 <- res_1000[complete.cases(res_1000), ]

  if (i == 1) {
res_1000_all <- res_1000
  } else {
res_1000_all <- rbind(res_1000_all, res_1000)
  }

if ((i/20 %% 1) == 0) {
  saveRDS(res_1000_all, "res_1000_all_stochas.Rds")
}

}

saveRDS(res_1000_all, "res_1000_all_stochas.Rds")
