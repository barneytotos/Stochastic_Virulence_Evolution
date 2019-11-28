#############################################################################################
## Because I need to run such a sparse response surface design, set up parameters manually ##
#############################################################################################

######
## First set are long runs from the top left
######

nt            <- 5e5
num_points    <- 1500
rptfreq       <- max(nt / num_points, 1) 
nrpt          <- nt %/% rptfreq
num_runs      <- 250
deterministic <- FALSE

## there will be 10 total rows for this batch of runs

params <- data.frame(
   nt                  = rep(nt, 9)
 , rptfreq             = rep(rptfreq, 9)
 , nrpt                = rep(nrpt, 9)
 , mut_var             = rep("beta", 9)
 , d                   = rep(0.01, 9)
 , mu                  = c(
     c(0.001, 0.005, 0.025, 0.005, 0.005)  ## 1 from paper: see keynote
   , rep(0.005, 2)                         ## 2 from paper: see keynote
   , rep(0.005, 2))                        ## 3 from paper: see keynote
 , mut_mean            = c(0)
 , mut_sd              = c(
     c(0.05, 0.05, 0.05, 0.01, 0.25)       ## 1 from paper: see keynote
   , rep(0.05, 2)                          ## 2 from paper: see keynote
   , rep(0.05, 2))                         ## 3 from paper: see keynote
 , alpha0              = rep(0.03, 9)
 , tune0               = rep(0.97, 9)
 , tol0                = 1
 , res0                = 1
 , mut_host_mean_shift = 1      
 , mut_host_sd_shift   = 1
 , mut_host_mu_shift   = 100000000000
 , mut_host_res_bias   = 0
 , host_dyn_only       = FALSE
 , power_c             = 0.75
 , power_exp           = 2
 , b_decay             = 2.3
 , b                   = 0.5
 , N                   = c(
    rep(600, 5)
  , c(200, 1800)
  , rep(600, 2))
 , birth_type          = "fill"
 , agg_eff_adjust      = TRUE
 , parasite_tuning     = TRUE
 , eff_scale           = c(
    rep(30, 5)
  , rep(30, 2)
  , c(10, 50))
 , R0_init             = 2
 , deterministic       = deterministic)

######
## Second and third sets are short runs from the bottom left and right
######

nt            <- 2e5
num_points    <- 600
rptfreq       <- max(nt / num_points, 1) 
nrpt          <- nt %/% rptfreq
num_runs      <- 250
deterministic <- FALSE

## there will be 10 total rows for this batch of runs

params2 <- data.frame(
   nt                  = rep(nt, 9)
 , rptfreq             = rep(rptfreq, 9)
 , nrpt                = rep(nrpt, 9)
 , mut_var             = rep("beta", 9)
 , d                   = rep(0.01, 9)
 , mu                  = c(
     c(0.001, 0.005, 0.025, 0.005, 0.005)   ## 1 from paper: see keynote
   , rep(0.005, 2)                          ## 2 from paper: see keynote
   , rep(0.005, 2))                         ## 3 from paper: see keynote
 , mut_mean            = c(0)
 , mut_sd              = c(
     c(0.05, 0.05, 0.05, 0.01, 0.25)        ## 1 from paper: see keynote
   , rep(0.05, 2)                           ## 2 from paper: see keynote
   , rep(0.05, 2))                          ## 3 from paper: see keynote
 , alpha0              = rep(0.97, 9)
 , tune0               = rep(0.03, 9)
 , tol0                = 1
 , res0                = 1
 , mut_host_mean_shift = 1      
 , mut_host_sd_shift   = 1
 , mut_host_mu_shift   = 100000000000
 , mut_host_res_bias   = 0
 , host_dyn_only       = FALSE
 , power_c             = 0.75
 , power_exp           = 2
 , b_decay             = 2.3
 , b                   = 0.5
 , N                   = c(
    rep(600, 5)
  , c(200, 1800)
  , rep(600, 2))
 , birth_type          = "fill"
 , agg_eff_adjust      = TRUE
 , parasite_tuning     = TRUE
 , eff_scale           = c(
    rep(30, 5)
  , rep(30, 2)
  , c(10, 50))
 , R0_init             = 2
 , deterministic       = deterministic)

params3 <- data.frame(
   nt                  = rep(nt, 4)
 , rptfreq             = rep(rptfreq, 4)
 , nrpt                = rep(nrpt, 4)
 , mut_var             = rep("beta", 4)
 , d                   = rep(0.01, 4)
 , mu                  = c(
#    c(0.001, 0.005, 0.025, 0.005, 0.005)   ## 1 from paper: see keynote
     rep(0.005, 2)                          ## 2 from paper: see keynote
   , rep(0.005, 2))                         ## 3 from paper: see keynote
 , mut_mean            = c(0)
 , mut_sd              = c(
#    c(0.05, 0.05, 0.05, 0.01, 0.25)        ## 1 from paper: see keynote
     rep(0.05, 2)                           ## 2 from paper: see keynote
   , rep(0.05, 2))                          ## 3 from paper: see keynote
 , alpha0              = rep(0.03, 4)
 , tune0               = rep(0.03, 4)
 , tol0                = 1
 , res0                = 1
 , mut_host_mean_shift = 1      
 , mut_host_sd_shift   = 1
 , mut_host_mu_shift   = 100000000000
 , mut_host_res_bias   = 0
 , host_dyn_only       = FALSE
 , power_c             = 0.75
 , power_exp           = 2
 , b_decay             = 2.3
 , b                   = 0.5
 , N                   = c(
 #  rep(600, 5)
    c(200, 1800)
  , rep(600, 2))
 , birth_type          = "fill"
 , agg_eff_adjust      = TRUE
 , parasite_tuning     = TRUE
 , eff_scale           = c(
#   rep(30, 5)
    rep(30, 2)
  , c(10, 50))
 , R0_init             = 2
 , deterministic       = deterministic)

######
## Rbind the three data frames
######

params <- rbind(params, params2, params3)
params <- params[rep(seq_len(nrow(params)), each = num_runs), ]

######
## Add runs to the data frame and random seeds
######

params <- transform(params, run = rep(seq(1, num_runs), 22))
params <- transform(params
  , param_num = seq(1, nrow(params))
  , seed      = sample(1:1e5, nrow(params), replace = FALSE))

## Update Nov 1: For trial run, just keep one parameter value set to make sure things are working
# params <- params[params$nt == 2e5, ][1:250, ]

######
## Run the sim
######

i = 1
for (i in 1:nrow(params)) {
  
  print(i / nrow(params))
  
  time_check <- print(system.time(
    res_1000 <- try(
      with(params
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
 , birth_type          = birth_type[i]
 , agg_eff_adjust      = agg_eff_adjust[i]
 , parasite_tuning     = parasite_tuning[i]
 , tradeoff_only       = TRUE               ### !!! new parameter to be added
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
          )))
      , silent = TRUE
      )))
  
  if (class(res_1000) != "try-error") {
 
  ## clean up run i (add parameters and remove NA if the system went extinct)   
res_1000 <- res_1000 %>% mutate(param_num = i, elapsed_time = time_check[3])  
res_1000 <- left_join(res_1000, params, by = "param_num")
res_1000 <- res_1000[complete.cases(res_1000), ]

  if (i == 1) {
res_1000_all <- res_1000
  } else {
res_1000_all <- rbind(res_1000_all, res_1000)
  }

}

if ((i/100 %% 1) == 0) {
  saveRDS(res_1000_all, "res_1000_all_stochas.Rds")
}

}


