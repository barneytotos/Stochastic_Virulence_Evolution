#############################################################
### Hypercube runs for efficiency model (to in file name) ###
#############################################################

######
## Set up parameter values to explore:
## 1) Time to equilibrium
## 2) Movement away from equilibrium / sd around equilibrium
## 3) Feedback of population size, mutation probability and sd on eco-evo feedbacks and 1 and 2 above
######

######
## Biological assumptions
## 1) Traits evolve either with 0 bias or with a bias towards inf host recovery rate
## 2) Parasites don't take a hit on efficiency when evolving in recovery or take a small hit
  ## another way to say this is 0 correlation in trait evolution or negative. No positive
######

######
## Set up parameter values to range over for tradeoff only model:
## A) start alpha -- starting parasite recovery rate
 # qunif(lhs, min = 0.01, max = 0.99)
## 1) mu          -- mutation rate
 # qunif(lhs, min = 0.005, max = 0.20)
## 2) mut_mean    -- mean in mutation (bias)
 # 0
## 3) mut_sd      -- sd in mutation
 # qunif(lhs, min = 0.05, max = 0.20)
## 4) power_c     -- height of tradeoff curve
 # qunif(lhs, min = 0.005, max = 0.1)
## 5) power_exp   -- slope of tradeoff curve
 # qunif(lhs, min = 1.5, max = 5.5)
## 6) gamma       -- host intrinsic recovery rate
 # qunif(lhs, min = 0.01, max = 0.4)
## 7) N           -- population size
 # qunif(lhs, min = 100, max = 2500)

## Additional parameters to include for efficiency model:
## B) start eff   -- starting parasite efficiency parameter
 # qunif(lhs, min = 0.01, max = 0.99)
## 8) eff_hit     -- size of negative correlation between parasite evo in alpha and efficiency
 # qunif(lhs, min = 0.01, max = 1)
##### 

######
## We will want to explore variation in these parameter values between:
## 1) Tradeoff curve with no second trait. 
  ## tradeoff_only = TRUE; parasite_tuning = FALSE; agg_eff_adjust = FALSE
## 2) Tradeoff curve with proportional beta
  ## tradeoff_only = FALSE; parasite_tuning = FALSE; agg_eff_adjust = FALSE
##### 

## Set up the parameter values over which to sample
num_runs      <- 1
deterministic <- FALSE
num_points    <- 1500

if (!file.exists("lhs_samps_eff.csv")) {
lhs     <- randomLHS(6000, 8) 
write.csv(lhs, "lhs_samps_eff.csv")
} else {
lhs     <- read.csv("lhs_samps_eff.csv")
lhs     <- lhs[, -1]
}

params        <- data.frame(
   parasite_tuning     = FALSE
 , tradeoff_only       = FALSE
 , agg_eff_adjust      = TRUE
## for efficiency runs only add this as a hypercube sample
 , eff_hit             = qunif(lhs[, 8], min = 0.00, max = 1.00) # 0.5
 , num_points          = num_points
 , mut_var             = "beta"
 , mu                  = qunif(lhs[, 1], min = 0.001, max = 0.10) # 0.01 # 
 , mut_mean            = 0
 , mut_sd              = qunif(lhs[, 2], min = 0.01, max = 0.30) # 0.1  # 
## Ignored under conditions of no tuning
 , alpha0              = qunif(lhs[, 3], min = 0.01, max = 0.99) # 0.03 # 
 , tune0               = 0.30 # qunif(lhs, min = 0.01, max = 0.99) # 
 , power_c             = qunif(lhs[, 4], min = 0.005, max = 0.1) # 0.01 # 
 , power_exp           = qunif(lhs[, 5], min = 1.5, max = 5.5) # 3    # 
 , N                   = round(qunif(lhs[, 6], min = 100, max = 2500)) # 600  # 
 , gamma0              = qunif(lhs[, 7], min = 0.01, max = 0.4) # 0.2
 , eff_scale           = 30
 , R0_init             = 2
 , deterministic       = deterministic
)

## duplicate these parameter values for each stochastic run 
num_param <- nrow(params)
params    <- params[rep(seq_len(nrow(params)), each = num_runs), ]
params    <- transform(params, run = rep(seq(1, num_runs), num_param))

## Add in number of time steps based on the mutation rate
params <- transform(params, nt = 0)
for (i in 1:nrow(params)) {
  if (params$mu[i] <= 0.005) {
      params$nt[i] <- 1e6
  } else if (params$mu[i] > 0.005 & params$mu[i] <= 0.01) {
      params$nt[i] <- 5e5
  } else {
      params$nt[i] <- 3e5
  }
}

params <- params %>% mutate(
  rptfreq = max(nt / num_points, 1)
) %>% mutate(
  nrpt    = nt %/% rptfreq
)

params <- transform(params
  , param_num = seq(1, nrow(params))
  , seed      = sample(1:1e5, nrow(params), replace = FALSE)
  , biology   = "efficiency"
  )

## Transform the params to include the optimum alpha
params <- transform(params, opt_alpha = 0)
for (i in 1:nrow(params)) {
params[i, ]$opt_alpha <- with(params[i, ],
  seq(0.01, 0.99, by = 0.01)[which.max(
  power_tradeoff(
  alpha = seq(0.01, 0.99, by = 0.01)
, c     = power_c
, curv  = power_exp) / 
    (1 - (1 - seq(0.01, 0.99, by = 0.01)) * (1 - gamma0))
)]
)
}

for (i in 1:nrow(params)) {
  
  print(i / nrow(params))
  
  time_check <- print(system.time(
    res_1000 <- try(
      with(params
        , run_sim(
   debug4              = F
 , nt                  = nt[i]
 , rptfreq             = rptfreq[i]
 , mut_var             = mut_var[i]
 , seed                = seed[i]
 , mu                  = mu[i]
 , gamma0              = gamma0[i]
 , alpha0              = alpha0[i]
 , tune0               = tune0[i]
 , mut_mean            = mut_mean[i]
 , mut_sd              = mut_sd[i]
 , power_c             = power_c[i]
 , power_exp           = power_exp[i]
 , N                   = N[i]
 , agg_eff_adjust      = agg_eff_adjust[i]
 , tradeoff_only       = tradeoff_only[i]
 , eff_hit             = eff_hit[i]
 , parasite_tuning     = parasite_tuning[i]
 , eff_scale           = eff_scale[i]
 , progress            = "text"
 , R0_init             = R0_init[i]
 , deterministic       = deterministic[i]
## Some defaults here for deterministic run. ALl parameters from here down are ignored if deterministic = FALSE
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

if ((i/50 %% 1) == 0) {
  temp_nam <- paste(paste("res_1000_all_stochas_eff", format(Sys.time(), "%a_%b_%d_%Y"), sep = "_"), ".Rds", sep = "")
  saveRDS(res_1000_all, temp_nam)
}

}


