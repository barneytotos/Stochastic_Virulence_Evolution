num_runs      <- 1
deterministic <- FALSE
num_points    <- 1000

## Generate way more than could possibly be needed
if (!file.exists("lhs_samps_eff.csv")) {
lhs     <- randomLHS(100000, 10) 
write.csv(lhs, "lhs_samps_eff.csv")
} else {
lhs     <- read.csv("lhs_samps_eff.csv")
lhs     <- lhs[, -1]
}

params.all  <- data.frame(
  no_tradeoff            = TRUE
 , nt_mut_var_pos_trait  = TRUE
 , pos_trait0            = 0.005
 , parasite_tuning       = TRUE
 , tradeoff_only         = FALSE
 , agg_eff_adjust        = FALSE
 , num_points            = num_points
 , nt                    = 2e5
 , mu                    = qunif(lhs[, 1], min = 0.001, max = 0.10)
 , mut_mean              = qunif(lhs[, 2], min = -1.0, max = 0.0)
 , mut_sd                = qunif(lhs[, 3], min = 0.01, max = 0.50) 
 , N                     = round(
                           qunif(lhs[, 6], min = 100, max = 2500))
 , gamma0                = qunif(lhs[, 7], min = 0.01, max = 0.4) 
 , R0_init               = 2
 , deterministic         = deterministic
)

num_complete <- 0
j            <- 1

while (num_complete < 500 | j < 10) {
  
## record the batch 
batch.rows <- (1000 * (j - 1) + 1):(1000 * (j - 1) + 1000)
params     <- params.all[batch.rows, ]  

## duplicate these parameter values for each stochastic run 
num_param <- nrow(params)
params    <- params[rep(seq_len(nrow(params)), each = num_runs), ]
params    <- transform(params, run = rep(seq(1, num_runs), num_param))

params <- params %>% mutate(
  rptfreq = max(nt / num_points, 1)
) %>% mutate(
  nrpt    = nt %/% rptfreq
)

params <- transform(params
  , param_num = seq(1, nrow(params))
  , seed      = sample(1:1e5, nrow(params), replace = FALSE)
  , biology   = "no tradeoff"
  )

for (i in 1:nrow(params)) {
  
  print(i / nrow(params))
  
  time_check <- print(system.time(
    res_1000 <- try(
      with(params
        , run_sim(
   no_tradeoff          = no_tradeoff[i]
 , nt_mut_var_pos_trait = nt_mut_var_pos_trait[i]
 , nt                   = nt[i]
 , rptfreq              = rptfreq[i]
 , seed                 = seed[i]
 , mu                   = mu[i]
 , gamma0               = gamma0[i]
 , pos_trait0           = pos_trait0[i]
 , mut_mean             = mut_mean[i]
 , mut_sd               = mut_sd[i]
 , power_c              = power_c[i]
 , power_exp            = power_exp[i]
 , N                    = N[i]
 , progress             = "none"
 , R0_init              = R0_init[i]
 , deterministic        = deterministic[i]
          ))
      , silent = TRUE
      )
    ))
  
  if (class(res_1000) != "try-error") {
    
  ## first check 
 
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

  ## Save completed runs with all of the chopped runs in batches of 100 completed runs
if (((num_complete/100) %% 1) == 0) {
  Sys.sleep(1) ## For whatever reason struggles without this, no idea why
  temp_nam <- paste(paste("batch_runs/res_1000_all_stochas_nt", paste(num_complete, format(Sys.time(), "%a_%b_%d_%Y"), sep = "_"), sep = "_"), ".Rds", sep = "")
  saveRDS(res_1000_all, temp_nam)
  Sys.sleep(1)
}
  
print(i/num_complete)
  
}

## Add to the batch counter
j <- j + 1

}


# res <- res_1000_all[res_1000_all$param_num == 4, ]
# class(res) <- c("hpevosim","data.frame")  
