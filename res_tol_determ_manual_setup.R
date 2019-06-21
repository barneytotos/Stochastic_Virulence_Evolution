####################################################################################
## Setup deterministic runs with similar parameter values to the stochastic model ##
####################################################################################

nt            <- 5e5
num_points    <- 1500
rptfreq       <- max(nt / num_points, 1) 
nrpt          <- nt %/% rptfreq
num_runs      <- 250
deterministic <- TRUE

#######
## Deterministic parameters from three different starting places, different 
 ## data frame for each parameter set
#######

params <- data.frame(
   nt                  = nt
 , rptfreq             = rptfreq
 , nrpt                = nrpt
 , mut_var             = "beta"
 , d                   = 0.01
 , mu                  = rep(c(0.025, 0.0625, 0.15625, 0.390625), 3)
 , mut_mean            = 0
 , mut_sd              = 1
 , alpha0              = 1 ## not used for deterministic run
 , tune0               = 1 ## not used for deterministic run
 , tol0                = 1
 , res0                = 1
 , mut_host_mean_shift = 1      
 , mut_host_sd_shift   = 1
 , mut_host_mu_shift   = 100000000000
 , mut_host_res_bias   = 0
 , host_dyn_only       = FALSE
 , power_c             = c(rep(20, 4), rep(8, 8))
 , power_exp           = 2
 , b_decay             = 2.3
 , b                   = 0.5
 , N                   = 1
 , balance_birth       = FALSE
 , stochastic_birth    = TRUE
 , fill_birth          = TRUE
 , agg_eff_adjust      = TRUE
 , parasite_tuning     = TRUE
 , eff_scale           = rep(c(10, 30, 50), each = 4)
 , R0_init             = 2
 , determ_length       = 1000
 , determ_timestep     = 2
 , lsoda_hini          = c(rep(0.02, 4), rep(0.05, 8))
 , Imat_seed1          = c(rep(72, 4), rep(80, 8)) ## tune
 , Imat_seed2          = c(rep(28, 4), rep(20, 8)) ## alpha
 , deterministic       = deterministic)

params2 <- data.frame(
   nt                  = nt
 , rptfreq             = rptfreq
 , nrpt                = nrpt
 , mut_var             = "beta"
 , d                   = 0.01
 , mu                  = rep(c(0.025, 0.0625, 0.15625, 0.390625), 3)
 , mut_mean            = 0
 , mut_sd              = 1
 , alpha0              = 1 ## not used for deterministic run
 , tune0               = 1 ## not used for deterministic run
 , tol0                = 1
 , res0                = 1
 , mut_host_mean_shift = 1      
 , mut_host_sd_shift   = 1
 , mut_host_mu_shift   = 100000000000
 , mut_host_res_bias   = 0
 , host_dyn_only       = FALSE
 , power_c             = c(rep(20, 4), rep(8, 8))
 , power_exp           = 2
 , b_decay             = 2.3
 , b                   = 0.5
 , N                   = 1
 , balance_birth       = FALSE
 , stochastic_birth    = TRUE
 , fill_birth          = TRUE
 , agg_eff_adjust      = TRUE
 , parasite_tuning     = TRUE
 , eff_scale           = rep(c(10, 30, 50), each = 4)
 , R0_init             = 2
 , determ_length       = 1000
 , determ_timestep     = 2
 , lsoda_hini          = c(rep(0.05, 4), rep(0.15, 8))
 , Imat_seed1          = c(rep(25, 4), rep(15, 8)) ## tune
 , Imat_seed2          = c(rep(75, 4), rep(85, 8)) ## alpha
 , deterministic       = deterministic)

params3 <- data.frame(
   nt                  = nt
 , rptfreq             = rptfreq
 , nrpt                = nrpt
 , mut_var             = "beta"
 , d                   = 0.01
 , mu                  = rep(c(0.025, 0.0625, 0.15625, 0.390625), 3)
 , mut_mean            = 0
 , mut_sd              = 1
 , alpha0              = 1 ## not used for deterministic run
 , tune0               = 1 ## not used for deterministic run
 , tol0                = 1
 , res0                = 1
 , mut_host_mean_shift = 1      
 , mut_host_sd_shift   = 1
 , mut_host_mu_shift   = 100000000000
 , mut_host_res_bias   = 0
 , host_dyn_only       = FALSE
 , power_c             = c(rep(20, 4), rep(8, 8))
 , power_exp           = 2
 , b_decay             = 2.3
 , b                   = 0.5
 , N                   = 1
 , balance_birth       = FALSE
 , stochastic_birth    = TRUE
 , fill_birth          = TRUE
 , agg_eff_adjust      = TRUE
 , parasite_tuning     = TRUE
 , eff_scale           = rep(c(10, 30, 50), each = 4)
 , R0_init             = 2
 , determ_length       = 1000
 , determ_timestep     = 4
 , lsoda_hini          = 0.20
 , Imat_seed1          = 5 ## tune
 , Imat_seed2          = 5 ## alpha
 , deterministic       = deterministic)

######
## Rbind the three data frames
######

params <- rbind(params, params2, params3)

######
## Add runs to the data frame and random seeds
######

params <- transform(params
  , param_num = seq(1, nrow(params))
  , seed      = sample(1:1e5, nrow(params), replace = FALSE))

#######
## Run the model and clean between each step
#######

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
 , determ_length       = determ_length[i]
 , determ_timestep     = determ_timestep[i]
 , lsoda_hini          = lsoda_hini[i]
          ## tune, alpha
 , Imat_seed           = c(Imat_seed1[i], Imat_seed2[i]))
        ), silent = TRUE
      )))
  
  if (class(res_1000) != "try-error") {
 
## Melt the deteministic matrix into a data frame and simplify
 ## plot on the R0 surface or not?
determ_med_R0 <- data.frame(
  Time     = seq(1, nrow(res_1000))
, med_R0   = 0
, dist_R0  = 0)
for (j in seq(1, nrow(res_1000), by = 1)) {
  
  mat_out <- matrix(data = c(unlist(res_1000[j, -c(1,2)]))
    , nrow = length(tuning)
    , ncol = length(alpha0)
    , byrow = FALSE
 #  , byrow = TRUE
    )
  determ_med_R0[j, ]$dist_R0 <- sum(mat_out * R0_surface)
  determ_med_R0[j, ]$med_R0  <- R0_surface[which(mat_out == max(mat_out), arr.ind = TRUE)]
  
  mat_out <- melt(mat_out) 
  mat_out <- transform(
    mat_out
  , agg  = rep(alpha0, length(tuning))
  , tune = rep(tuning, each = length(alpha0))
  , Time = res_1000$time[j])
  
  if (i == j) {
    mat_out_all <- mat_out
  } else {
    mat_out_all <- rbind(mat_out_all, mat_out)
  }

}

run_name <- paste("res_1000_all_determ", i, sep = "_")
run_name <- paste(run_name , ".Rds", sep = "")
## Ya, ugly. Update this....
run_name <- paste("/Users/Morgan/Documents/Research/McMaster/Bolker_Projects_in_Progress/nearlyneutralHP/MPK_vir_evo/res_out_determ/", run_name, sep = "")

saveRDS(list(mat_out_all, determ_med_R0, params[i, ]), run_name)
 
}

}

#########
### check some runs
#########

determ_med_R0 <- data.frame(
  Time     = seq(1, nrow(res_1000))
, med_R0   = 0
, dist_R0  = 0)
for (j in seq(1, nrow(res_1000), by = 1)) {
  
  mat_out <- matrix(data = c(unlist(res_1000[j, -c(1, 2)]))
    , nrow = length(tuning)
    , ncol = length(alpha0)
    , byrow = FALSE
 #   , byrow = TRUE
    )
  determ_med_R0[j, ]$dist_R0 <- sum(mat_out * R0_surface)
  determ_med_R0[j, ]$med_R0  <- R0_surface[which(mat_out == max(mat_out), arr.ind = TRUE)]
  
  mat_out <- melt(mat_out) 
  mat_out <- transform(
    mat_out
  , agg  = rep(alpha0, length(tuning))
  , tune = rep(tuning, each = length(alpha0))
  , Time = res_1000$time[j])
  
  if (j == 1) {
    mat_out_all <- mat_out
  } else {
    mat_out_all <- rbind(mat_out_all, mat_out)
  }

}

mat_out_all_temp <- droplevels(subset(mat_out_all
  ,   Time == 10  | Time == 20 | Time == 28 | Time == 32 | Time == 46 |
      Time == 60  | Time == 70 | Time == 74 | Time == 100 | Time == 150 |
      Time == 200 
    #  Time == 260 | Time == 300 | Time == 800
  ))
  
  names(mat_out_all_temp)[3] <- "R0"

ggplot(R0_surface_p, aes(tune, agg, z = R0/max(R0))) +
  scale_fill_distiller(palette = "Spectral", na.value = "white") +
  geom_raster(aes(fill = R0/max(R0))) +
  stat_contour(bins = 10) +
  xlab("Parasite Aggressiveness") + 
  ylab("Parasite Tuning") + 
  guides(fill = guide_colorbar(title = "Relative R0")) +
#  scale_x_continuous(
#      breaks = seq(1, length(alpha0), by = 10)
#    , labels = alpha0[seq(1, length(alpha0), by = 10)]) +
#  scale_y_continuous(
#      breaks = seq(1, length(tuning), by = 10)
#    , labels = tuning[seq(1, length(tuning), by = 10)]) +
  theme(legend.key.size = unit(.45, "cm")) +
  geom_vline(xintercept = mut_link_p$linkfun(.205), lwd = 0.5, linetype = "dotted") +
    stat_contour(
     data = mat_out_all_temp
    , aes(group = Time)
    , bins = 5
    , color = "white"
    , alpha = 1.00)
