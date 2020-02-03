library(testthat)
## run from head
source("R/utils.R")

nt <- 1e4
num_points <- 50
rptfreq       <- max(nt / num_points, 1) 
num_runs      <- 1
deterministic <- FALSE ## run the model via reaction-diffusion deterministic?

source("funs_SIS.R")
params <- list(
               no_tradeoff          = FALSE
             , nt_mut_var_pos_trait = TRUE
             , pos_trait0           = 0.005
             , nt                   = nt
             , rptfreq              = rptfreq
             , deterministic        = deterministic
             , parasite_tuning      = FALSE
             , tradeoff_only        = TRUE
             , agg_eff_adjust       = FALSE
             , mut_var              = "beta"
             , d                    = 0.01
               ## Need to convert this to a rate of diffusion if deterministic == TRUE
             , mu                   = if (!deterministic) 0.01 else 0.002
             , mut_mean             = 0
             , mut_sd               = 0.15
             , neg_trait0           = 0.01
             , tune0                = 0.03
             , tol0                 = 1
             , res0                 = 1
             , gamma0               = 0.2
             , run                  = 1
             , mut_host_mean_shift  = 1      
             , mut_host_sd_shift    = 1
             , mut_host_mu_shift    = 100000000000 # need some huge number here. Not used in this iteration anyway
             , mut_host_res_bias    = 0
             , host_dyn_only        = FALSE
             , power_c              = 0.1
             , power_exp            = 3 # 2
             , b_decay              = 2.3
             , b                    = 0.5
             , N                    = 400
             , birth_type           = "fill"
             , eff_scale            = 30
             , param_num = 1
             , seed = 101
               )
sum_vars <- c("num_S","num_I","num_I_strains","mean_negtrait","sd_negtrait",
              "mean_postrait","sd_postrait","shannon")

run_sim_params <- params[intersect(names(params),names(formals(run_sim)))]
res <- do.call(run_sim,run_sim_params)

expect_equal(nrow(res),50)
expect_equal(ncol(res),18)
expect_true(all(res$pop_size==400))

m <- colMeans(tail(res[sum_vars],10))
## dput(round(m,8))
expect_equal(m,
             c(num_S = 89.8, num_I = 310.2, num_I_strains = 15.9,
               mean_negtrait = 0.1216931, sd_negtrait = 0.01837436,
               mean_postrait = 0.0493278, sd_postrait = 0.0025132, 
               shannon = 1.79367083)
)

plot.evosim(res)

run_sim_params_nt <- transform(run_sim_params,
                                    no_tradeoff = TRUE,
                                    mut_mean = -0.05)
res_nt <- do.call(run_sim,run_sim_params_nt)

plot.evosim(res_nt)

m_nt <- colMeans(tail(res_nt[sum_vars],10))

## dput(round(m_nt,8))
expect_equal(m_nt,
             c(num_S = 69.3, num_I = 330.7, num_I_strains = 17.2,
               ## fixed at orig values
               mean_negtrait = 0.01, sd_negtrait = 0,
               ## arbitrarily close to 1
               mean_postrait = 0.98985252, sd_postrait = 0.00909884, 
               shannon = 1.99599423)
             )

