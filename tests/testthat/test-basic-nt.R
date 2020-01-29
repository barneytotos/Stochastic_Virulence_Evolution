library(testthat)
## stub for testing
nt <- 1e4
num_points <- 50
rptfreq       <- max(nt / num_points, 1) 
num_runs      <- 1
deterministic <- FALSE ## run the model via reaction-diffusion deterministic?

source("funs_SIS.R")
params <- list(
               no_tradeoff          = TRUE
             , nt_mut_var_pos_trait = TRUE
             , pos_trait0           = 0.005
             , nt                   = nt
             , rptfreq              = rptfreq
             , deterministic        = FALSE
             , parasite_tuning      = FALSE
             , tradeoff_only        = TRUE
             , agg_eff_adjust       = FALSE
             , mut_var              = "beta"
             , d                    = 0.01
               ## Need to convert this to a rate of diffusion if deterministic == TRUE
             , mu                   = if (!deterministic) 0.01 else 0.002
             , mut_mean             = -0.05
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

# debug(run_sim); debug(do_mut)
run_sim_params <- params[intersect(names(params),names(formals(run_sim)))]
res <- do.call(run_sim,run_sim_params)


##
expect_equal(nrow(res),50)
expect_equal(ncol(res),18)          
expect_true(all(res$pop_size==400))

## MPK: naming has changed
sum_vars <- c("num_S","num_I","num_I_strains","mean_negtrait","sd_negtrait",
              "mean_postrait","sd_postrait","shannon")
m <- colMeans(tail(res[sum_vars],10))

## MPK: names changed values have not
expect_equal(m,
             c(
                     num_S = 69.3
                   , num_I = 330.7
                   , num_I_strains = 17.2
                   , mean_negtrait = 0.01
                   , sd_negtrait = 0.0
                   , mean_postrait = 0.9898525183219263
                   , sd_postrait = 0.009098842919324621
                   , shannon = 1.995994234250763))

## plot (not part of a real pipeline)
alpha_vars <- res[,c("mean_postrait")]
alpha_vars$mean_postrait <- log(alpha_vars$mean_postrait)
matplot(res$time,
       alpha_vars,
       type="l",
       xlab="time",
       ylab="alpha",
       lwd=c(2,1,1,1))
