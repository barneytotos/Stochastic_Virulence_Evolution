##' default parameters
##' @export
stoch_params0 <- list(
               no_tradeoff          = FALSE
             , nt_mut_var_pos_trait = TRUE
             , pos_trait0           = 0.03
             , nt                   = 1e4
             , rptfreq              = 200
             , deterministic        = FALSE
             , parasite_tuning      = FALSE
             , tradeoff_only        = TRUE
             , agg_eff_adjust       = FALSE
             , mut_var              = "beta"
             , d                    = 0.01
               ## Need to convert this to a rate of diffusion if deterministic == TRUE (0.002)
             , mu                   = 0.01
             , mut_mean             = 0
             , mut_sd               = 0.15
             , neg_trait0           = 0.01
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
