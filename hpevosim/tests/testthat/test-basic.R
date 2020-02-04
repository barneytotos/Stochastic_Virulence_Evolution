library(testthat)
## package will be loaded at this point

## FIXME: Ran out of time to work on this for now. I know these should be elsewhere but 
 ## dont have time now to figure out how the package is loading packages or w/e
library(deSolve)
library(ReacTran)

sum_vars <- c("num_S","num_I","num_I_strains","mean_negtrait","sd_negtrait",
              "mean_postrait","sd_postrait","shannon")

run_sim_params <- stoch_params0[intersect(names(stoch_params0),names(formals(run_sim)))]
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

plot(res)  ## goes to look for plot.evosim function

## re-running without tradeoff
run_sim_params_nt <- transform(run_sim_params,
                                    no_tradeoff = TRUE,
                                    mut_mean = -0.05)
res_nt <- do.call(run_sim,run_sim_params_nt)

plot(res_nt)

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

## running again, this time with a deterministic model. 
run_sim_params_determ <- transform(run_sim_params_nt,
                                    deterministic = TRUE,
                                    mut_mean = -0.02,
                                    mu = 0.1,
                                    gamma0 = 0.02)

res_det <- do.call(run_sim,run_sim_params_determ)

## FIXME: Need to improve the way these results are returned and the clarity of these checks.
m_det <- colSums(res_det)
names(m_det) <- NULL
expect_equal(which(m_det != 0), c(1,2,seq(5403, 5502, by = 1)))

m_det <- m_det[which(m_det != 0)]

expect_equal(m_det[1:5], c(4100.0000, 9.8869994595, 0.5694942941, 0.6238255844, 0.7240922323))

res_det_fp <- matrix(data = unlist(res_det[40, -c(1, 2)]), nrow = 100, ncol = 100, byrow = T)
plot(rev(seq(0.00, 0.99, by = 0.01)), res_det_fp[55, ], xlab = "Beta", ylab = "Proportion") 

