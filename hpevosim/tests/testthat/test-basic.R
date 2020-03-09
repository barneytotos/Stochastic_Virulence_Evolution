## library(hpevosim)
library(testthat)

sum_vars <- c("num_S","num_I","num_I_strains","mean_negtrait","sd_negtrait",
              "mean_postrait","sd_postrait","shannon")

run_sim_params <- stoch_params0[intersect(names(stoch_params0),names(formals(run_sim)))]

## reduce nt for faster tests
run_sim_params <- transform(run_sim_params, nt = 5e3)

context("Testing--Stochastic: Tradeoff Only")
res <- do.call(run_sim,run_sim_params)

expect_equal(nrow(res),25)
expect_equal(ncol(res),25)
expect_true(all(res$pop_size==400))

m <- colMeans(tail(res[sum_vars],10))
## dput(round(m,8))
expect_equal(m,
             c(num_S = 87.3, num_I = 312.7, num_I_strains = 17.8,
               mean_negtrait = 0.10641000, sd_negtrait = 0.01894388,  ## gamma
               mean_postrait = 0.04694959, sd_postrait = 0.00211817,    ## beta
               shannon = 1.96662822),
               tolerance=1e-5
)

## goes to look for plot.evosim function
plot(res)  

## re-running without tradeoff
context("Testing--Stochastic: No Tradeoff")
run_sim_params_nt <- transform(run_sim_params,
                                    no_tradeoff = TRUE,
                                    mut_mean = -0.05, 
                                    pos_trait0 = 0.01)

res_nt <- do.call(run_sim,run_sim_params_nt)

plot(res_nt)

m_nt <- colMeans(tail(res_nt[sum_vars],10))

## dput(round(m_nt,8))
expect_equal(m_nt,
             c(num_S = 68.5, num_I = 331.5, num_I_strains = 13.6,
               ## fixed at orig values
               mean_negtrait = 0.01, sd_negtrait = 0,
               ## arbitrarily close to 1
               mean_postrait = 0.94835541, sd_postrait = 0.01208515, 
               shannon = 1.59122209),
               tolerance=1e-5
             )

## running again, this time with a deterministic model.
## evolve only in beta
run_sim_params_determ <- transform(run_sim_params_nt,
                                    deterministic = TRUE,
                                    mut_mean = -0.02,
                                    mu = 0.1,
                                    gamma0 = 0.02,
                                    pos_trait0 = 0.3,
                                    neg_trait0 = 0.05,
                                    determ_length = 100,
                                    lsoda_hini = 1
  )

## Sped up, but now the plots in tests don't show much

context("Testing--Deterministic: No Tradeoff")
res_det  <- do.call(run_sim,run_sim_params_determ)
## Throwing away the parameters, but don't need them for now
## FIXME?: Better integration of parameters 
res_det  <- res_det[["hpevosim_determ.out"]]
## FIXME: progress bar? [Put inside ODE? Doesn't seem to be an option in ode()]

## library(ggplot2)
## ggplot(res_det[["hpevosim_determ.out"]],
##        aes(time,postrait,z=abundance)) + geom_raster()
## abund <- with(res_det[[1]],array(abundance,c(41,100,100)))
## image(abund[,,51])

## Slightly awkward way to do this. Probably best to export? Idea here was that 
 ## the user-specified starting trait would be slotted into the closest bin for the 
  ## deterministic model, but a bit awkward to rely on the output. Just calc from start, but then have to check
## abs(round(...))

init_gamma_inf <- res_det %>% 
  dplyr::filter(time == 0, abundance > 0) %>% 
  dplyr::select(negtrait) %>%
  unlist()

## This round() is annoying
expect_equal(
 (res_det %>% 
  dplyr::filter(time == 40, round(negtrait, 7) == round(init_gamma_inf, 7)) %>% 
  dplyr::summarize(
    total_I   = sum(abundance)
  , mean_beta = weighted.mean(beta, abundance)
  ) %>% unlist())
, c(total_I = 0.7727768933, mean_beta = 0.3043389318),
               tolerance=1e-5
  )

whichtimes <- c(40, 100)

## FIXME: Convert x and y to actual trait values
ggplot(res_det %>% dplyr::filter(time == whichtimes[1])
  , aes(x = postrait_index, y = negtrait_index, z = abundance)) + 
   scale_fill_gradient(low = "white", high = "red4") +
   geom_raster(aes(fill = abundance)) +
  xlab("Transmission Rate (index)") + 
  ylab("Recovery Rate (index)")

(ggplot(res_det
       %>% dplyr::filter(round(negtrait, 7) == round(init_gamma_inf, 7),
       (time == whichtimes[1] | time == whichtimes[2]))
     , aes(beta, abundance, colour = as.factor(time)))
    + geom_line()
    + scale_color_brewer(palette = "Dark2", name = "Time")
    + xlab("Transmission Rate")
    + ylab("Density")
)

ggplot(res_det %>% dplyr::filter(time == whichtime)
  , aes(postrait, abundance)) + 
    geom_line() +
  xlab("Transmission Rate") + 
  ylab("Density")
  
## Plotting sqrt() so that we can see it better
ggplot(res_det %>% dplyr::filter(round(negtrait, 7) == round(init_gamma_inf, 7))
  , aes(x = time, y = postrait_index, z = abundance)) +
   scale_fill_gradient(low = "white", high = "red4") +
    geom_raster(aes(fill = sqrt(abundance))) +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) +
    xlab("Time") + 
    ylab("Transmission Rate (index)")

library(rgl)

## Problems with "increasing x and y values expected" using persp3d. Don't feel like taking the time to figure it out right now
with(res_det %>% filter(round(negtrait, 7) == round(init_gamma_inf, 7)) %>% droplevels()
  , plot3d(
  x = time, y = postrait_index, z = abundance, xlab = "time", ylab = "beta prob",
  zlim = c(0, 0.4), col = "black"
))

## running again, this time with a tradeoff curve.
## evolve only in gamma, beta gets pulled along according to the tradeoff curve
run_sim_params_determ <- transform(run_sim_params_determ,
                                    power_c =    2,
                                    no_tradeoff = FALSE, 
                                    tradeoff_only = TRUE
  )

context("Testing--Deterministic: Tradeoff only")
res_det  <- do.call(run_sim,run_sim_params_determ)
res_det  <- res_det[["hpevosim_determ.out"]]
## Semi-confusing because of differences between models, but with tradeoff only "postrait" isn't evolving, 
 ## so pull out the starting value of postrait, then later can plot negtrait and beta, which is just 
  ## a function of negtrait in this model
init_beta_inf <- res_det %>% 
  filter(time == 0, abundance > 0) %>% 
  dplyr::select(postrait) %>%
  unlist()

expect_equal(
 (res_det %>% 
  filter(time == 40, round(postrait, 7) == round(init_beta_inf, 7)) %>% 
  summarize(
    total_I   = sum(abundance)
  , mean_beta = weighted.mean(beta, abundance)
  ) %>% unlist())
, c(total_I = 0.9070393943, mean_beta = 0.7178273596),
               tolerance=1e-5
  )

whichtimes <- c(40, 100)

## When evolving only in one trait (along the tradeoff curve) the state-space plot is pretty meaningless, so skip it
 ## Transmission rate has been pre-calculated inside of tidy.deSolve inside of run_sim using the tradeoff curve parameters
ggplot(res_det %>% filter(time == whichtimes, round(postrait, 7) == round(init_beta_inf, 7))
  , aes(beta, abundance, colour = as.factor(time))) + 
    geom_line() +
  scale_color_brewer(palette = "Dark2", name = "Time") +
  xlab("Transmission Rate") + 
  ylab("Density") 
  
  ## Once again, plotting sqrt() so that we can see it better
ggplot(res_det %>% filter(round(postrait, 7) == round(init_beta_inf, 7))
  , aes(x = time, y = negtrait_index, z = abundance)) +
   scale_fill_gradient(low = "white", high = "red4") +
   geom_raster(aes(fill = sqrt(abundance))) +
  xlab("Time") + 
  ylab("Recovery Rate (index)")

## FIXME: I have some vision of doing this so wanted to get a start on it, but it isn't correct yet

## Recover the tradeoff curve from the output
res_det.tc <- res_det %>% group_by(negtrait) %>% summarize(beta = mean(beta))
## And summarize quantiles through time
res_det.q  <- res_det %>% group_by(time) %>% 
  summarize(
    negtrait = weighted.mean(negtrait, abundance)
  , beta     = weighted.mean(beta, abundance)
    )

## Plot movement of the distribution on the tradeoff curve
ggplot(res_det.tc %>% filter(negtrait < 0.1)
  , aes(negtrait, beta)) + 
  geom_line() +
  geom_path(data = res_det.q, colour = "red4", lwd = 1) +
  xlab("Recovery Rate") +
  ylab("Transmission Rate") 

## Run for efficiency model
run_sim_params_determ <- transform(run_sim_params_determ,
                                    no_tradeoff   = F, 
                                    tradeoff_only = F,
                                    determ_length = 200,
                                    gamma0        = 0.2, ## increasing "background recovery" so that the negtrait optimum isn't so close to 0
                                    mut_mean      = 0.00
   )

## NOTE: Some things to explore for efficiency model
 ## a few negative mut_mean tried so far seems to pull gamma away from tradeoff optimum as efficiency increases which
  ## is leading to higher beta. Probably if this was run long enough it would escape this (I believe in the thesis this
   ## was a partial cause of the roundabout path to the optimum), but we figured netural gamma is probably best for this
    ## scenario so putting this down for now

context("Testing--Deterministic: Efficiency")
res_det  <- do.call(run_sim,run_sim_params_determ)
## Throwing away the parameters, but don't need them for now
## FIXME?: Better integration of parameters 
res_det  <- res_det[["hpevosim_determ.out"]]
## FIXME: progress bar? [Put inside ODE? Doesn't seem to be an option in ode()]

## This round() is annoying
expect_equal(
 (res_det %>% 
  filter(time == 40) %>% 
  summarize(
    total_I   = sum(abundance)
  , mean_beta = weighted.mean(beta, abundance)
  ) %>% unlist())
, c(total_I = 0.05965713701, mean_beta = 0.24803734939),
               tolerance=1e-5
  )

whichtimes <- c(40, 100, 200)

## sorta ugly
opt_negtrait <- which.max(power_R0(alpha = unique(res_det$negtrait), c = run_sim_params_determ$power_c, curv = run_sim_params_determ$power_exp, gamma = run_sim_params_determ$gamma0, N = 1))

## FIXME: Convert x and y to actual trait values.
 ## Would also like to plot on same plot.
ggplot(res_det %>% filter(time == whichtimes[1] | time == whichtimes[2] | time == whichtimes[3])
  , aes(x = postrait_index, y = negtrait_index, z = abundance)) + 
   scale_fill_gradient(low = "white", high = "red4") +
   geom_raster(aes(fill = abundance)) +
  facet_wrap(~time) +
  geom_hline(aes(yintercept = opt_negtrait), linetype = "dashed", lwd = 0.4) +
  xlab("Efficiency (index)") + 
  ylab("Recovery Rate (index)")

res_det.mb <- res_det %>% group_by(time) %>%
  summarize(
    mean_gamma = weighted.mean(x = negtrait, w = abundance)
  , mean_beta  = weighted.mean(x = beta, w = abundance)
  , mean_eff   = weighted.mean(x = postrait, w = abundance))

grid.arrange(
  ggplot(res_det.mb, aes(time, mean_gamma)) + geom_line()
, ggplot(res_det.mb, aes(time, mean_beta)) + geom_line()
, ggplot(res_det.mb, aes(time, mean_eff)) + geom_line()
, ncol = 1
)

## FIXME: Need a cleaner way to plot distribution of beta through time
  
## Finally, run for tuning model
#run_sim_params_determ <- transform(run_sim_params_determ,
#                                    parasite_tuning = T
#   )
#
#res_det  <- do.call(run_sim,run_sim_params_determ)
#res_det  <- res_det[["hpevosim_determ.out"]]
#
#expect_equal(
# (res_det %>% 
#  filter(time == 40) %>% 
#  summarize(
#    total_I   = sum(abundance)
#  , mean_beta = weighted.mean(beta, abundance)
#  ) %>% unlist())
#, c(total_I = 0.6160694006, mean_beta = 0.7408605573),
#               tolerance=1e-5
#  )
#
#whichtimes <- c(40, 200, 400)
#
## sorta ugly
#opt_negtrait <- which.max(power_R0(alpha = unique(res_det$negtrait)
#  , c = run_sim_params_determ$power_c
#  , curv = run_sim_params_determ$power_exp, gamma = run_sim_params_determ$gamma0, N = 1))
#
#ggplot(res_det %>% filter(time == whichtimes[1] | time == whichtimes[2] | time == whichtimes[3])
#  , aes(x = postrait_index, y = negtrait_index, z = abundance)) + 
#   scale_fill_gradient(low = "white", high = "red4") +
#   geom_raster(aes(fill = abundance)) +
#  facet_wrap(~time) +
#  geom_hline(aes(yintercept = opt_negtrait), linetype = "dashed", lwd = 0.4) +
#  xlab("Efficiency (index)") + 
#  ylab("Recovery Rate (index)")
#
#res_det.mb <- res_det %>% group_by(time) %>%
#  summarize(
#    mean_gamma  = weighted.mean(x = negtrait, w = abundance)
#  , mean_beta   = weighted.mean(x = beta, w = abundance)
#  , mean_tuning = weighted.mean(x = postrait, w = abundance)
#  , mean_eff    = weighted.mean(x = effic, w = abundance))
#
#grid.arrange(
#  ggplot(res_det.mb, aes(time, mean_gamma)) + geom_line()
#, ggplot(res_det.mb, aes(time, mean_beta)) + geom_line()
#, ggplot(res_det.mb, aes(time, mean_tuning)) + geom_line()
#, ggplot(res_det.mb, aes(time, mean_eff)) + geom_line()  
#, ncol = 2
#)

