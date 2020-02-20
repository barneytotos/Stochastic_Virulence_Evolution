## library(hpevosim)
library(testthat)

## FIXME: add context() statements
## FIXME: shorten runs to speed up tests?

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
## evolve only in beta
## imatseed
run_sim_params_determ <- transform(run_sim_params_nt,
                                    deterministic = TRUE,
                                    mut_mean = -0.02,
                                    mu = 0.1,
                                    gamma0 = 0.02,
                                    pos_trait0 = 0.3,
                                    neg_trait0 = 0.05
#  , no_tradeoff   = F
#  , tradeoff_only = F
  )

res_det  <- do.call(run_sim,run_sim_params_determ)
## FIXME: progress bar? [Put inside ODE? Doesn't seem to be an option in ode()]

## Slightly awkward way to do this. Probably best to export? Idea here was that 
 ## the user-specified starting trait would be slotted into the closest bin for the 
  ## deterministic model, but a bit awkward to rely on the output. Just calc from start, but then have to check
   ## abs(round(...))
init_gamma_inf <- res_det %>% 
  filter(time == 0, abundance > 0) %>% 
  dplyr::select(negtrait) %>%
  unlist()

## This round() is annoying
expect_equal(
 (res_det %>% 
  filter(time == 40, round(negtrait, 7) == round(init_gamma_inf, 7)) %>% 
  summarize(
    total_I   = sum(abundance)
  , mean_beta = sum(abundance * beta)
  ) %>% unlist())
, c(total_I = 0.7727769, mean_beta = 0.2351861)
  )

whichtime <- 80

## FIXME: Convert x and y to actual trait values
ggplot(res_det %>% filter(time == whichtime)
  , aes(x = postrait_index, y = negtrait_index, z = abundance)) + 
   scale_fill_gradient(low = "white", high = "red4") +
   geom_raster(aes(fill = abundance)) +
  xlab("Transmission Rate") + 
  ylab("Recovery Rate")

ggplot(res_det %>% filter(time == whichtime)
  , aes(postrait, abundance)) + 
    geom_line() +
  xlab("Transmission Rate") + 
  ylab("Density")
  
## Plotting sqrt() so that we can see it better
ggplot(res_det %>% filter(round(negtrait, 7) == round(init_gamma_inf, 7))
  , aes(x = time, y = postrait_index, z = abundance)) +
   scale_fill_gradient(low = "white", high = "red4") +
   geom_raster(aes(fill = sqrt(abundance))) +
  xlab("Time") + 
  ylab("Recovery Rate")

#### Some old stuff >> to incorporate

## normalize; if you normalize by column (i.e. distribution across time)
## it looks weird and interesting but maybe meaningless)
bmat_flip_norm <- sweep(bmat_flip,MARGIN=1,rowSums(bmat_flip),"/")
image(y=as.numeric(colnames(bmat_flip)),
      z=bmat_flip_norm, xlab="time",ylab="beta prob",
      zlim=c(0,0.4) ## restrict scale so not overwhelmed by initial
                    ## concentrations (could use log10(0.01+x))
      )

if (FALSE) {
    library(rgl)
    persp3d(y=as.numeric(colnames(bmat_flip)),
            z=bmat_flip_norm, xlab="time",ylab="beta prob",
            zlim=c(0,0.4), ## restrict scale so not overwhelmed by initial
            ## concentrations (could use log10(0.01+x))
            col="gray"
            )
}

