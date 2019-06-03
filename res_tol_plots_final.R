#######################################################################################
### Plots for both determinisitc and stochastic runs, final output for results for  ###
#######################################################################################

## Build the R0 surface with the parameters of a given run parameter set
eff_scale <- stochas.res_subset$eff_scale[1]

## Function that I use for plotting
eff_calc    <- function (x, y, eff_scale) {
  exp(-(x - y)^2 / eff_scale)
}

############
## Beta and R0 surface
############
mut_link_h <- make.link("log")
mut_link_p <- make.link("logit")
alpha0     <- c(seq(mut_link_p$linkfun(0.01), mut_link_p$linkfun(0.99), length = 100))
tuning     <- c(seq(mut_link_p$linkfun(0.01), mut_link_p$linkfun(0.99), length = 100))
tol0       <- 1
res0       <- 1
startvals  <- calc_startvals_determ(
  alpha0     = alpha0
, tuning     = tuning
, res0       = res0
, tol0       = tol0
, gamma0     = 0.2
, d          = 0.01
, N          = 4000
, power_c    = 2
, power_exp  = 2
, mut_link_h = mut_link_h
, mut_link_p = mut_link_p
, eff_scale  = eff_scale)

beta_surface <- startvals$joint_beta
R0_surface   <- sweep(beta_surface, 2
  , (
#    mut_link_p$linkinv(alpha0 - log(tol0) - log(res0)) + 0.2 + 0.01
   1 - (1 - mut_link_p$linkinv(alpha0 - log(tol0) - log(res0))) * (1 - 0.2) * (1 - 0.01)
    )
  , FUN = "/")
R0_surface_p <- melt(R0_surface)
R0_surface_p <- transform(
    R0_surface_p
  , agg  = rep(alpha0, length(tuning))
  , tune = rep(tuning, each = length(alpha0)))
names(R0_surface_p)[3] <- "R0"

##############################
## Setup single run results ##
##############################

## In case multiple runs were saved for some reason because of batching
  ## for stochas.res_N200 runs 51-79 were done 2x. remove 28 runs
#stochas.res_temp[stochas.res_temp$run > 50 & stochas.res_temp$run < 80, ]$run <- rep(seq(251, (250 + (87000 / 1500)), by = 1), each = 1500)
#stochas.res_temp <- stochas.res_temp %>% filter(stochas.res_temp$run < 279)

######
## Stochastic run
######

## Summarize stochastic run
stochas.res_subset <- transform(stochas.res_subset, R0 = 0)
stochas.res_subset_s <- stochas.res_subset %>%
  group_by(time) %>%
  summarize(
    q05_mean_plalpha = quantile(mean_plalpha, 0.05)
  , q25_mean_plalpha = quantile(mean_plalpha, 0.25)
  , q50_mean_plalpha = quantile(mean_plalpha, 0.50)
  , q75_mean_plalpha = quantile(mean_plalpha, 0.75)
  , q95_mean_plalpha = quantile(mean_plalpha, 0.95)
  , q05_mean_plbeta  = quantile(mean_plbeta, 0.05)
  , q25_mean_plbeta  = quantile(mean_plbeta, 0.25)
  , q50_mean_plbeta  = quantile(mean_plbeta, 0.50)
  , q75_mean_plbeta  = quantile(mean_plbeta, 0.75)
  , q95_mean_plbeta  = quantile(mean_plbeta, 0.95)) %>%
  mutate(R0 = 0)

## Could consider melting and using aes(group = variable) instead to make plot code cleaner
#stochas.res_subset_s <- melt(stochas.res_subset_s, c("time", "R0"))

######
## Deterministic run
######

###########################################################################################
## At one point I thought that the deterministic data frame has the names flipped, but I am
## not so certain about that anymore
###########################################################################################
#names(determ.res_subset)[c(4, 5)] <- c("tune", "agg")

## Summarize deterministic run into the median and quantiles
### ****** IMPORTANT ****** R0 here is not R0, it is the density, it is just called R0 so that
### it can be plotted on the R0 surface
## Also ******** Important ********* recal that plbeta means tune here

#######
## Cleanly calculate the mean and sd for both agg and tune
#######
determ.res_temp_agg <- determ.res_subset %>%
  group_by(Time, agg) %>%
  ## marginal
  summarize(marg_agg = sum(R0)) %>%
  group_by(Time) %>%
  summarize(
    mean   = sum(agg * marg_agg)
  , st_dev = sqrt(sum(agg^2 * marg_agg) - sum(agg * marg_agg)^2)
    )
determ.res_temp_agg <- transform(determ.res_temp_agg, parameter = "agg")

determ.res_temp_tune <- determ.res_subset %>%
  group_by(Time, tune) %>%
  ## marginal
  summarize(marg_tune = sum(R0)) %>%
  group_by(Time) %>%
  summarize(
    mean   = sum(tune * marg_tune)
  , st_dev = sqrt(sum(tune^2 * marg_tune) - sum(tune * marg_tune)^2)
    )
determ.res_temp_tune <- transform(determ.res_temp_tune, parameter = "tune")

mat_out_all_temp_s <- rbind(determ.res_temp_agg, determ.res_temp_tune)
mat_out_all_temp_s <- mat_out_all_temp_s %>% mutate(R0 = 0)

## For now have to adjust this manually.
 ## e.g. add mut_sd 
mat_out_all_temp_s <- transform(mat_out_all_temp_s
  , mut_sd = 0.0625
  )

#######
## ....OR.... Jump through a number of hoops to calculate the median and quantiles
#######

## pleasant way to filter with a series
#years_wanted <- seq(10, 1000, by = 30)
#mat_out_all_temp <- determ.res_subset %>% filter(Time %in% years_wanted)
  
#mat_out_all_temp_s_tune <- mat_out_all_temp %>%
#  group_by(Time, tune) %>%
#  ## marginal
#  summarize(marg_agg = sum(R0)) %>%
#  group_by(Time) %>%
#  ## cumulative sum in order to get quantiles
#  mutate(cum_sum  = cumsum(marg_agg) / max(cumsum(marg_agg))) %>%
#  ## find the quantiles
#  summarize(
#    q05_mean_plalpha = tune[which(cum_sum > 0.05)[1]]
#  , q25_mean_plalpha = tune[which(cum_sum > 0.25)[1]]
#  , q50_mean_plalpha = tune[which(cum_sum > 0.50)[1]]
#  , q75_mean_plalpha = tune[which(cum_sum > 0.75)[1]]
#  , q95_mean_plalpha = tune[which(cum_sum > 0.95)[1]]
#  )

#mat_out_all_temp_s_agg <- mat_out_all_temp %>%
#  group_by(Time, agg) %>%
#  ## marginal
#  summarize(marg_tune = sum(R0)) %>%
#  group_by(Time) %>%
#  ## cumulative sum in order to get quantiles
#  mutate(cum_sum  = cumsum(marg_tune) / max(cumsum(marg_tune))) %>%
#  ## find the quantiles
#  summarize(
#    q05_mean_plbeta = agg[which(cum_sum > 0.05)[1]]
#  , q25_mean_plbeta = agg[which(cum_sum > 0.25)[1]]
#  , q50_mean_plbeta = agg[which(cum_sum > 0.50)[1]]
#  , q75_mean_plbeta = agg[which(cum_sum > 0.75)[1]]
#  , q95_mean_plbeta = agg[which(cum_sum > 0.95)[1]]
#  )
  
#mat_out_all_temp_s <- cbind(mat_out_all_temp_s_agg, mat_out_all_temp_s_tune[, -1])
#mat_out_all_temp_s <- mat_out_all_temp_s %>% mutate(R0 = 0)

######
## Then restrict to a smaller range of times for the R0 surface plots
######

## Have to set this up by hand for each depending on the speed of the deterministic simulation
mat_out_all_temp <- droplevels(subset(determ.res_subset
    , Time == 2 |
      Time == 10 |
      Time == 50 |
      Time == 100 |
      Time == 180 |
      Time == 300 |
      Time == 650 |
      Time == 1000
  ))
  
names(mat_out_all_temp)[3] <- "R0"

############
## Full R0 surface with all results
############

ggright <- ggplot(R0_surface_p, aes(tune, agg, z = R0/max(R0))) +
  scale_fill_distiller(palette = "Greys", na.value = "white") +
 # geom_raster(aes(fill = R0/max(R0))) +
  geom_raster(aes(fill = R0)) +
  stat_contour(bins = 10, colour = "black") +
  xlab("Parasite log(replication rate)") + 
  ylab("Parasite Tuning") + 
  guides(fill = guide_colorbar(title = "Relative R0")) +
  theme(legend.key.size = unit(.45, "cm")) +
  geom_vline(xintercept = mut_link_p$linkfun(.208), lwd = 0.5, linetype = "dashed") +
  geom_vline(xintercept = mut_link_p$linkfun(.26), lwd = 0.5, linetype = "dotted") +  
  geom_path(
     data = stochas.res_subset
     , aes(mean_plalpha, mean_plbeta
      , group = run)
    , lwd  = 0.2, alpha = 0.3
    , colour = "steelblue4") +
  geom_path(
     data = stochas.res_subset_s
   , aes(q05_mean_plalpha, q05_mean_plbeta)
    , lwd  = 0.7, colour = "steelblue3") +
  geom_path(
     data = stochas.res_subset_s
    , aes(q25_mean_plalpha, q25_mean_plbeta)
    , lwd  = 1.0, colour = "steelblue3") +
  geom_path(
     data = stochas.res_subset_s
    , aes(q50_mean_plalpha, q50_mean_plbeta)
    , lwd  = 1.3, colour = "steelblue3") +
  geom_path(
     data = stochas.res_subset_s
    , aes(q75_mean_plalpha, q75_mean_plbeta)
    , lwd  = 1.0, colour = "steelblue3") +
  geom_path(
     data = stochas.res_subset_s
    , aes(q95_mean_plalpha, q95_mean_plbeta)
    , lwd  = 0.7, colour = "steelblue3") +
  stat_contour(
     data = mat_out_all_temp
    , aes(agg, tune, group = Time)
    , bins  = 5
  #  , color = "firebrick3"
    , color = "orchid4"
    , lwd   = 0.8
    , alpha = 1.0) +
  geom_path(
     data = grad_ascent.res_subset
    , aes(alpha, tune
#      , group = run_num
      )
   , lwd  = 2.0
#  , alpha = 0.5
    , colour = "black") +
#  geom_path(
#     data = grad_ascent1
#    , aes(alpha, tune
#      )
#  , lwd  = 2.0
#    , colour = "black") +  
  geom_path(
    data = grad_ascent_sing.res_subset
    , aes(alpha , tune)
    , lwd  = 2.0, colour = "black", linetype = "dashed")

############
## Translating movement on this space into classic tradeoff space
############
  
## The tradeoff curve for ggplots tracking progress on tradeoff space
power_trade_dat <- data.frame(
    alpha = seq(0.01, 1.0, by = 0.01)
  , beta  = power_tradeoff(alpha = seq(0.01, 1.0, by = 0.01), c = 0.75, curv = 2)
  , R0    = power_tradeoff(alpha = seq(0.01, 1.0, by = 0.01), c = 0.75, curv = 2) / ( seq(0.01, 1.0, by = 0.01) + 0.2 + 0.01 )
  )
power_trade_dat_prob <- data.frame(
    alpha = seq(0.01, 1.0, by = 0.01)
  , beta  = power_tradeoff(alpha = seq(0.01, 1.0, by = 0.01), c = 0.75, curv = 2)
  , R0    = power_tradeoff(alpha = seq(0.01, 1.0, by = 0.01), c = 0.75, curv = 2) / 
    (1 - (1 - seq(0.01, 1.0, by = 0.01)) * (1 - 0.20) * (1 - 0.01))
  )
power_trade_dat <- power_trade_dat %>%
  mutate(
    beta_rel = beta / max(beta)
  , R0_rel   = R0 / max(R0)
  )
power_trade_dat_prob <- power_trade_dat_prob %>%
  mutate(
    beta_rel = beta / max(beta)
  , R0_rel   = R0 / max(R0)
  )

######
## If jumping through some hoops for median, probably need to smooth for the plot, 
## so likely just best to use mean and sd
######

## Not super scientific, but smooth to fill in the gaps of the deterministic model because of the
  ## breaks in the time frame
#smoothed_q50_mean_plalpha <- predict(smooth.spline(mat_out_all_temp_s$q50_mean_plalpha, spar = 0.4), x = seq(1, 34, by = 0.5))
#smoothed_q50_mean_plbeta  <- predict(smooth.spline(mat_out_all_temp_s$q50_mean_plbeta, spar = 0.4), x = seq(1, 34, by = 0.5))
#determ_smooth <- data.frame(
#  q50_mean_plalpha = smoothed_q50_mean_plalpha$y
#, q50_mean_plbeta  = smoothed_q50_mean_plbeta$y
#)

#######
## Plot for the virulence beta result
#######

## First make a data frame to house all of the results for efficient plotting 
alpha_beta_plot <- rbind(
  stochas.res_subset_s   %>% dplyr::select(q50_mean_plalpha, q50_mean_plbeta) %>% mutate(model = "DTS")
  ## Is using the smoothed quantiles use this
#, determ_smooth          %>% dplyr::select(q50_mean_plalpha, q50_mean_plbeta) %>% mutate(model = "RD")
  ## otherwise use this
, data.frame(
    q50_mean_plalpha = mat_out_all_temp_s[mat_out_all_temp_s$parameter == "agg", ]$mean
  , q50_mean_plbeta  = mat_out_all_temp_s[mat_out_all_temp_s$parameter == "tune", ]$mean
  , model            = "RD"
)
, with(grad_ascent.res_subset, data.frame(q50_mean_plalpha = alpha, q50_mean_plbeta = tune, model = "AD"))
)
alpha_beta_plot <- transform(alpha_beta_plot, model = factor(model, levels = c("DTS", "RD", "AD")))

alpha_beta_plot_t <- alpha_beta_plot %>% filter(model == "AD")

ggtop <- ggplot(power_trade_dat, aes(alpha, beta_rel)) + 
  geom_line(lwd = 0.5, linetype = "dashed") + 
    geom_path(
    data = alpha_beta_plot
  , aes(
    plogis(q50_mean_plalpha)
  , power_tradeoff(plogis(q50_mean_plalpha), 0.75, 2) * exp(-(q50_mean_plalpha - q50_mean_plbeta)^2 / eff_scale) / max(power_trade_dat$beta)
     , colour = model 
    )
    , lwd    = 1
 #   , lwd    = 0
    ) +
  geom_vline(xintercept = .21, lwd = 0.5, linetype = "dashed") +
  geom_vline(xintercept = .26, lwd = 0.5, linetype = "dotted") +  
    scale_colour_manual(
       values = 
        c("steelblue3", "firebrick3", "black")
  #      c("white")
      , labels = 
        c("Discrete Time Stochastic", "Reaction Diffusion", "Adaptive Dynamics")
  #       c("Adaptive Dynamics")
      , name  = "Model") +
#  guides(colour = FALSE, linetype = FALSE) +
  theme(
    legend.position = c(0.65, 0.25)
  , legend.key.size = unit(.60, "cm")
  , legend.title = element_text(size = 14)
  , legend.text  = element_text(size = 12)) +
    xlab("") + ylab("Parasite Relative Beta")
    
ggbottom <- ggplot(power_trade_dat, aes(alpha, R0_rel)) + 
  geom_line(lwd = 0.5, linetype = "dashed") + 
  geom_line(data = power_trade_dat_prob, lwd = 0.5, linetype = "dotted") +
    geom_path(
    data = alpha_beta_plot
  , aes(
    plogis(q50_mean_plalpha)
  , (power_tradeoff(plogis(q50_mean_plalpha), 0.75, 2) * exp(-(q50_mean_plalpha - q50_mean_plbeta)^2 / eff_scale) /
      (plogis(q50_mean_plalpha) + 0.2 + 0.01)) / max(power_trade_dat$R0)
    , colour = model
    )
    , lwd    = 1
   # , lwd    = 0
    ) +
  geom_vline(xintercept = .21, lwd = 0.5, linetype = "dashed") +
  geom_vline(xintercept = .26, lwd = 0.5, linetype = "dotted") +  
    scale_colour_manual(
       values = 
        c("steelblue3", "firebrick3", "black")
    #    c("black")
      , labels = 
        c("Discrete Time Stochastic", "Reaction Diffusion", "Adaptive Dynamics")
   #      c("Adaptive Dynamics")
      ) +
    theme(legend.position = "none") +
    xlab("Parasite Virulence") + ylab("Parasite Relative R0")

# grid.arrange(ggtop, ggbottom, ncol = 1)

grid.arrange(
  grobs         = list(ggtop, ggbottom, ggright)
#, widths        = c(1, 1, 2)
, layout_matrix = rbind(c(1, 3, 3), c(2, 3, 3))
  )

########################################################
## Results across parameter values, not on R0 surface ## 
########################################################

#########
### Stochastic results
#########

## Summary of all runs
stochas.res_temp_s <- stochas.res_subset %>%
  group_by(time
    , mut_sd
    ) %>%
  summarize(
    q50_mean_plalpha = quantile(mean_plalpha, 0.50)
  , q50_mean_plbeta  = quantile(mean_plbeta, 0.50)
  , q50_sd_plalpha   = quantile(sd_plalpha, 0.50)
  , q50_sd_plbeta    = quantile(sd_plbeta, 0.50)
  , sd_mean_plalpha  = sd(mean_plalpha, 0.50)
  , sd_mean_plbeta   = sd(mean_plbeta, 0.50)) %>%
  mutate(R0 = 0)

stochas.res_temp_s <- transform(stochas.res_temp_s,
  mut_sd = factor(mut_sd, levels = c("0.25", "0.05", "0.01"))
  )

######
## Stochastic plot over a single parameter
######

unique_vals <- unique(stochas.res_temp_s$mut_sd)
unique_cols <- c("aquamarine4", "coral4", "cyan2")
unique_lins <- c("solid", "dashed", "dotted")

## Temorary deterministic plot to add to the stochastic plot
temp_determ_plot <- mat_out_all_temp_s[mat_out_all_temp_s$parameter == "agg", ]
#temp_determ_plot$time <- seq(min(unique(stochas.res_temp_s$time)), max(unique(stochas.res_temp_s$time)), length = nrow(temp_determ_plot))
temp_determ_plot$time <- temp_determ_plot$Time * 500
names(temp_determ_plot)[2] <- "q50_mean_plalpha"

ggplot(stochas.res_temp_s, aes(time, plogis(q50_mean_plalpha))) + 
## Not sure how to do this with colour and linetype. Just do one at a time
  geom_path(
    data = stochas.res_subset[stochas.res_subset$mut_sd == unique_vals[3], ]
  , aes(time, plogis(mean_plalpha) 
  , group    = as.factor(run)
    )  
    , alpha = 0.075
    , colour   = unique_cols[1]
 #   , linetype = unique_lins[1]
    ) +
  geom_path(
    data = stochas.res_subset[stochas.res_subset$mut_sd == unique_vals[2], ]
  , aes(time, plogis(mean_plalpha) 
  , group    = as.factor(run)
    )  
    , alpha = 0.075
    , colour   = unique_cols[2]
#    , linetype = unique_lins[3]
    ) +
  geom_path(
    data = stochas.res_subset[stochas.res_subset$mut_sd == unique_vals[1], ]
  , aes(time, plogis(mean_plalpha) 
  , group    = as.factor(run)
    )  
    , alpha = 0.075
    , colour   = unique_cols[3]
 #   , linetype = unique_lins[2]
    ) +
  geom_path(aes(
    colour   = as.factor(mut_sd)
#  , linetype = as.factor(mut_sd)
    )
  , lwd = 1) +
  scale_colour_manual(
    values = unique_cols
  , name   = "Mutational Standard 
Deviation") +
#  scale_linetype_manual(
#    values = unique_lins
#  , name   = "Mutational Standard 
#Deviation") +
  ## Add the deterministic result
  geom_path(
    data = temp_determ_plot
  , aes(time, plogis(q50_mean_plalpha))) +
  geom_ribbon(
    data = temp_determ_plot
  , aes(time
    , ymin = plogis(q50_mean_plalpha - st_dev)
    , ymax = plogis(q50_mean_plalpha + st_dev))
  , alpha = 0.3
  ) +
  scale_x_log10(labels=trans_format('log10',math_format(10^.x))) +
  xlab("Time (steps)") + 
  ylab("Median parasite virulence") +
  theme(
    legend.position = c(0.85, 0.88)
  , legend.key.size = unit(.60, "cm")
  , legend.title = element_text(size = 14)
  , legend.text  = element_text(size = 12))
  
#########
### Determinisitc results
#########

## Requires mat_out_all_temp_s created earlier in this script

## Thought that these were flipped, but not confident: Eeek! The deterministic data frame has the names flipped!!!
# names(determ.res_temp)[c(4, 5)] <- c("tune", "agg")

#determ_smooth <- data.frame(
#    q05_mean_plalpha = predict(smooth.spline(mat_out_all_temp_s$q05_mean_plalpha, spar = 0.4), x = seq(1, 34, by = 0.5))$y
#  , q25_mean_plalpha = predict(smooth.spline(mat_out_all_temp_s$q25_mean_plalpha, spar = 0.4), x = seq(1, 34, by = 0.5))$y
#  , q50_mean_plalpha = predict(smooth.spline(mat_out_all_temp_s$q50_mean_plalpha, spar = 0.4), x = seq(1, 34, by = 0.5))$y
#  , q75_mean_plalpha = predict(smooth.spline(mat_out_all_temp_s$q75_mean_plalpha, spar = 0.4), x = seq(1, 34, by = 0.5))$y
#  , q95_mean_plalpha = predict(smooth.spline(mat_out_all_temp_s$q95_mean_plalpha, spar = 0.4), x = seq(1, 34, by = 0.5))$y
#  , q05_mean_plbeta  = predict(smooth.spline(mat_out_all_temp_s$q05_mean_plbeta, spar = 0.4), x = seq(1, 34, by = 0.5))$y
#  , q25_mean_plbeta  = predict(smooth.spline(mat_out_all_temp_s$q25_mean_plbeta, spar = 0.4), x = seq(1, 34, by = 0.5))$y
#  , q50_mean_plbeta  = predict(smooth.spline(mat_out_all_temp_s$q50_mean_plbeta, spar = 0.4), x = seq(1, 34, by = 0.5))$y
#  , q75_mean_plbeta  = predict(smooth.spline(mat_out_all_temp_s$q75_mean_plbeta, spar = 0.4), x = seq(1, 34, by = 0.5))$y
#  , q95_mean_plbeta  = predict(smooth.spline(mat_out_all_temp_s$q95_mean_plbeta, spar = 0.4), x = seq(1, 34, by = 0.5))$y
#  )

## go back to determ tidy and reload the needed deterministic run
#determ_smooth <- transform(determ_smooth
#  , Time = seq(mat_out_all_temp_s$Time[1], mat_out_all_temp_s$Time[length(mat_out_all_temp_s$Time)], by = 15)
#  , mu   = rep(determ.res_subset[[3]]$mu, nrow(determ_smooth))
#)

## Need to calculate sd of the RD distribution
## **** For now copying code, but should add to the place far above where I calculate the
 ## mean of the distribution

## Plot just the deterministic results
ggplot(mat_out_all_temp_s, aes(Time, mean)) + 
  geom_ribbon(aes(
    Time
    , ymax = mean + st_dev
    , ymin = mean - st_dev
    , fill     = as.factor(parameter)
    )
    , alpha = 0.2) +
  geom_path(
    aes(
      linetype = as.factor(parameter)
    , colour   = as.factor(parameter)
      )
   ) + 
  scale_colour_manual(
    values = unique_cols[1:2]
  , name   = "Parasite Trait"
  , labels = c("log(Replication rate)", "Tuning")) +
  scale_linetype_manual(
  values = unique_lins[1:2]
  , name   = "Parasite Trait"
  , labels = c("log(Replication rate)", "Tuning")) +
  scale_fill_manual(
  values = unique_cols[1:2]
  , name   = "Parasite Trait"
  , labels = c("log(Replication rate)", "Tuning")) +
  xlab("Time steps") + 
  ylab("Parasite trait value") +
  theme(
  #  legend.position = c(750, 2)
    legend.key.size = unit(.60, "cm")
  , legend.title = element_text(size = 14)
  , legend.text  = element_text(size = 12))

ggplot(mat_out_all_temp_s, aes(Time, st_dev)) + 
  geom_path(
    aes(
      linetype = as.factor(parameter)
    , colour   = as.factor(parameter))
    ) + 
  scale_colour_manual(
    values = unique_cols[1:2]
  , name   = "Parasite Trait"
  , labels = c("log(Replication rate)", "Tuning")) +
    scale_linetype_manual(
  values = unique_lins[1:2]
  , name   = "Parasite Trait"
  , labels = c("log(Replication rate)", "Tuning")) +
  xlab("Time steps)") + 
  ylab("SD in parasite virulence")


#########
### Other summaries of stochastic results for second section of the results
#########

## check what subset I am working with currently
with(stochas.res_subset, list(mu = unique(mu), mut_sd = unique(mut_sd), alpha0 = unique(alpha0), tune0 = unique(tune0), N = unique(N)))

## (1) Variation in trajectory: q50_sd_plalpha (also could consider showing on the tradeoff curve?)
## (2) Time to equilibrium: first time q50_mean_plalpha == 0.21
## (3) Transient virulence departure: maximum q50_mean_plalpha and cumulative time above q50_mean_plalpha == 0.21
## (4) Fluctuations from equilibrium: size of q50_sd_plalpha after first time q50_mean_plalpha == 0.21

#######
## Summary of all parameter values all runs
#######

stochas.res_temp_s <- stochas.res_subset %>%
  group_by(
      time
    , mu
    , mut_sd
    , N
    , eff_scale
    ) %>%
  summarize(
    q50_mean_plalpha = quantile(mean_plalpha, 0.50)
  , q50_mean_plbeta  = quantile(mean_plbeta, 0.50)
  , q50_sd_plalpha   = quantile(sd_plalpha, 0.50)
  , q50_sd_plbeta    = quantile(sd_plbeta, 0.50)
  , q50_sd_alpha     = quantile(sd_alpha, 0.50)
  , q50_sd_beta      = quantile(sd_beta, 0.50)    
  , sd_mean_plalpha  = sd(mean_plalpha, 0.50)
  , sd_mean_plbeta   = sd(mean_plbeta, 0.50)
  , num_I_strains_05 = quantile(num_I_strains, 0.05)  
  , num_I_strains_50 = quantile(num_I_strains, 0.50)
  , num_I_strains_95 = quantile(num_I_strains, 0.95) 
    ) %>%
  mutate(R0 = 0)

stochas.res_temp_ss <- stochas.res_temp_s %>%
   group_by(
      mu
    , mut_sd
    , N
    , eff_scale
    ) %>%
  summarize(
## (2) Time to equilibrium: first time q50_mean_plalpha == opt_vir (somewhat hard because it passes over and comes back. Use tune)
    tte     = min(which((q50_mean_plalpha < mut_link_p$linkfun(0.26)) & (q50_mean_plbeta < 0)))
  , 
## (2.2) Also need to check the first time that it passes by, to calculate the width of the peak (time above equilibrium)
    tte_first     = min(which((q50_mean_plalpha > mut_link_p$linkfun(0.26))))
  ,     
## (3) Transient virulence departure: maximum q50_mean_plalpha and cumulative time above q50_mean_plalpha == 0.21
    max_vir = max(q50_mean_plalpha)
  , cum_vir = length(which(q50_mean_plalpha > mut_link_p$linkfun(0.26)))
  , num_I_strains = median(num_I_strains_50))

## check what the unique times are for each of these parameter combinations, then convert tte to an actual time
stochas.res_temp_s_times <- stochas.res_temp_s %>%
   group_by(
      mu
    , mut_sd
    , N
    , eff_scale
    ) %>%
  summarize(
    max_time = max(time)
  , numtimes = n())

## update tte
stochas.res_temp_ss <- left_join(stochas.res_temp_ss, stochas.res_temp_s_times, by = c("mu", "mut_sd", "N", "eff_scale"))
stochas.res_temp_ss <- stochas.res_temp_ss %>%
  group_by(mu, mut_sd, N, eff_scale) %>%
  mutate(
    tte_time       = seq(1, max_time, length = numtimes)[tte]
  , tte_first_time = seq(1, max_time, length = numtimes)[tte_first]
  , cum_vir_time   = cum_vir * max_time / numtimes
    ) %>%
  ## convert to generations of the host (dividing by natural death rate of the host, which in this model is 
   ## also the birth rate (at least for now))
  mutate(
    tte_time_gen     = tte_time / (1/0.01)
  , tte_first_time_gen  = tte_first_time / (1/0.01)
  , cum_vir_time_gen = cum_vir_time / (1/0.01)
  )

## (4) Fluctuations from equilibrium: size of q50_sd_plalpha after first time q50_mean_plalpha == 0.21
  ## after the adjustments above, this can be calculated using the first time to equilibrium

## First have to left_join the time to equilibrium value, then filter the times, then summarize
stochas.res_temp_s <- left_join(stochas.res_temp_s, stochas.res_temp_ss[, c(1, 2, 3, 4, 12, 13) ], by = c("mu", "mut_sd", "N", "eff_scale"))
## For plots where I am interested in a reduced time period
stochas.res_temp_s_red <- stochas.res_temp_s %>% filter(time > tte_time)

stochas.res_temp_ss2 <- stochas.res_temp_s_red %>%
  group_by(mu, mut_sd, N, eff_scale) %>%
 summarize(
    from_equil_q50_sd_plalpha = quantile(q50_sd_plalpha, 0.50)
  , from_equil_q50_sd_plbeta  = quantile(q50_sd_plbeta, 0.50)
    )

stochas.res_temp_ss <- left_join(stochas.res_temp_ss, stochas.res_temp_ss2, by = c("mu", "mut_sd", "N", "eff_scale"))

######
### Plot and table for reporting 1:4
######

## Table from stochas.res_temp_ss
## Plots from stochas.res_temp_s, annotated with stochas.res_temp_ss
stochas.res_temp_s <- transform(stochas.res_temp_s,
  mut_sd    = factor(mut_sd, levels = c("0.25", "0.05", "0.01"))
, eff_scale = factor(eff_scale, levels = c("30", "10", "50")))

unique_vals <- unique(stochas.res_temp_s$mut_sd)
unique_cols <- c("aquamarine4", "coral4", "cyan2")
unique_lins <- c("solid", "dashed", "dotted")

ggplot(stochas.res_temp_s, aes(time, plogis(q50_mean_plalpha))) + 
  geom_path(aes(
    colour   = as.factor(mut_sd)
  , linetype = as.factor(eff_scale)
    )
  , lwd = 1) +
  geom_ribbon(
    aes(time
  , ymin   = plogis(q50_mean_plalpha - sd_mean_plalpha * 2)
  , ymax   = plogis(q50_mean_plalpha + sd_mean_plalpha * 2)
## Post-hoc wrong way to do this, but an _ok_ :( placeholder for now
#   , ymin   = plogis(q50_mean_plalpha) - q50_sd_alpha 
#   , ymax   = plogis(q50_mean_plalpha) + q50_sd_alpha  
  , colour = as.factor(mut_sd)
  , fill   = as.factor(mut_sd)
    )  
    , alpha = 0.2, lwd = 0) +
  scale_colour_manual(
    values = unique_cols
  , name   = "Mutational Standard 
Deviation") +
  scale_fill_manual(
    values = unique_cols
  , name   = "Mutational Standard 
Deviation") +
  scale_linetype_manual(
    values = unique_lins
  , name   = "Efficiency Scale") +
  guides(colour = FALSE, linetype = FALSE) +
  ## Add the deterministic result
  scale_x_log10(labels=trans_format('log10',math_format(10^.x))) +
  xlab("Time (steps)") + 
  ylab("Median parasite virulence 
(with variation in median among stochastic runs)") +
#  ylab("Median parasite virulence 
#(with variation in virulence among strains within runs)") +
  facet_grid(N ~ mu) +
  theme(
    legend.position = c(0.85, 0.83)
  , legend.key.size = unit(.50, "cm")
  , legend.title = element_text(size = 20)
  , legend.text  = element_text(size = 20)
  , axis.text.x = element_text(size = 24)
  , axis.text.y = element_text(size = 24)
  , axis.title.x = element_text(size = 24)
  , axis.title.y = element_text(size = 24)
  , panel.margin.x = unit(0.5, "lines")
  , strip.text.x = element_text(size = 24)
  , strip.text.y = element_text(size = 24)) 
