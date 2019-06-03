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
R0_surface   <- sweep(beta_surface, 2, (mut_link_p$linkinv(
    alpha0 - log(tol0) - log(res0)) +
    0.2 + 0.01), FUN = "/")
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

## Temp R0 column for plotting on the surface. May actually want to plot calculated R0 later
determ.res_subset  <- determ.res_subset[[1]]
###########################################################################################
## At one point I thought that the deterministic data frame has the names flipped, but I am
## not so certain about that anymore
###########################################################################################
#names(determ.res_subset)[c(4, 5)] <- c("tune", "agg")

########
## First gather quantiles for a larger range of times
########
names(determ.res_subset)[3] <- "R0"

## Summarize deterministic run into the median and quantiles
### ****** IMPORTANT ****** R0 here is not R0, it is the density, it is just called R0 so that
### it can be plotted on the R0 surface
## Also ******** Important ********* recal that plbeta means tune here

## pleasant way to filter with a series
years_wanted <- seq(10, 1000, by = 30)
mat_out_all_temp <- determ.res_subset %>% filter(Time %in% years_wanted)
  
mat_out_all_temp_s_tune <- mat_out_all_temp %>%
  group_by(Time, tune) %>%
  ## marginal
  summarize(marg_agg = sum(R0)) %>%
  group_by(Time) %>%
  ## cumulative sum in order to get quantiles
  mutate(cum_sum  = cumsum(marg_agg) / max(cumsum(marg_agg))) %>%
  ## find the quantiles
  summarize(
    q05_mean_plalpha = tune[which(cum_sum > 0.05)[1]]
  , q25_mean_plalpha = tune[which(cum_sum > 0.25)[1]]
  , q50_mean_plalpha = tune[which(cum_sum > 0.50)[1]]
  , q75_mean_plalpha = tune[which(cum_sum > 0.75)[1]]
  , q95_mean_plalpha = tune[which(cum_sum > 0.95)[1]]
  )

mat_out_all_temp_s_agg <- mat_out_all_temp %>%
  group_by(Time, agg) %>%
  ## marginal
  summarize(marg_tune = sum(R0)) %>%
  group_by(Time) %>%
  ## cumulative sum in order to get quantiles
  mutate(cum_sum  = cumsum(marg_tune) / max(cumsum(marg_tune))) %>%
  ## find the quantiles
  summarize(
    q05_mean_plbeta = agg[which(cum_sum > 0.05)[1]]
  , q25_mean_plbeta = agg[which(cum_sum > 0.25)[1]]
  , q50_mean_plbeta = agg[which(cum_sum > 0.50)[1]]
  , q75_mean_plbeta = agg[which(cum_sum > 0.75)[1]]
  , q95_mean_plbeta = agg[which(cum_sum > 0.95)[1]]
  )
  
mat_out_all_temp_s <- cbind(mat_out_all_temp_s_agg, mat_out_all_temp_s_tune[, -1])
mat_out_all_temp_s <- mat_out_all_temp_s %>% mutate(R0 = 0)

######
## Then restrict to a smaller range of times for the R0 surface plots
######

## Have to set this up by hand for each depending on the speed of the deterministic simulation
mat_out_all_temp <- droplevels(subset(determ.res_subset
    , Time == 2 |
      Time == 10 |
      Time == 30 |
      Time == 50 |
      Time == 100 |
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
  geom_raster(aes(fill = R0/max(R0))) +
  stat_contour(bins = 10, colour = "black") +
  xlab("Parasite log(replication rate)") + 
  ylab("Parasite Tuning") + 
  guides(fill = guide_colorbar(title = "Relative R0")) +
  theme(legend.key.size = unit(.45, "cm")) +
  geom_vline(xintercept = mut_link_p$linkfun(.208), lwd = 0.5, linetype = "dotted") +
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
    , aes(group = Time)
    , bins  = 5
    , color = "firebrick3"
    , lwd   = 0.8
    , alpha = 1.0) +
  geom_path(
     data = grad_ascent.res_subset
    , aes(alpha, tune)
    , lwd  = 2.0, colour = "black") + 
  geom_path(
    data = grad_ascent_sing.res_subset
    , aes(alpha , tune)
    , lwd  = 2.0, colour = "black", linetype = "dashed")

############
## Translating movement on this space into classic tradeoff space
############
  #### ggplots tracking progress on tradeoff space
  power_trade_dat <- data.frame(
    alpha = seq(0.01, 1.0, by = 0.01)
  , beta  = power_tradeoff(alpha = seq(0.01, 1.0, by = 0.01), c = 0.75, curv = 2)
  , R0    = power_tradeoff(alpha = seq(0.01, 1.0, by = 0.01), c = 0.75, curv = 2) / ( seq(0.01, 1.0, by = 0.01) + 0.2 + 0.01 )
  )
power_trade_dat <- power_trade_dat %>%
  mutate(
    beta_rel = beta / max(beta)
  , R0_rel   = R0 / max(R0)
  )

## Not super scientific, but smooth to fill in the gaps of the deterministic model because of the
  ## breaks in the time frame
smoothed_q50_mean_plalpha <- predict(smooth.spline(mat_out_all_temp_s$q50_mean_plalpha, spar = 0.4), x = seq(1, 34, by = 0.5))
smoothed_q50_mean_plbeta  <- predict(smooth.spline(mat_out_all_temp_s$q50_mean_plbeta, spar = 0.4), x = seq(1, 34, by = 0.5))
determ_smooth <- data.frame(
  q50_mean_plalpha = smoothed_q50_mean_plalpha$y
, q50_mean_plbeta  = smoothed_q50_mean_plbeta$y
)

#######
## Plot for the virulence beta result
#######

## First make a data frame to house all of the results for efficient plotting 
alpha_beta_plot <- rbind(
  stochas.res_subset_s   %>% dplyr::select(q50_mean_plalpha, q50_mean_plbeta) %>% mutate(model = "DTS")
, determ_smooth          %>% dplyr::select(q50_mean_plalpha, q50_mean_plbeta) %>% mutate(model = "RD")
, with(grad_ascent.res_subset, data.frame(q50_mean_plalpha = alpha, q50_mean_plbeta = tune, model = "AD"))
)
alpha_beta_plot <- transform(alpha_beta_plot, model = factor(model, levels = c("DTS", "RD", "AD")))

ggtop <- ggplot(power_trade_dat, aes(alpha, beta_rel)) + 
  geom_line(lwd = 0.5, linetype = "dotted") + 
    geom_path(
    data = alpha_beta_plot
  , aes(
    plogis(q50_mean_plalpha)
  , power_tradeoff(plogis(q50_mean_plalpha), 0.75, 2) * exp(-(q50_mean_plalpha - q50_mean_plbeta)^2 / eff_scale) / max(power_trade_dat$beta)
     , colour = model 
    )
    , lwd    = 1
    ) +
    scale_colour_manual(values = c("steelblue3", "firebrick3", "black")
      , labels = c("Discrete Time Stochastic", "Reaction Diffusion", "Adaptive Dynamics")
      , name  = "Model") +
  theme(
    legend.position = c(0.65, 0.25)
  , legend.key.size = unit(.60, "cm")
  , legend.title = element_text(size = 14)
  , legend.text  = element_text(size = 12)) +
    xlab("") + ylab("Parasite Relative Beta")
    
ggbottom <- ggplot(power_trade_dat, aes(alpha, R0_rel)) + 
  geom_line(lwd = 0.5, linetype = "dotted") + 
    geom_path(
    data = alpha_beta_plot
  , aes(
    plogis(q50_mean_plalpha)
  , (power_tradeoff(plogis(q50_mean_plalpha), 0.75, 2) * exp(-(q50_mean_plalpha - q50_mean_plbeta)^2 / eff_scale) /
      (plogis(q50_mean_plalpha) + 0.2 + 0.01)) / max(power_trade_dat$R0)
    , colour = model
    )
    , lwd    = 1
    ) +
    scale_colour_manual(values = c("steelblue3", "firebrick3", "black")
      , labels = c("Discrete Time Stochastic", "Reaction Diffusion", "Adaptive Dynamics")) +
    theme(legend.position = "none") +
    xlab("Parasite Virulence") + ylab("Parasite Relative R0")

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

## Pick a single subset to plot
stochas.res_temp <- stochas.res_Ncheck

stochas.res_temp_s <- stochas.res_temp %>%
  group_by(time, N) %>%
  summarize(
    q50_mean_plalpha = quantile(mean_plalpha, 0.50)
  , q50_mean_plbeta  = quantile(mean_plbeta, 0.50)
  , q50_sd_plalpha   = quantile(sd_plalpha, 0.50)
  , q50_sd_plbeta    = quantile(sd_plbeta, 0.50)
  , sd_mean_plalpha  = sd(mean_plalpha, 0.50)
  , sd_mean_plbeta   = sd(mean_plbeta, 0.50)) %>%
  mutate(R0 = 0)

ggplot(stochas.res_temp_s, aes(time / nrpt, plogis(q50_mean_plalpha))) + 
  geom_line(aes(
    colour   = as.factor(N)
  , linetype = as.factor(N)
    )) + 
  xlab("Time (steps)") + 
  ylab("Median parasite virulence among all strains across all runs")

ggplot(stochas.res_temp_s[stochas.res_temp_s$time > 300000, ], aes(time / nrpt, plogis(q50_mean_plalpha))) + 
  geom_line(aes(
    colour   = as.factor(N)
  , linetype = as.factor(N)
    )) + 
  xlab("Time (steps)") + 
  ylab("Median parasite virulence among all strains across all runs")

ggplot(stochas.res_temp_s, aes(time / nrpt, q50_sd_plalpha)) + 
  geom_line(aes(
    colour   = as.factor(N)
  , linetype = as.factor(N)
    )) +
  xlab("Time (steps)") + 
  ylab("Median SD in parasite virulence among all strains across all runs")

ggplot(stochas.res_temp_s[stochas.res_temp_s$time > 300000, ], aes(time / nrpt, q50_sd_plalpha)) + 
  geom_line(aes(
    colour   = as.factor(N)
  , linetype = as.factor(N)
    )) +
  xlab("Time (steps)") + 
  ylab("Median SD in parasite virulence among all strains across all runs")

ggplot(stochas.res_temp_s, aes(time / nrpt, sd_mean_plalpha)) + 
  geom_line(aes(
    colour   = as.factor(N)
  , linetype = as.factor(N)
    )) +
  xlab("Time (steps)") + 
  ylab("SD in mean parasite virulence across all runs")

#########
### Determinisitc results
#########

determ.res      <- readRDS("res_out/res_out_RD/res_1000_all_determ_6.Rds")
determ.res_temp <- determ.res[[1]]
## Eeek! The deterministic data frame has the names flipped!!!
names(determ.res_temp)[c(4, 5)] <- c("tune", "agg")
determ.res_temp_s <- determ.res_temp %>% 
  group_by(Time) %>%
  summarize(
    mean_agg = sum(agg * value)
  , sd_agg   = sd(agg * value))
determ.res_temp_s <- transform(determ.res_temp_s
  , mu = determ.res[[3]]$mu)

# Combine some of the deterministic runs with only a single different parameter value
determ.res_temp_s_tot <- determ.res_temp_s
determ.res_temp_s_tot <- rbind(determ.res_temp_s_tot, determ.res_temp_s)

## Plot just the deterministic results
ggplot(determ.res_temp_s, aes(Time, mean_agg)) + 
  geom_line(aes(linetype = as.factor(mu)), lwd = 1) + 
  xlab("Time (steps)") + 
  ylab("Median parasite virulence") 

ggplot(determ.res_temp_s_tot, aes(Time, sd_agg)) + 
  geom_line(lwd = 1, aes(linetype = as.factor(mu))) + 
  xlab("Time (steps)") + 
  ylab("SD in parasite virulence")

## Put the deterministic result together with the stochastic result
ggplot(stochas.res_temp_s, aes(time / nrpt, plogis(q50_mean_plalpha))) + 
  geom_line(aes(colour = as.factor(mut_sd), linetype = as.factor(mut_sd))) + 
  geom_line(data = determ.res_temp_s_tot, aes(Time/5, plogis(mean_agg), linetype = as.factor(mu))) + 
  xlab("Time (steps)") + 
  ylab("Median parasite virulence among all strains across all runs")

  
