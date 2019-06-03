########################################################
### Plots for both determinisitc and stochastic runs ###
########################################################

############
## Beta and R0 surface
############
mut_link_h <- make.link("log")
mut_link_p <- make.link("logit")
#alpha0    <- c(seq(0.01, 0.99, by = 0.01), 0.999)
#tuning    <- c(seq(0.01, 0.99, by = 0.01), 0.999)
alpha0     <- c(seq(
        mut_link_p$linkfun(0.01), mut_link_p$linkfun(0.99)
        , length = 100))
tuning     <- c(seq(
        mut_link_p$linkfun(0.01), mut_link_p$linkfun(0.99)
        , length = 100))
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
, power_c    = 8
, power_exp  = 2
, mut_link_h = mut_link_h
, mut_link_p = mut_link_p
, eff_scale  = 30)

beta_surface <- startvals$joint_beta
R0_surface   <- sweep(beta_surface, 2, (mut_link_p$linkinv(
#  mut_link_p$linkfun(alpha0) - log(tol0) - log(res0)) +
    alpha0 - log(tol0) - log(res0)) +
    0.2 + 0.01), FUN = "/")
R0_surface_p <- melt(R0_surface)
R0_surface_p <- transform(
    R0_surface_p
  , agg  = rep(alpha0, length(tuning))
  , tune = rep(tuning, each = length(alpha0)))
names(R0_surface_p)[3] <- "R0"

### Plot just the surface  
ggplot(R0_surface_p, aes(
  tune
, agg
  , z = R0
  )) +
  scale_fill_distiller(palette = "Spectral", na.value = "white") +
  geom_raster(aes(
 #   fill = R0/max(R0)
     fill = R0
    )) +
  stat_contour() +
  xlab("Parasite Aggressiveness") + 
  ylab("Parasite Tuning") + 
  guides(fill = guide_colorbar(title = "Relative R0"))


################################################################################################
## Stochastic ##################################################################################
################################################################################################

if (deterministic == FALSE) {

############
## Array of plots from a single run of the stochastic model
############
source("res_tol_stochsatic_sing_run.R")

############
## Plot only a subset of the runs from the stochastic model
############
  
## Add a blank column so that the datapoints will show up on the surface 
res_1000_all     <- transform(res_1000_all, R0 = 0)
res_1000_stochas <- res_1000_all %>% filter(alpha0 == 0.05, mut_sd == 0.05, R0_init == 2)
res_1000_stochas <- droplevels(res_1000_stochas)
  
############
## Movement of the stochastic results on the surface
############
  
ggplot(R0_surface_p, aes(Var1, Var2, z = R0/max(R0))) +
  scale_fill_distiller(palette = "Spectral", na.value = "white") +
  geom_raster(aes(fill = R0/max(R0))) +
  stat_contour() +
  xlab("Parasite Aggressiveness") + 
  ylab("Parasite Tuning") + 
  guides(fill = guide_colorbar(title = "Relative R0")) +
  scale_x_continuous(
      breaks = seq(1, length(alpha0), by = 10)
    , labels = alpha0[seq(1, length(alpha0), by = 10)]) +
  scale_y_continuous(
      breaks = seq(1, length(tuning), by = 10)
    , labels = plogis(tuning)[seq(1, length(tuning), by = 10)]) +
  theme(legend.key.size = unit(.45, "cm")) +
  geom_point(
    data = res_1000_stochas[seq(1, nrow(res_1000_stochas), by = 5), ]
    , aes(plogis(mean_plalpha) * 100, plogis(mean_plbeta) * 100)
    , shape = 19
    , lwd  = 1.0) +
  geom_path(
     data = res_1000_stochas[seq(1, nrow(res_1000_stochas), by = 1), ]
    , aes(plogis(mean_plalpha) * 100, plogis(mean_plbeta) * 100)
    , lwd  = 0.5)

} else {

################################################################################################
## Deterministic ###############################################################################
################################################################################################

############
## Movement on the surface of a single deterministic run
## Run this after each deterministic run to build the single data frame of contours
############
  
## plot on the R0 surface or not?
on_R0_surface <- TRUE
do_plot       <- FALSE

determ_med_R0 <- data.frame(
  Time     = seq(1, nrow(res_1000))
, med_R0   = 0
, dist_R0  = 0)
for (i in seq(1, nrow(res_1000), by = 1)) {
  
  mat_out <- matrix(data = c(unlist(res_1000[i, -c(1,2)]))
    , nrow = length(tuning)
    , ncol = length(alpha0)
    , byrow = FALSE
 #  , byrow = TRUE
    )
  determ_med_R0[i, ]$dist_R0 <- sum(mat_out * R0_surface)
  determ_med_R0[i, ]$med_R0  <- R0_surface[which(mat_out == max(mat_out), arr.ind = TRUE)]
  
  mat_out <- melt(mat_out) 
  mat_out <- transform(
    mat_out
  , agg  = rep(alpha0, length(tuning))
  , tune = rep(tuning, each = length(alpha0))
  , Time = res_1000$time[i])
  
  if (i == 1) {
    mat_out_all <- mat_out
  } else {
    mat_out_all <- rbind(mat_out_all, mat_out)
  }

if (do_plot == TRUE) {
    
if (on_R0_surface == FALSE) {
  
names(mat_out)[3] <- "Proportion"
  
ggplot(mat_out, aes(Var1, Var2, z = Proportion)) +
  stat_contour(geom = "polygon", aes(fill = Proportion), color = "white", alpha = 0.25) + 
  scale_fill_distiller(palette = "Spectral", na.value = "white") +
#  stat_contour(geom = "polygon", aes(fill = grad)) + 
  geom_tile(aes(fill = Proportion)) + 
  stat_contour(bins = 15) + 
  xlab("Parasite Aggressiveness") + 
  ylab("Parasite Tuning") + 
  guides(fill = guide_colorbar(title = "Proportion")) +
  scale_x_continuous(
      breaks = seq(1, length(alpha0), by = 10)
    , labels = alpha0[seq(1, length(alpha0), by = 10)]) +
  scale_y_continuous(
      breaks = seq(1, length(tuning), by = 10)
    , labels = plogis(tuning)[seq(1, length(tuning), by = 10)]) +
  theme(legend.key.size = unit(.45, "cm"))

} else {
  
  names(mat_out)[3] <- "R0"
  
  ggplot(R0_surface_p, aes(Var1, Var2, z = R0/max(R0))) +
  scale_fill_distiller(palette = "Spectral", na.value = "white") +
  geom_raster(aes(fill = R0/max(R0))) +
  stat_contour() +
  xlab("Parasite Aggressiveness") + 
  ylab("Parasite Tuning") + 
  guides(fill = guide_colorbar(title = "Relative R0")) +
  scale_x_continuous(
      breaks = seq(1, length(alpha0), by = 10)
    , labels = alpha0[seq(1, length(alpha0), by = 10)]) +
  scale_y_continuous(
      breaks = seq(1, length(tuning), by = 10)
    , labels = plogis(tuning)[seq(1, length(tuning), by = 10)]) +
  theme(legend.key.size = unit(.45, "cm")) +
  stat_contour(
     data = mat_out
    , bins = 5
    , color = "white"
    , alpha = 0.80)
  
}

Sys.sleep(2)
ggsave(filename = paste(paste("gg", i, sep = "_"), ".pdf"))

}
  
}

############
## Summary of a single deterministic run
############

summary.out <- data.frame(
  mean_beta  = rep(0, nrow(res_1000))
, sd_beta    = rep(0, nrow(res_1000))
, mean_R0    = rep(0, nrow(res_1000))
, sd_R0      = rep(0, nrow(res_1000))
, mean_alpha = rep(0, nrow(res_1000))
, sd_alpha   = rep(0, nrow(res_1000))
, mean_tune  = rep(0, nrow(res_1000))
, sd_tune    = rep(0, nrow(res_1000)))

for (i in seq(1, nrow(res_1000), by = 1)) {
  mat_out   <- matrix(data = c(unlist(res_1000[i, -c(1,2)])), nrow = length(tuning), ncol = length(alpha0))
  
  summary.out[i, ] <- transform(
    summary.out[i, ]
    , mean_beta  = sum(mat_out * beta_surface)/sum(mat_out)
    , sd_beta    = sd(mat_out * beta_surface)
    , mean_R0    = sum(mat_out * R0_surface)/sum(mat_out)
    , sd_R0      = sd(mat_out * R0_surface)
    , mean_alpha = sum(sweep(mat_out, 1, alpha0, FUN = "*"))/sum(mat_out)
    , sd_alpha   = sd(sweep(mat_out , 1, alpha0, FUN = "*"))
    , mean_tune  = sum(sweep(mat_out, 2, tuning, FUN = "*"))/sum(mat_out)
    , sd_tune    = sd(sweep(mat_out , 2, tuning, FUN = "*")))
  
}

summary.out <- transform(summary.out
  , Time = seq(1, nrow(summary.out))
  , R0   = rep(0, nrow(summary.out)))

############
## Plot the movement of the deterministic model on the R0 surface
############

if (on_R0_surface == FALSE) {

ggplot(R0_surface_p, aes(Var1, Var2, z = R0/max(R0))) +
  scale_fill_distiller(palette = "Spectral", na.value = "white") +
  geom_raster(aes(fill = R0/max(R0))) +
  stat_contour() +
  xlab("Parasite Aggressiveness") + 
  ylab("Parasite Tuning") + 
  guides(fill = guide_colorbar(title = "Relative R0")) +
  scale_x_continuous(
      breaks = seq(1, length(alpha0), by = 10)
    , labels = alpha0[seq(1, length(alpha0), by = 10)]) +
  scale_y_continuous(
      breaks = seq(1, length(tuning), by = 10)
    , labels = tuning[seq(1, length(tuning), by = 10)]) +
  theme(legend.key.size = unit(.45, "cm")) +
    geom_point(
    data = summary.out[seq(1, nrow(summary.out), by = 3), ]
    , aes(mean_alpha * 100, mean_tune * 100
     , size = sd_alpha + sd_tune
#    , shape = as.factor(mut_mean)
#    , colour = as.factor(tol0)
      )
      , shape = 1
  #  , lwd  = 2.5
      )#  +
 # geom_vline(xintercept = 51, lwd = 0.5, linetype = "dashed")
  
} else {
  
  mat_out_all_temp <- droplevels(subset(mat_out_all
   , Time == 20 | Time == 40 | Time == 100 | Time == 200 | Time == 500))
  
  names(mat_out_all_temp)[3] <- "R0"
  
ggplot(R0_surface_p, aes(Var1, Var2, z = R0/max(R0))) +
  scale_fill_distiller(palette = "Spectral", na.value = "white") +
  geom_raster(aes(fill = R0/max(R0))) +
  stat_contour() +
  xlab("Parasite Aggressiveness") + 
  ylab("Parasite Tuning") + 
  guides(fill = guide_colorbar(title = "Relative R0")) +
  scale_x_continuous(
      breaks = seq(1, length(alpha0), by = 10)
    , labels = alpha0[seq(1, length(alpha0), by = 10)]) +
  scale_y_continuous(
      breaks = seq(1, length(tuning), by = 10)
    , labels = tuning[seq(1, length(tuning), by = 10)]) +
  theme(legend.key.size = unit(.45, "cm")) +
  stat_contour(
     data = mat_out_all_temp
    , aes(group = Time)
    , bins = 5
    , color = "white"
    , alpha = 1.00)
#+ geom_vline(xintercept = 21, lwd = 0.5, linetype = "dashed")  
  
}

############
## gganimate
############

full_contour <- TRUE

if (full_contour == FALSE) {

ggplot(R0_surface_p, aes(Var1, Var2, z = R0/max(R0))) +
  scale_fill_distiller(palette = "Spectral", na.value = "white") +
  geom_raster(aes(fill = R0/max(R0))) +
  stat_contour() +
  xlab("Parasite Aggressiveness") + 
  ylab("Parasite Tuning") + 
  guides(fill = guide_colorbar(title = "Relative R0")) +
  scale_x_continuous(
      breaks = seq(1, length(alpha0), by = 10)
    , labels = alpha0[seq(1, length(alpha0), by = 10)]) +
  scale_y_continuous(
      breaks = seq(1, length(tuning), by = 10)
    , labels = plogis(tuning)[seq(1, length(tuning), by = 10)]) +
  theme(legend.key.size = unit(.45, "cm")) +
    geom_point(
    data = summary.out[seq(1, nrow(summary.out), by = 5), ]
    , aes(mean_alpha * 100, plogis(mean_tune) * 100
     , size = sd_alpha + sd_tune
#    , shape = as.factor(mut_mean)
#    , colour = as.factor(tol0)
      )
      , shape = 1
  #  , lwd  = 2.5
      ) +
  transition_states(
    Time,
    transition_length = 1,
    state_length = 1
   )
  
} else {
  
  names(mat_out_all)[3] <- "R0"
  
  mat_out_all_temp <- droplevels(subset(mat_out_all
   , Time == 10 | Time == 40 | Time == 100 | Time == 200 | Time == 300 | Time == 400 | Time == 500))
  
  ggplot(R0_surface_p, aes(Var1, Var2, z = R0/max(R0))) +
  scale_fill_distiller(palette = "Spectral", na.value = "white") +
  geom_raster(aes(fill = R0/max(R0))) +
  stat_contour() +
  xlab("Parasite Aggressiveness") + 
  ylab("Parasite Tuning") + 
  guides(fill = guide_colorbar(title = "Relative R0")) +
  scale_x_continuous(
      breaks = seq(1, length(alpha0), by = 10)
    , labels = alpha0[seq(1, length(alpha0), by = 10)]) +
  scale_y_continuous(
      breaks = seq(1, length(tuning), by = 10)
    , labels = plogis(tuning)[seq(1, length(tuning), by = 10)]) +
  theme(legend.key.size = unit(.45, "cm")) +
  stat_contour(
     data = mat_out_all_temp
    , bins = 5
    , color = "white"
    , alpha = 0.50) +
    transition_states(
    Time,
    transition_length = 1,
    state_length = 1
   )  
  
  
}

}

##########################################################################################
## Both ##################################################################################
##########################################################################################

############
## Plot both the deterministic and the stochastic models results on the same contour plot
############ 

  mat_out_all_temp <- droplevels(subset(mat_out_all
   , Time == 20 |
      Time == 85 |
      Time == 90 |
      Time == 200 | 
      Time == 600))
  
  names(mat_out_all_temp)[3] <- "R0"

ggplot(R0_surface_p, aes(Var1, Var2, z = R0/max(R0))) +
  scale_fill_distiller(palette = "Spectral", na.value = "white") +
  geom_raster(aes(fill = R0/max(R0))) +
  stat_contour() +
  xlab("Parasite Aggressiveness") + 
  ylab("Parasite Tuning") + 
  guides(fill = guide_colorbar(title = "Relative R0")) +
  scale_x_continuous(
      breaks = seq(1, length(alpha0), by = 10)
    , labels = alpha0[seq(1, length(alpha0), by = 10)]) +
  scale_y_continuous(
      breaks = seq(1, length(tuning), by = 10)
    , labels = tuning[seq(1, length(tuning), by = 10)]) +
  theme(legend.key.size = unit(.45, "cm")) +
  stat_contour(
     data = mat_out_all_temp
    , aes(group = Time)
    , bins = 5
    , color = "white"
    , alpha = 1.00) +
  geom_vline(xintercept = 48.5, lwd = 0.5, linetype = "dashed") +
  geom_point(
    data = res_1000_stochas[
      c(
        seq(1, nrow(res_1000_stochas) / 6, by = 1)
      , seq(nrow(res_1000_stochas) / 6, nrow(res_1000_stochas) / 2, by = 20)
      , seq(nrow(res_1000_stochas) / 2, nrow(res_1000_stochas), by = 40))
      , ]
    , aes(plogis(mean_plalpha) * 100, plogis(mean_plbeta) * 100, size = sd_plalpha + sd_plbeta)
    , shape = 1
  # , lwd  = 1.0
    ) +
  geom_path(
     data = res_1000_stochas[seq(1, nrow(res_1000_stochas), by = 1), ]
    , aes(plogis(mean_plalpha) * 100, plogis(mean_plbeta) * 100)
    , lwd  = 0.5) +
  geom_path(
     data = res_1000_stochas2[seq(1, nrow(res_1000_stochas2), by = 1), ]
    , aes(plogis(mean_plalpha) * 100, plogis(mean_plbeta) * 100)
    , lwd  = 0.5) +
  geom_path(
     data = res_1000_stochas3[seq(1, nrow(res_1000_stochas3), by = 1), ]
    , aes(plogis(mean_plalpha) * 100, plogis(mean_plbeta) * 100)
    , lwd  = 0.5) +
  geom_path(
     data = res_1000_stochas4[seq(1, nrow(res_1000_stochas4), by = 1), ]
    , aes(plogis(mean_plalpha) * 100, plogis(mean_plbeta) * 100)
    , lwd  = 0.5) 


############
## Plot all stochastic runs
############ 

mat_out_all_temp <- droplevels(subset(mat_out_all
   , Time == 20 |
      Time == 50 |
      Time == 60 |
      Time == 90 |
      Time == 200 |
      Time == 400 |      
      Time == 800))
  
  names(mat_out_all_temp)[3] <- "R0"
  
res_1000_stochas_s <- transform(res_1000_stochas_s, R0 = 0)  

ggplot(R0_surface_p, aes(Var1, Var2, z = R0/max(R0))) +
  scale_fill_distiller(palette = "Spectral", na.value = "white") +
  geom_raster(aes(fill = R0/max(R0))) +
  stat_contour() +
  xlab("Parasite Aggressiveness") + 
  ylab("Parasite Tuning") + 
  guides(fill = guide_colorbar(title = "Relative R0")) +
  scale_x_continuous(
      breaks = seq(1, length(alpha0), by = 10)
    , labels = alpha0[seq(1, length(alpha0), by = 10)]) +
  scale_y_continuous(
      breaks = seq(1, length(tuning), by = 10)
    , labels = tuning[seq(1, length(tuning), by = 10)]) +
  theme(legend.key.size = unit(.45, "cm")) +
  stat_contour(
     data = mat_out_all_temp
    , aes(group = Time)
    , bins = 5
    , color = "white"
    , alpha = 1.00) +
  geom_vline(xintercept = 21, lwd = 0.5, linetype = "dashed") +
  geom_path(
     data = res_1000_stochas[seq(1, nrow(res_1000_stochas), by = 1), ]
    , aes(plogis(mean_plalpha) * 100, plogis(mean_plbeta) * 100
      , group = run)
    , lwd  = 0.2, alpha = 0.2)

############
## Plot a summary of multiple runs
############ 

mat_out_all_temp <- droplevels(subset(mat_out_all
   , Time == 20 |
      Time == 50 |
      Time == 60 |
      Time == 90 |
      Time == 110 |
      Time == 200 |
      Time == 400 |      
      Time == 800))
  
  names(mat_out_all_temp)[3] <- "R0"
  
res_1000_stochas_s <- transform(res_1000_stochas_s, R0 = 0)  

ggplot(R0_surface_p, aes(Var1, Var2, z = R0/max(R0))) +
  scale_fill_distiller(palette = "Dark2", na.value = "white") +
  geom_raster(aes(fill = R0/max(R0))) +
  stat_contour() +
  xlab("Parasite Aggressiveness") + 
  ylab("Parasite Tuning") + 
  guides(fill = guide_colorbar(title = "Relative R0")) +
  scale_x_continuous(
      breaks = seq(1, length(alpha0), by = 10)
    , labels = alpha0[seq(1, length(alpha0), by = 10)]) +
  scale_y_continuous(
      breaks = seq(1, length(tuning), by = 10)
    , labels = tuning[seq(1, length(tuning), by = 10)]) +
  theme(legend.key.size = unit(.45, "cm")) +
  geom_vline(xintercept = 21, lwd = 0.5, linetype = "dashed") +
  geom_path(
     data = res_1000_stochas[seq(1, nrow(res_1000_stochas), by = 1), ]
    , aes(plogis(mean_plalpha) * 100, plogis(mean_plbeta) * 100
      , group = run
      )
    , lwd  = 0.2, alpha = 0.4) +
  geom_path(
     data = res_1000_stochas_s[seq(1, nrow(res_1000_stochas_s), by = 1), ]
    , aes(plogis(q05_mean_plalpha) * 100, plogis(q05_mean_plbeta) * 100)
    , lwd  = 0.4, colour = "white") +
  geom_path(
     data = res_1000_stochas_s[seq(1, nrow(res_1000_stochas_s), by = 1), ]
    , aes(plogis(q25_mean_plalpha) * 100, plogis(q25_mean_plbeta) * 100)
    , lwd  = 0.7, colour = "white") +
  geom_path(
     data = res_1000_stochas_s[seq(1, nrow(res_1000_stochas_s), by = 1), ]
    , aes(plogis(q50_mean_plalpha) * 100, plogis(q50_mean_plbeta) * 100)
    , lwd  = 1.1, colour = "white") +
  geom_path(
     data = res_1000_stochas_s[seq(1, nrow(res_1000_stochas_s), by = 1), ]
    , aes(plogis(q75_mean_plalpha) * 100, plogis(q75_mean_plbeta) * 100)
    , lwd  = 0.7, colour = "white") +
  geom_path(
     data = res_1000_stochas_s[seq(1, nrow(res_1000_stochas_s), by = 1), ]
    , aes(plogis(q95_mean_plalpha) * 100, plogis(q95_mean_plbeta) * 100)
    , lwd  = 0.4, colour = "white") +
  stat_contour(
     data = mat_out_all_temp
    , aes(group = Time)
    , bins = 5
    , color = "white"
    , alpha = 1.00)


##########################################################################################
## Both, also including the AD solution###################################################
##########################################################################################

mat_out_all_temp <- droplevels(subset(mat_out_all
  ,   Time == 10 |
  #    Time == 20 |
  #    Time == 28 |
  #    Time == 32 |
      Time == 46 |
  #    Time == 60 |
#      Time == 70 |
#      Time == 74 |
      Time == 100 |
#      Time == 200 |
     Time == 260 |
#      Time == 300 |
      Time == 400
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
  geom_path(
     data = res_1000_stochas[seq(1, nrow(res_1000_stochas), by = 1), ]
#    , aes(plogis(mean_plalpha) * 100, plogis(mean_plbeta) * 100
     , aes(mean_plalpha, mean_plbeta
      , group = run
       )
    , lwd  = 0.2, alpha = 0.3) +
  geom_path(
     data = res_1000_stochas_s[seq(1, nrow(res_1000_stochas_s), by = 1), ]
  #  , aes(plogis(q05_mean_plalpha) * 100, plogis(q05_mean_plbeta) * 100)
   , aes(q05_mean_plalpha, q05_mean_plbeta)
    , lwd  = 0.4, colour = "gray80") +
  geom_path(
     data = res_1000_stochas_s[seq(1, nrow(res_1000_stochas_s), by = 1), ]
 #   , aes(plogis(q25_mean_plalpha) * 100, plogis(q25_mean_plbeta) * 100)
    , aes(q25_mean_plalpha, q25_mean_plbeta)
    , lwd  = 0.7, colour = "gray80") +
  geom_path(
     data = res_1000_stochas_s[seq(1, nrow(res_1000_stochas_s), by = 1), ]
#    , aes(plogis(q50_mean_plalpha) * 100, plogis(q50_mean_plbeta) * 100)
    , aes(q50_mean_plalpha, q50_mean_plbeta)
    , lwd  = 1.1, colour = "gray80") +
  geom_path(
     data = res_1000_stochas_s[seq(1, nrow(res_1000_stochas_s), by = 1), ]
 #   , aes(plogis(q75_mean_plalpha) * 100, plogis(q75_mean_plbeta) * 100)
    , aes(q75_mean_plalpha, q75_mean_plbeta)
    , lwd  = 0.7, colour = "gray80") +
  geom_path(
     data = res_1000_stochas_s[seq(1, nrow(res_1000_stochas_s), by = 1), ]
#    , aes(plogis(q95_mean_plalpha) * 100, plogis(q95_mean_plbeta) * 100)
    , aes(q95_mean_plalpha, q95_mean_plbeta)
    , lwd  = 0.4, colour = "gray80") +
  stat_contour(
     data = mat_out_all_temp
    , aes(group = Time)
    , bins = 5
    , color = "white"
    , alpha = 1.00) +
  geom_path(
     data = grad_ascent[[1]]
    , aes(alpha, tune)
    , lwd  = 1.0, colour = "white") #+
 # geom_path(
#     data = grad_ascentF[[1]]
#    , aes(alpha , tune)
#    , lwd  = 2.0, colour = "white", linetype = "dotted")
#  geom_path(
#     data = grad_max_range[grad_max_range$priority_move == "Tune", ]
#    , aes(alpha_loc / 10, tune_loc / 10)
#    , lwd  = 1.5, colour = "white", linetype = "dashed") +
#  geom_path(
#     data = grad_max_range[grad_max_range$priority_move == "Agg", ]
#    , aes(alpha_loc / 10, tune_loc / 10)
#    , lwd  = 1.5, colour = "white", linetype = "dashed")  

##########################################################################################
## Compare just R0 of AD and reaction diffusion model ####################################
##########################################################################################

## Timescale of these two are hard because of how they are set up (diffusion is likely
 ## comparable in some way to a mutational step, but it isn't immediately clear how), 
  ## so force them to line up to simply check the curvature to the peak

grad_ascent   <- transform(grad_ascent, Time = seq(1, nrow(grad_ascent) / 5, length = nrow(grad_ascent)))
determ_med_R0 <- transform(determ_med_R0, time = seq(1, nrow(determ_med_R0) * 5, length = nrow(determ_med_R0)))

ggplot(determ_med_R0, aes(time, med_R0)) + 
  geom_line(lwd = 1) +
#  geom_line(aes(Time, med_R0, lwd = 0.5)) +
  geom_line(data = grad_ascent
    , aes(time, R0), lwd = 1, linetype = "dotted", col = "blue")

mat_out_all_temp <- droplevels(subset(mat_out_all, Time == 38))

ggplot(mat_out_all_temp, aes(Var1, Var2, z = value)) +
  scale_fill_distiller(palette = "Dark2", na.value = "white") +
  geom_raster(aes(fill = value)) +
  stat_contour() +
  xlab("Parasite Aggressiveness") + 
  ylab("Parasite Tuning")

########################################################################################
## Check the gradient at each step for the AD model ####################################
########################################################################################

grad_ascent[[2]] <- transform(grad_ascent[[2]], step = seq(1, nrow(grad_ascent[[2]]), by = 1))

grad_ascent_check <- melt(grad_ascent[[2]], "step")
grad_ascent_check <- grad_ascent_check %>% filter(step != 1)

ggplot(grad_ascent_check, aes(step, value)) + geom_line(aes(colour = variable))

grad_ascent_check2 <- grad_ascent_check %>% filter(step > 460 & variable != "down" & variable != "diagDR")

ggplot(grad_ascent_check2, aes(step, value)) + geom_line(aes(colour = variable))

