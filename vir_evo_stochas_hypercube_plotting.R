#########################################################
### Some initial plotting of these hypercube results. ###
#########################################################

 ## Remember we want:
## 1) Time to equilibrium. Treat this as first passage time. 
 ## First need to add equilibrium to the output. Added to params
 ## !*! May want some sort of eps to check if the mean has come within?
## 2) Transient movement when approaching equilibrium
 ## Not quite sure about this one but for plotting purposes treat as sd prior to equilibrium for now
## 3) SD around equilibrium
 ## Mean SD after first passage
## 4) ... 

## Jump through some hoops to determine when the chains reached equilibrium
res_1000_all_s <- res_1000_all %>%
  mutate(
    higher_start = ifelse(alpha0 > opt_alpha, 1, 0)
  , lower_start  = ifelse(alpha0 < opt_alpha, 1, 0)
  , higher_equil = ifelse(mean_alpha > opt_alpha, 1, 0)
  , lower_equil  = ifelse(mean_alpha < opt_alpha, 1, 0)) %>%
  mutate(
  first_pass_setup = ifelse(lower_start > 1
      , ifelse(mean_alpha > opt_alpha, 1, 0)
      , ifelse(mean_alpha < opt_alpha, 1, 0))
  ) %>% 
  group_by(param_num) %>%
  mutate(
    time_point = row_number()
  ) %>%
  mutate(
    when_equil = min(which(first_pass_setup == 1))
  ) %>% 
  mutate(
    at_equil = ifelse(time_point >= when_equil, 1, 0)
  ) %>% 
  mutate(
    cum_time = cumsum(rptfreq)
  )

## Then calc sd prior and after equilibirum. Could also consider looking at some sort of
 ## transient max or min or something
res_1000_all_s.a <- res_1000_all_s %>%
  filter(at_equil == 1) %>%
  group_by(param_num) %>%
  summarize(
    first_equil    = mean(when_equil)
  , alpha.sd_after = mean(sd_alpha)
  , alpha.dr.after = sd(mean_alpha)
  , beta.sd_after  = mean(sd_beta)
  , beta.dr.after  = sd(mean_beta)
  , time_to_equil  = min(cum_time))

res_1000_all_s.b <- res_1000_all_s %>%
  filter(at_equil == 0) %>%
  group_by(param_num) %>%
  summarize(
    alpha.sd_before = mean(sd_alpha)
  , alpha.dr.before = sd(mean_alpha)
  , beta.sd_before = mean(sd_beta)
  , beta.dr.before  = sd(mean_beta))

## Melt first to get the variables of interest with names and values and then add back parameter values
res_1000_all_s.gg <- left_join(res_1000_all_s.a, res_1000_all_s.b, "param_num")
## Have time to equilibrium so dont need first time sample (#) of equilibrium
res_1000_all_s.gg <- res_1000_all_s.gg[, -grep("first_equil", names(res_1000_all_s.gg))]
## take log of first_equil then melt
res_1000_all_s.gg <- transform(res_1000_all_s.gg, time_to_equil = log(time_to_equil))
res_1000_all_s.gg <- melt(res_1000_all_s.gg, "param_num")

## melt the parameter values as well
if (tradeoff_only == TRUE) {
params.s    <- params %>% dplyr::select(param_num, mu, mut_sd, alpha0, power_c, power_exp, N, gamma0)
} else {
params.s    <- params %>% dplyr::select(param_num, mu, mut_sd, alpha0, power_c, power_exp, N, gamma0, eff_hit)  
}
params.melt <- melt(params.s, "param_num")

## Add pack the parameter values for plotting
res_1000_all_s.gg        <- left_join(res_1000_all_s.gg, params.melt, by = "param_num")
names(res_1000_all_s.gg) <- c("param_num", "Out.Name", "Out.Value", "Param.Name", "Param.Value")
res_1000_all_s.gg        <- as.data.frame(res_1000_all_s.gg)

## add in a column (one parameters) to color by 
res_1000_all_s.gg <- left_join(res_1000_all_s.gg, params.s[, c(1, 2)], "param_num")

## remove the rows with NA where the chains get to equil before the first time point
res_1000_all_s.gg <- res_1000_all_s.gg[complete.cases(res_1000_all_s.gg), ]

library(viridis)
## gg pairs plotting of the hypercube results
res_1000_all_s.gg2 <- res_1000_all_s.gg %>% 
  filter(Param.Name != "alpha0")

after.rows  <- grep("after", res_1000_all_s.gg2$Out.Name)
before.rows <- grep("before", res_1000_all_s.gg2$Out.Name)

res_1000_all_s.gg2.a <- res_1000_all_s.gg2[after.rows, ]
res_1000_all_s.gg2.a <- droplevels(res_1000_all_s.gg2.a)

res_1000_all_s.gg2.b <- res_1000_all_s.gg2[before.rows, ]
res_1000_all_s.gg2.b <- droplevels(res_1000_all_s.gg2.b)

ggplot(res_1000_all_s.gg2.a
  , aes(Param.Value, Out.Value)) +
    geom_point(aes(colour = mu)) +
    facet_grid(Out.Name ~ Param.Name, scale = "free") +
    scale_color_viridis() 

ggplot(res_1000_all_s.gg2.b
  , aes(Param.Value, Out.Value)) +
    geom_point(aes(colour = mu)) +
    facet_grid(Out.Name ~ Param.Name, scale = "free") +
    scale_color_viridis() 
