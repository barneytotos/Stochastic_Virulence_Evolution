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

### load previous results
## res_1000_all <- readRDS("Saved_Results/res_1000_all_stochas_eff_0_Tue_Jan_21_2020.rds")
## res_1000_all1 <- readRDS("Saved_Results/res_1000_all_stochas_eff_300_Thu_Jan_23_2020.rds")
## res_1000_all <- rbind(res_1000_all, res_1000_all1)

## Jump through some hoops to determine when the chains reached equilibrium
res_1000_all_s <- res_1000_all %>%
  mutate(
    higher_start = ifelse(alpha0 > opt_alpha, 1, 0)
  , lower_start  = ifelse(alpha0 < opt_alpha, 1, 0)
  , higher_equil = ifelse(mean_alpha > opt_alpha, 1, 0)
  , lower_equil  = ifelse(mean_alpha < opt_alpha, 1, 0)) %>%
  mutate(
  first_pass_setup = ifelse(lower_start == 1
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

## melt the parameter values as well
if (tradeoff_only) {
params.s    <- params %>% dplyr::select(param_num, mu, mut_sd, alpha0, power_c, power_exp, N, gamma0)
} else {
# params.s    <- params %>% dplyr::select(param_num, mu, mut_sd, alpha0, power_c, power_exp, N, gamma0, eff_hit)  
## Not sure why these are needed, as all of the params are also stored in the sim output data frame
params.s    <-  res_1000_all_s %>% dplyr::select(
  param_num
, mu
, mut_sd
, alpha0
, power_c
, power_exp
, N
, gamma0
# , eff_hit
, opt_alpha
) 
}

if (!exists("params.melt")) {
 params.melt <- melt(params.s, "param_num")
 params.melt <- params.melt[!duplicated(params.melt), ]
}

## Then calc sd prior and after equilibirum. Could also consider looking at some sort of
 ## transient max or min or something
res_1000_all_s.a <- res_1000_all_s %>%
  filter(at_equil == 1) %>%
  group_by(param_num) %>%
  summarize(
    first_equil         = mean(when_equil)
  , alpha.sd_after      = mean(sd_alpha)
  , alpha.dr.after      = sd(mean_alpha)
  , beta.sd_after       = mean(sd_beta)
  , beta.dr.after       = sd(mean_beta)
  , time_to_equil.after = min(cum_time)
  , max.alpha.after     = max(mean_alpha)
  , shannon.div.after   = mean(shannon))

res_1000_all_s.b <- res_1000_all_s %>%
  filter(at_equil == 0) %>%
  group_by(param_num) %>%
  summarize(
    alpha.sd_before   = mean(sd_alpha)
  , alpha.dr.before   = sd(mean_alpha)
  , beta.sd_before    = mean(sd_beta)
  , beta.dr.before    = sd(mean_beta)
  , max.alpha.before   = max(mean_alpha)
  , shannon.div.before = mean(shannon))

## Melt first to get the variables of interest with names and values and then add back parameter values
res_1000_all_s.gg <- left_join(res_1000_all_s.a, res_1000_all_s.b, "param_num")
## Have time to equilibrium so dont need first time sample (#) of equilibrium
res_1000_all_s.gg <- res_1000_all_s.gg[, -grep("first_equil", names(res_1000_all_s.gg))]
## take log of first_equil then melt
res_1000_all_s.gg <- transform(res_1000_all_s.gg, time_to_equil = log(time_to_equil.after))
res_1000_all_s.gg <- melt(res_1000_all_s.gg, "param_num")

## Add pack the parameter values for plotting
res_1000_all_s.gg        <- left_join(res_1000_all_s.gg, params.melt, by = "param_num")
names(res_1000_all_s.gg) <- c("param_num", "Out.Name", "Out.Value", "Param.Name", "Param.Value")
res_1000_all_s.gg        <- as.data.frame(res_1000_all_s.gg)

## add in a column (one parameters) to color by 
param.f_col       <- params.s[, c(1, grep("power_c", colnames(params.s))[1])]
param.f_col       <- param.f_col[!duplicated(param.f_col), ]

res_1000_all_s.gg <- left_join(res_1000_all_s.gg, param.f_col, "param_num")

## remove the rows with NA where the chains get to equil before the first time point
res_1000_all_s.gg <- res_1000_all_s.gg[complete.cases(res_1000_all_s.gg), ]

## gg pairs plotting of the hypercube results
res_1000_all_s.gg2 <- res_1000_all_s.gg %>% 
  filter(Param.Name != "alpha0")

after.rows  <- grep("after", res_1000_all_s.gg2$Out.Name)
before.rows <- grep("before", res_1000_all_s.gg2$Out.Name)

res_1000_all_s.gg2.a <- res_1000_all_s.gg2[after.rows, ]
res_1000_all_s.gg2.a <- droplevels(res_1000_all_s.gg2.a)

res_1000_all_s.gg2.b <- res_1000_all_s.gg2[before.rows, ]
res_1000_all_s.gg2.b <- droplevels(res_1000_all_s.gg2.b)

ggplot(res_1000_all_s.gg2.a[
  res_1000_all_s.gg2.a$Out.Name != "time_to_equil.after", ]
  , aes(Param.Value, Out.Value)) +
    geom_point(aes(
   #   colour = mu
      colour = power_c
      )) +
    facet_grid(Out.Name ~ Param.Name, scale = "free") +
    scale_color_viridis() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(res_1000_all_s.gg2.b
  , aes(Param.Value, Out.Value)) +
    geom_point(aes(colour = mu)) +
    facet_grid(Out.Name ~ Param.Name, scale = "free") +
    scale_color_viridis() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(res_1000_all_s.gg2.a[
  res_1000_all_s.gg2.a$Out.Name == "time_to_equil.after", ]
  , aes(Param.Value, Out.Value)) +
    geom_point(aes(colour = mu)) +
    facet_grid(Out.Name ~ Param.Name, scale = "free") +
    scale_color_viridis() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))

#####
## Plot just a few full runs
#####

ggplot(res_1000_all_s[res_1000_all_s$param_num < 22, ]
  , aes(time, mean_alpha)) +
  geom_path() +
  facet_wrap(~ param_num) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_x_log10()
  
  


