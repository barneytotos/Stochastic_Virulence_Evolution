##############################################
### Plotting for new batch sims stochastic ###
##############################################

## Remember we want:

## 1) Time to equilibrium. Treat this as first passage time. 
 ## First need to add equilibrium to the output. Added to params
 ## !*! May want some sort of eps to check if the mean has come within?
## 2) Transient movement when approaching equilibrium
 ## Not quite sure about this one but for plotting purposes treat as sd prior to equilibrium for now
## 3) SD around equilibrium
 ## Mean SD after first passage
## 4) ... 

endsteps <- 50

### load previous results
res_all  <- readRDS("batch_runs/eff_FALSE_100_Fri_Mar_06_2020.Rds")

plot_runs <- c(4, 7, 8)

### Can first explore just single runs to check to see if the model is working. What to check will vary by the model chosen
gg_check <- ggplot((res_all %>% filter(param_num <= 18, time < 20000))
  , aes(time, mean_postrait)) + geom_path() + 
  geom_path(aes(time, mean_negtrait), col = "blue", lwd = 1) +
  facet_wrap(~param_num)
gg_check

## Record the time point number
res_all_s <- res_all %>%
  group_by(param_num) %>%
  mutate(
    time_point = row_number()
  , cum_time   = cumsum(rptfreq)
  )

## Determine where the starting values were relative to the optimum. Some of these will/will not be relevant 
 ## depending on what traits will be evolving. Subset by those later, but calculate all up front
  ## Will throw a few warnings for model.choice == "nt" but w/e
res_all_s <- res_all_s %>%
  mutate(
    higher_negtrait_start = ifelse(neg_trait0 > opt_negtrait, 1, 0)
  , lower_negtrait_start  = ifelse(neg_trait0 < opt_negtrait, 1, 0)
  , higher_negtrait_equil = ifelse(mean_negtrait > neg_trait0, 1, 0)
  , lower_negtrait_equil  = ifelse(mean_negtrait < neg_trait0, 1, 0)
  , higher_postrait_start = ifelse(pos_trait0 > opt_postrait, 1, 0)
  , lower_postrait_start  = ifelse(pos_trait0 < opt_postrait, 1, 0)
  , higher_postrait_equil = ifelse(mean_postrait > pos_trait0, 1, 0)
  , lower_postrait_equil  = ifelse(mean_postrait < pos_trait0, 1, 0)    
    ) %>%
  mutate(
  first_pass_setup_negtrait = ifelse(lower_negtrait_start == 1
      , ifelse(mean_negtrait > opt_negtrait, 1, 0)
      , ifelse(mean_negtrait < opt_negtrait, 1, 0))
, first_pass_setup_postrait = ifelse(lower_postrait_start == 1
      , ifelse(mean_postrait > opt_postrait, 1, 0)
      , ifelse(mean_postrait < opt_postrait, 1, 0))
  ) %>% 
  group_by(param_num) %>%
  mutate(
    time_point = row_number()
  ) %>%
  mutate(
    when_equil_negtrait = min(which(first_pass_setup_negtrait == 1))
  , when_equil_postrait = min(which(first_pass_setup_postrait == 1))
  ) %>% 
  mutate(
    at_equil_negtrait = ifelse(time_point >= when_equil_negtrait, 1, 0)
  , at_equil_postrait = ifelse(time_point >= when_equil_postrait, 1, 0)
  ) %>% 
  mutate(
    cum_time = cumsum(rptfreq)
  )

## Subset to the last X number of time points for summary
res_all_s.a <- res_all_s %>%
  filter(
    time_point >= (max(time_point) - endsteps)
  )

res_all_s.a <- res_all_s.a %>%
  summarize(
    mean_postrait       = mean(mean_postrait)
  , sd_postrait         = mean(sd_postrait)
  , shannon.div         = mean(shannon)
  , mu                  = mean(mu)
  , mut_mean            = mean(mut_mean)
  , mut_sd              = mean(mut_sd)
  , N                   = mean(N)
  , gamma0              = mean(gamma0)
    )

res_all_s.a.res <- res_all_s.a %>% dplyr::select(param_num, mean_postrait, sd_postrait, shannon.div)  
res_all_s.a.par <- res_all_s.a %>% dplyr::select(param_num, mu, mut_mean, mut_sd, N, gamma0)

## Melt
res_all_s.a.res.gg <- melt(res_all_s.a.res, c("param_num"))
res_all_s.a.par.gg <- melt(res_all_s.a.par, c("param_num"))

## Combine 
res_all_s.gg       <- left_join(
    res_all_s.a.res.gg
  , res_all_s.a.par.gg
  , by = "param_num")

names(res_all_s.gg) <- c("param_num", "Out.Name", "Out.Value", "Param.Name", "Param.Value")

## add in a column (one parameters) to color by 
param.f_col       <- res_all_s.a[, c(1, grep("mu", colnames(res_all_s.a))[1])]
param.f_col       <- param.f_col[!duplicated(param.f_col), ]

res_all_s.gg <- left_join(res_all_s.gg, param.f_col, "param_num")

## Plot the hypercube results
ggplot(res_all_s.gg, aes(Param.Value, Out.Value)) +
    geom_point(aes(colour = mu)) +
    facet_grid(Out.Name ~ Param.Name, scale = "free") +
    scale_color_viridis() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, size = 14)
      )


