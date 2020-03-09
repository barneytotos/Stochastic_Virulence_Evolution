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

endsteps     <- 50
model.choice <- "nt"

### load previous results
res_all  <- readRDS("batch_runs/nt_FALSE_500_Sun_Mar_08_2020.Rds")

plot_runs <- c(4, 7, 8)

### Can first explore just single runs to check to see if the model is working. What to check will vary by the model chosen
gg_check <- ggplot((res_all %>% filter(param_num <= 12, time < 20000))
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
  , mean_negtrait       = mean(mean_negtrait)
  , sd_negtrait         = mean(sd_negtrait)
  , shannon.div         = mean(shannon)
  , mu                  = mean(mu)
  , mut_mean            = mean(mut_mean)
  , mut_sd              = mean(mut_sd)
  , N                   = mean(N)
  , gamma0              = mean(gamma0)
    )

res_all_s.a.res <- res_all_s.a %>% dplyr::select(
  param_num
, mean_postrait
, sd_postrait
, mean_negtrait
, sd_negtrait
, shannon.div
  )  

if (model.choice != "nt") {

res_all_s.a.par <- res_all_s.a %>% dplyr::select(
  param_num
, mu
, mut_mean
, mut_sd
, N
, gamma0
  )

} else {
  
res_all_s.a <- left_join(res_all_s.a
  , (res_all %>% dplyr::select(param_num, nt_mut_var_pos_trait) %>% filter(!duplicated(param_num)))
  , by = "param_num")
  
res_all_s.a.par <- res_all_s.a %>% dplyr::select(
  param_num
, mu
, mut_mean
, mut_sd
, N
, gamma0
, nt_mut_var_pos_trait
  )  
  
}

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

if (model.choice != "nt") {

## Plot the hypercube results
ggplot(res_all_s.gg, aes(Param.Value, Out.Value)) +
    geom_point(aes(colour = mu)) +
    facet_grid(Out.Name ~ Param.Name, scale = "free") +
    scale_color_viridis() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, size = 14)
      )

} else {
  
## probably could do this all in one dplyr call but w/e
pos_only <- res_all_s.a.par[with(res_all_s.a.par, which(nt_mut_var_pos_trait == TRUE)), ]$param_num
neg_only <- res_all_s.a.par[with(res_all_s.a.par, which(nt_mut_var_pos_trait == FALSE)), ]$param_num
  
ggplot(
  (res_all_s.gg %>% filter(param_num %in% pos_only))
  , aes(Param.Value, Out.Value)) +
    geom_point(aes(colour = mu)) +
    facet_grid(Out.Name ~ Param.Name, scale = "free") +
    scale_color_viridis() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, size = 14)
      )
  
ggplot(
  (res_all_s.gg %>% filter(param_num %in% neg_only))
  , aes(Param.Value, Out.Value)) +
    geom_point(aes(colour = mu)) +
    facet_grid(Out.Name ~ Param.Name, scale = "free") +
    scale_color_viridis() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, size = 14)
      )
  
}

