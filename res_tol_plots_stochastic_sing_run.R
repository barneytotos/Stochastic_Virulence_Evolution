############
## Array of plots from a single run of the stochastic model
############

#### ggplots of run
if (params$parasite_tuning[1] == FALSE) {
  gg1.1 <- ggplot(res_1000_all, aes(time / nrpt, plogis(mean_plbeta)))  + geom_line(aes(colour = fill_birth, linetype = mut_mean)) +
    xlab("Time (steps)") + ylab("Parasite Efficiency") 
  gg1.2 <- ggplot(res_1000_all, aes(time / nrpt, power_tradeoff(plogis(mean_plalpha), params$power_c[1], params$power_exp[1]) * 
        plogis(mean_plbeta)))                                             + geom_line(aes(colour = fill_birth, linetype = mut_mean)) +
      xlab("Time (steps)") + ylab("Parasite Realized Beta")
} else {
  gg1.1 <- ggplot(res_1000_all, aes(time / nrpt
    , exp(-(mean_plbeta - mean_plalpha)^2 / eff_scale)))  +
      geom_line(aes(colour = fill_birth, linetype = mut_mean)) +
      xlab("Time (steps)") + ylab("Parasite Efficiency")
  gg1.1b <- ggplot(res_1000_all, aes(time / nrpt, mean_plbeta)) +
      geom_line(aes(colour = fill_birth, linetype = mut_mean)) +
      xlab("Time (steps)") + ylab("Parasite Tuning")
  gg1.2 <- ggplot(res_1000_all, aes(time / nrpt, power_tradeoff(plogis(mean_plalpha), params$power_c[1], params$power_exp[1]) * 
        exp(-(mean_plbeta - mean_plalpha)^2 / eff_scale))) + geom_line(aes(colour = fill_birth, linetype = mut_mean)) +
      xlab("Time (steps)") + ylab("Parasite Realized Beta")
}
    gg1.3 <- ggplot(res_1000_all, aes(time / nrpt, plogis(mean_plalpha))) + geom_line(aes(colour = fill_birth, linetype = mut_mean)) + 
      xlab("Time (steps)") + ylab("Parasite alpha")
if (params$parasite_tuning[1] == FALSE) {
    gg1.4 <- ggplot(res_1000_all, aes(time / nrpt
      , (( power_tradeoff(plogis(mean_plalpha), params$power_c[1], params$power_exp[1]) * plogis(mean_plbeta) ) / 
      ( plogis(mean_plalpha) + 0.2 + params$d[1] ))))                     + geom_line(aes(colour = fill_birth, linetype = mut_mean)) + 
      xlab("Time (steps)") + ylab("Parasite R0")c
} else {
    gg1.4 <- ggplot(res_1000_all, aes(time / nrpt
      , (( power_tradeoff(plogis(mean_plalpha), params$power_c[1], params$power_exp[1]) * 
          exp(-(mean_plbeta - mean_plalpha)^2 / eff_scale)) / 
      ( plogis(mean_plalpha) + 0.2 + params$d[1] ))))                     + geom_line(aes(colour = fill_birth, linetype = mut_mean)) + 
      xlab("Time (steps)") + ylab("Parasite R0")  
}
    
    ### For tracking actual mortality rate and transmission rate in a co-evolving population, just take the 
     ### mean of the matrix of all alpha and the mean of the matrix of all beta???
    
    ## also need to track parasite alpha with host effects
    gg1.2b <- ggplot(res_1000_all, aes(time / nrpt, mean_alpha))     + geom_line(aes(colour = fill_birth, linetype = mut_mean)) + 
      xlab("Time (steps)") + ylab("Overall alpha")
    ## and parasite beta with host effects
    gg1.3b <- ggplot(res_1000_all, aes(time / nrpt, mean_beta))      + geom_line(aes(colour = fill_birth, linetype = mut_mean)) + 
      xlab("Time (steps)") + ylab("Overall beta")     
    ## and parasite R0 in the evolved host population
    gg1.4b <- ggplot(res_1000_all, aes(time / nrpt
      , (mean_beta / ( mean_alpha + 0.2 + params$d[1] ))))           + geom_line(aes(colour = fill_birth, linetype = mut_mean)) + 
      xlab("Time (steps)") + ylab("Overall R0")     
      
    gg1.5 <- ggplot(res_1000_all, aes(time / nrpt, mean_hlres))      + geom_line(aes(colour = fill_birth, linetype = mut_mean)) +
      xlab("Time (steps)") + ylab("Host Resistance")
    gg1.6 <- ggplot(res_1000_all, aes(time / nrpt, mean_hltol))      + geom_line(aes(colour = fill_birth, linetype = mut_mean)) + 
      xlab("Time (steps)") + ylab("Host Tolerance")
  
    gg2.1 <- ggplot(res_1000_all, aes(time / nrpt, num_I_strains))   + geom_line(aes(colour = fill_birth, linetype = mut_mean)) + 
      xlab("Time (steps)") + ylab("Number of I Strains")
    gg2.2 <- ggplot(res_1000_all, aes(time / nrpt, num_S_strains))   + geom_line(aes(colour = fill_birth, linetype = mut_mean)) + 
      xlab("Time (steps)") + ylab("Number of S Strains")
   
    gg3.1 <- ggplot(res_1000_all, aes(time / nrpt, num_S))           + geom_line(aes(colour = fill_birth, linetype = mut_mean)) + 
      xlab("Time (steps)") + ylab("Number of S Individuals")
    gg3.2 <- ggplot(res_1000_all, aes(time / nrpt, num_I))           + geom_line(aes(colour = fill_birth, linetype = mut_mean)) + 
      xlab("Time (steps)") + ylab("Number of I Individuals")
    gg3.3 <- ggplot(res_1000_all, aes(time / nrpt, pop_size))        + geom_line(aes(colour = fill_birth, linetype = mut_mean)) + 
      xlab("Time (steps)") + ylab("Total population Size")
  
  #### ggplots tracking progress on tradeoff space
  power_trade_dat <- data.frame(
    alpha = seq(0.01, 1.0, by = 0.01)
  , beta  = power_tradeoff(alpha = seq(0.01, 1.0, by = 0.01), c = params$power_c[1], curv = params$power_exp[1])
  , R0    = power_tradeoff(alpha = seq(0.01, 1.0, by = 0.01), c = params$power_c[1], curv = params$power_exp[1]) / 
      ( seq(0.01, 1.0, by = 0.01) + 0.2 + params$d[1] )
  )
  
  ### As it currently stands, these are not incorporating host tolerance into the mix
if (params$parasite_tuning[1] == FALSE) {
  
  gg4.1 <- ggplot(power_trade_dat, aes(alpha, beta)) + geom_line() + 
    geom_point(data = res_1000_all, aes(
      plogis(mean_plalpha)
    , power_tradeoff(plogis(mean_plalpha), params$power_c[1], params$power_exp[1]) * plogis(mean_plbeta)
      )) + 
    geom_point(data = res_1000_all, aes(
      plogis(mean_plalpha)
    , power_tradeoff(plogis(mean_plalpha), params$power_c[1], params$power_exp[1]) * plogis(mean_plbeta)
    , colour = fill_birth
    , shape = mut_mean
    ))
  
  gg4.2 <- ggplot(power_trade_dat, aes(alpha, R0)) + geom_line() + 
    geom_point(data = res_1000_all, aes(
      plogis(mean_plalpha)
    , ( power_tradeoff(plogis(mean_plalpha), params$power_c[1], params$power_exp[1]) * plogis(mean_plbeta) ) / 
      ( plogis(mean_plalpha) + 0.2 + params$d[1] )  
      )) + 
    geom_point(data = res_1000_all, aes(
      plogis(mean_plalpha)
    , ( power_tradeoff(plogis(mean_plalpha), params$power_c[1], params$power_exp[1]) * plogis(mean_plbeta) ) /
      ( plogis(mean_plalpha) + 0.2 + params$d[1] )
    , colour = fill_birth
    , shape = mut_mean
    ))
  
} else {
  
    gg4.1 <- ggplot(power_trade_dat, aes(alpha, beta)) + geom_line() + 
    geom_point(data = res_1000_all, aes(
      plogis(mean_plalpha)
    , power_tradeoff(plogis(mean_plalpha), params$power_c[1], params$power_exp[1]) *  exp(-(mean_plbeta - mean_plalpha)^2 / eff_scale)
      )) + 
    geom_point(data = res_1000_all, aes(
      plogis(mean_plalpha)
    , power_tradeoff(plogis(mean_plalpha), params$power_c[1], params$power_exp[1]) *  exp(-(mean_plbeta - mean_plalpha)^2 / eff_scale)
    , colour = fill_birth
    , shape = mut_mean
    ))
  
  gg4.2 <- ggplot(power_trade_dat, aes(alpha, R0)) + geom_line() + 
    geom_point(data = res_1000_all, aes(
      plogis(mean_plalpha)
    , ( power_tradeoff(plogis(mean_plalpha), params$power_c[1], params$power_exp[1]) *  exp(-(mean_plbeta - mean_plalpha)^2 / eff_scale) ) / 
      ( plogis(mean_plalpha) + 0.2 + params$d[1] )  
      )) + 
    geom_point(data = res_1000_all, aes(
      plogis(mean_plalpha)
    , ( power_tradeoff(plogis(mean_plalpha), params$power_c[1], params$power_exp[1]) *  exp(-(mean_plbeta - mean_plalpha)^2 / eff_scale) ) /
      ( plogis(mean_plalpha) + 0.2 + params$d[1] )
    , colour = fill_birth
    , shape = mut_mean
    ))
  
}
  
#  grid.arrange(gg1.2, gg4.1, gg3.1, gg1.1, gg4.2, gg3.2, gg1.3, gg1.4, gg2.1, ncol = 3, nrow = 3)
grid.arrange(gg1.2, gg4.1, gg1.2b, gg1.1, gg4.2, gg1.3b, gg1.3, gg1.4, gg1.4b, ncol = 3, nrow = 3)
re