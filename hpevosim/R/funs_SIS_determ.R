######
## ODE models
######

## The same as the RD functions in funs_SIR.R but with no death 

### Naming:
## wt vs nt = with vs no tradeoff
## wm vs nm = with vs no mutation

### Parameters:
## beta  : transmission of parasite strain X (linked to recovery below through tradeoff -- or not)
## gamma0: intrinsic host recovery
## gamma : recovery "boost" of an infected host with parasite strain X


## First models with tradeoff

evosim_determ.wt.nm <- function (t, y, parms) {
  
    with(c(as.list(parms))  , {
  
  S <- y[1] 
  I <- y[2:(N_strains + 1)]
  
    ## dS / dt (Susceptible hosts)
   dS <- - sum(beta * I) * S + sum(I * (gamma + gamma0))
    
    ## dI / dt (Infected hosts of one of many types)
   dI <- beta * I * S - I * (gamma + gamma0)
    
  list(c(dS, dI))
  
    }
      
      )
    
}
## SIS with mutation
evosim_determ.wt.wm <- function (t, y, parms) {
  
    with(c(as.list(parms))  , {
  
  S <- y[1] 
  I <- y[2:(N_strains + 1)]
  
    ## dS / dt (Susceptible hosts)
   dS <- - sum(beta * I) * S + sum(I * (gamma + gamma0))
    
    ## dI / dt (Infected hosts of one of many types)
   dI <- beta * I * S - I * (gamma + gamma0) + 
    tran.1D(
       I
     , D  = mut_rate
     , v  = biased_mut
     , dx = 1
     , flux.down = 0
     , flux.up   = 0)$dC
    
  list(c(dS, dI))
  
    }
      
      )    
      
}

## Second, models without tradeoff

evosim_determ.nt.nm <- function (t, y, parms) {
  
    with(c(as.list(parms))  , {
  
  S <- y[1] 
  I <- y[2:(N_strains + 1)]
  
    ## dS / dt (Susceptible hosts)
   dS <- - sum(beta * I) * S + sum(I * (gamma + gamma0))
    
    ## dI / dt (Infected hosts of one of many types)
   dI <- beta * I * S - I * (gamma + gamma0)
    
  list(c(dS, dI))
  
    }
      
      )
    
}
## SIS with mutation
evosim_determ.nt.wm <- function (t, y, parms) {
  
    with(c(as.list(parms))  , {
  
  S <- y[1] 
  I <- y[2:(N_strains + 1)]
  
    ## dS / dt (Susceptible hosts)
   dS <- - sum(beta * I) * S + sum(I * (gamma + gamma0))
    
    ## dI / dt (Infected hosts of one of many types)
   dI <- beta * I * S - I * (gamma + gamma0) + 
    tran.1D(
       I
     , D  = mut_rate
     , v  = biased_mut
     , dx = 1
     , flux.down = 0
     , flux.up   = 0)$dC
    
  list(c(dS, dI))
  
    }
      
      )    
      
}

## Deterministic version for calculating start vals
calc_startvals_determ              <- function (neg_trait0, pos_trait0, tuning, N, power_c, power_exp, mut_link_p, eff_scale, no_tradeoff, parasite_tuning) {

## No resistance in SIS, but leaving structure for now...
negtrait_r     <- mut_link_p$linkinv(neg_trait0) #- mut_link_h$linkfun(res0))

## Also no tolerance in SIS, but leaving structure for now...
negtrait_rt    <- negtrait_r # - mut_link_h$linkfun(tol0))

if (!no_tradeoff) {
if (parasite_tuning) {
## Symmetrical matrix of efficiencies associated with combination of each parasite trait
     effic       <- outer(neg_trait0, tuning, FUN = eff_calc, eff_scale = eff_scale)
## From these calculate beta of each of the strains
     joint_beta  <- sweep(effic, 2, power_tradeoff(alpha = negtrait_r, c = power_c, curv = power_exp), FUN = "*")
} else {
  if (!tradeoff_only) {
   # effic <- matrix(data = tuning, nrow = 1, ncol = length(neg_trait0))
    effic <- mut_link_p$linkinv(tuning)
  } else {
    ## Should refer to number of starting strains
   # effic <- matrix(data = 1, nrow = 1, ncol = length(neg_trait0))
    effic <- rep(1, length(neg_trait0))
  }
  
joint_beta <- outer(effic, power_tradeoff(alpha = negtrait_r, c = power_c, curv = power_exp))
  
} 
  
} else {
  
effic       <- rep(1, length(neg_trait0))
## Slightly funky here, but use neg_trait0 because it is set up with the same breaks as pos_trait
joint_beta  <- outer(effic, mut_link_p$linkinv(neg_trait0))

}

return(list(
  intrinsic_postrait = effic
, tuning             = tuning
, joint_postrait     = joint_beta
, joint_negtrait_rt  = negtrait_rt
  ))

}
## For plotting, want to calculate the surface on a linear scale, which is easier to look at
calc_startvals_determ_for_plotting <- function (neg_trait0, tuning, N, power_c, power_exp, mut_link_p, eff_scale) {

## No resistance in SIS, but leaving structure for now...
alpha_r     <- mut_link_p$linkinv(mut_link_p$linkfun(neg_trait0)) #- mut_link_h$linkfun(res0))

## Also no tolerance in SIS, but leaving structure for now...
alpha_rt    <- mut_link_p$linkinv(mut_link_p$linkfun(alpha_r)) #- mut_link_h$linkfun(tol0))

## calculate efficiency based on the matching of tuning and aggressiveness
eff_calc    <- function (x, y, eff_scale) {
  exp(-(x - y)^2 / eff_scale)
}
## Symmetrical matrix of efficiencies associated with combination of each parasite trait
effic       <- outer(mut_link_p$linkfun(neg_trait0), mut_link_p$linkfun(tuning), FUN = eff_calc, eff_scale = eff_scale)

## From these calculate beta of each of the strains
joint_beta  <- sweep(effic, 2, power_tradeoff(alpha = alpha_r, c = power_c, curv = power_exp), FUN = "*")

return(list(
  intrinsic_beta = effic
, tuning         = tuning
, joint_beta     = joint_beta
, joint_alpha    = alpha_rt
  ))

}
