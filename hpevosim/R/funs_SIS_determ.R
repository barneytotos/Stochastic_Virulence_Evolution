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

## Create a negtrait matrix from the negtrait vector
negtrait_mat <- matrix(
  data = rep(negtrait_rt, each = length(tuning))
, nrow = length(tuning), ncol = length(neg_trait0))
  
} 
  
} else {
  
effic       <- rep(1, length(neg_trait0))
## Slightly funky here, but use neg_trait0 because it is set up with the same breaks as pos_trait
joint_beta  <- outer(effic, mut_link_p$linkinv(neg_trait0))

## Create a negtrait matrix from the negtrait vector
negtrait_mat <- matrix(
  data = rep(negtrait_rt, each = length(tuning))
, nrow = length(tuning), ncol = length(neg_trait0)
, byrow = T)

}

return(list(
  intrinsic_postrait = effic
, tuning             = tuning
, joint_postrait     = joint_beta
, joint_negtrait     = negtrait_mat
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

run_sim <- function(
   no_tradeoff          = TRUE
 , nt_mut_var_pos_trait = TRUE    ## For nearly neutral model are we tracking evolution in transmission (pos_trait) (TRUE) or recovery (neg_trait) (FALSE)
 , pos_trait0           = 0.005 
    ## These first parameters are all about what type of simulation to run
     ## Note: hosts are always set to evolve. Crudely can force them never to evolve by just setting the probability to effectively 0
 , deterministic        = FALSE   ## Run an advection diffusion version?
 , parasite_tuning      = TRUE    ## Reformulation where parasite efficiency is defined as matching tuning and aggressiveness
 , tradeoff_only        = FALSE   ## Stepping back to ignore tuning. Parasite just evolving according to the tradeoff curve
 , agg_eff_adjust       = FALSE   ## For efficiency model but not tuning model. Does an increase in parasite aggressiveness decrease efficiency (as a cost)
    ## These next parameters all control how the simulation is run
 , R0_init              = 2       ## >1, not actually used for tuning model
 , neg_trait0           = 0.03    ## Either host recovery rate to the parasite (SIS model) or parasite-induced host mortality (SIR model). Default to main thesis result start
 , tune0                = 0.97    ## Starting tuning. Give a value, but only used if parasite_tuning == TRUE. Default to main thesis result start
 , gamma0               = 0.2     ## Background host recovery (immune pressure or however you want to think of it)
                                  ## (Can set to zero if using non-tradeoff model)
 , N                    = 200     ## Host population size, integer > 0
 , mu                   = 0.01    ## Mutation probability, > 0
 , mut_type             = "shift" ## Type of mutation supported. Only shift viable in this version of the code
 , mut_mean             = -1      ## Mutation mean, < 0 (for sensibility). Ignored for parasite_tuning == TRUE
 , mut_sd               = 0.5     ## Mutation sd, > 0
 , mut_link_p           = NULL    ## Default of logit scale setup in the function
# , mut_link_h           = NULL   ## Default of logit scale setup in the function
# , mut_host_sd_shift    = 1      ## **1 for identical sd to the parasite
# , mut_host_mean_shift  = 1      ## **1 for identical mean to the parasite
# , mut_host_mu_shift    = 2      ## Proportion less frequent host mutation is than parasite mutation, only used if hosts are evolving 
# , mut_host_res_bias    = 0.5    ## 0.5 means equal probability of evolving resistance or tolerance, larger proportion favors more resistance evolution
# , res0                = 1       ## Starting host mean resistance value
# , res0_sd             = 0       ## Variation in resistance among starting host strains (if > 1 host strain). Not yet used, but plan to
# , tol0                = 1       ## Starting host mean tolerance value
# , tol0_sd             = 0       ## Variation in tolerance among starting host strains (if > 1 host strain). Not yet used, but plan to
 , power_c              = 0.002   ## Power law tradeoff scaling
 , power_exp            = 3       ## Power law tradeoff exponent
 , eff_hit              = 1       ## Size of the negative correlation between parasite alpha and efficiency evolution (only used if agg_eff_adjust == TRUE)
 , Imat                 = NULL    ## Setup within function in this version, don't adjust
 , eff_scale            = 50      ## Weighting of the matching of parasite tuning and aggressiveness
 , nt                   = 100000  ## Length of simulation (time steps)
 , rptfreq              = max(nt / 500, 1) ## How often the state of the system is saved
 , seed                 = NULL    ## Can set seed if desired
 , progress             = "bar"   ## Progress bar?
 , debug                = FALSE   ## A debug option for stopping at various points in the sim. See code to check where.
 , debug2               = FALSE   ## A debug option for stopping at various points in the sim. See code to check where.
 , debug3               = FALSE   ## A debug option for stopping at various points in the sim. See code to check where.
 , debug4               = FALSE
 , debug4_val           = 1
    ## A few deterministic parameters
 , determ_length        = 200     ## length of time to run the deterministic model
 , determ_timestep      = 5       ## Lsoda parameter for RD model
 , lsoda_hini           = NULL    ## Lsoda parameter for RD model
 , Imat_seed            = NULL    ## Set up inside function, keep as null
  ) {

    if (round(N)!=N) {
        warning("rounding N")
        N <- round(N)
    }
    if (mut_mean>0) {
        warning("positive mutation bias")
    }
    stopifnot(
      R0_init > 0
    , N       > 0
    , mu      > 0
    , mut_sd  > 0
    , (nt/rptfreq) %% 1 == 0)

    dfun <- function(lab="") {
        if (debug) {
            cat(lab,"\n")
            print(state)
        }
        return(NULL)
    }

    if (!is.null(seed)) set.seed(seed)

   ## Set up mutation link
      ## Positive with log link can push gamma overboard when mutation in gamma is on average disadvantageous.
        if (is.null(mut_link_p)) mut_link_p <- make.link("cloglog")

    ## If deterministic == TRUE, set up the full array of possible strains from the start (for diffusion via mutation)
    if (!deterministic) {

        ## Based on my new setup it seems better to directly calculate efficiency and starting beta as a function of defined tuning
        ## and aggressiveness starting values instead of going backwards from R0

        startvals  <- calc_startvals(neg_trait0, tune0, N,
                                     power_c, power_exp,
                                     mut_link_p, eff_scale, parasite_tuning,
                                     tradeoff_only, no_tradeoff, pos_trait0)
    #    pos_trait0 <- startvals$joint_postrait
        
     if (is.null(Imat)) {
   ## start at equilibrium I ...
      Imat <- max(1, round(N*(1 - 1/R0_init)))
      Imat <- as.matrix(Imat)
    }

    Svec  <- N-Imat
    
    } else {

      neg_trait0     <- c(seq(
        mut_link_p$linkfun(0.01), mut_link_p$linkfun(0.99)
        , length = 100))
      
      ## Placeholder if the tuning model isn't being used
      tuning     <- c(seq(
        mut_link_p$linkfun(0.01), mut_link_p$linkfun(0.99)
        , length = 100))
    
      startvals  <- calc_startvals_determ(
         neg_trait0, pos_trait0, tuning, N, power_c
       , power_exp, mut_link_p, eff_scale, no_tradeoff, parasite_tuning)
      
      ## Also adjust Imat to capture how many total strains there are
      Imat       <- matrix(data = 0, ncol = length(neg_trait0), nrow = length(tuning))
      
        ## Later can set a parameter for which strains are the starting strains
        ## Imat[nrow(Imat), 1] <- N - sum(Svec)
        ## Imat[1, 1] <- N - sum(Svec)
      
      ## With the deterministic rate model need to start at a high enough beta to overcome gamma0
       ## need to make this a parameter, but for now this seems ok...
       Svec <- N - 10
      
       Imat[Imat_seed[1], Imat_seed[2]] <- N - Svec

      # Convert to proportions
       Imat <- Imat / N
       Svec <- Svec / N
       N <- 1
    }

    ## The trait that defines evolution of parasite efficiency will differ depending on the
    ##  model formulation
    ## Set up this way so that the rest of the code relies upon the same structure,
    ## and all the knobs can use the same code with fewer modifications
    if (no_tradeoff) {
      lpostrait  <- mut_link_p$linkfun(startvals$joint_postrait)
    } else {
    if (parasite_tuning || (!parasite_tuning && !tradeoff_only)) {
        if (!deterministic) {
            lpostrait  <- mut_link_p$linkfun(startvals$tuning)     
        } else {
            lpostrait  <- startvals$tuning 
        }
    } else {
            lpostrait  <- mut_link_p$linkfun(startvals$joint_postrait)
    }
    }

    ## negtait being either recovery (normally gamma) or in the case of SIR, alpha (intrinsic parasite mortality pressure)
        lnegtrait  <- mut_link_p$linkfun(startvals$joint_negtrait)

    ## Initial trait vectors for the host genotypes.
    ## Assumes all hosts start with identical trait value (for now)
    ##    hrtraitvec <- rep(mut_link_h$linkfun(res0), length(res0))
    ##    httraitvec <- rep(mut_link_h$linkfun(tol0), length(tol0))

    ## Without host resistance and tolerance evolution neg_trait0 is just what was defined
 #   neg_trait0     <-  startvals$joint_negtrait

    ## Parameters structure (parallel vectors), so these can
     ## Be modified via function and passed back ...
    state      <- list(
#      beta       = as.matrix(pos_trait0)
#    , alpha      = as.matrix(neg_trait0)
      pos_trait   = as.matrix(startvals$joint_postrait)
    , neg_trait   = as.matrix(startvals$joint_negtrait)      
    , lpostrait  = lpostrait
    , lnegtrait  = lnegtrait
#   , hrtraitvec = hrtraitvec
#   , httraitvec = httraitvec
    , Imat       = Imat
    , Svec       = Svec)
    
    dfun("init")

    nrpt <- nt %/% rptfreq

    ## FIXME: add, somehow, to documentation (is there a way to not repeat ourselves?)
     ## Added tracking host responses
    res <- as.data.frame(matrix(
      NA, nrow = nrpt, ncol = 18
    , dimnames = list(
      NULL
    , c("time"
      , "num_S"
      , "num_S_strains"
      , "num_I"
      , "num_I_strains"
      , "pop_size"
#     , paste0(c("mean_hl", "sd_hl"), "res")
#     , paste0(c("mean_hl", "sd_hl"), "tol")
#     , paste0(c("mean_pl", "sd_pl"), mut_var)
      , paste0(c("mean_pl", "sd_pl"), "negtrait")
      , "mean_negtrait"
      , "sd_negtrait"
      , "median_negtrait"
      , "lower_negtrait"
      , "upper_negtrait"
      , "mean_postrait"
      , "sd_postrait"
      , "total_mutations"
      , "total_extinctions"
      , "shannon"
      ))))
    
    mut_counter <- 0
    num_extinct <- 0

    ## Run the stochastic version
    if (!deterministic) {

    for (i in 1:nrpt) {
            for (j in 1:rptfreq) {

               ## Debug code chunk that can be cut and paste to wherever there is a problem
                ## Fill i and j in manually after checking output after an error is encountered
                if (debug2) {
                  print(paste(i, j, sep = "  -  "))
                  print(str(state$Svec))
                  assign("state_check", state, .GlobalEnv)
                  if (i == 20 & j == 20) browser()
                }

                ## [Step 1]: Infection.
                ## Prob of escaping infection completely
                uninf      <- with(state,
                       rbinom_mat(
                  n    = length(Svec)
                , size = Svec
                , prob = prod((1-pos_trait)^Imat)
                , nrow = nrow(Svec)
                , ncol = ncol(Svec))
                )
                
                ## Division of new infections among strains
                ## 'prob' is internally normalized
                newinf     <- with(state, get_inf(Svec = Svec, uninf = uninf, Imat = Imat, beta = pos_trait))
             
                ## sanity check
                stopifnot(sum(rowSums(newinf) + uninf) == sum(state$Svec))

                ## [Step 2]: Recovery of I.
                recover    <- rbinom_mat(
                  n    = length(c(state$Imat))
                , size = c(state$Imat)
                  ## total host recovery based on some background probability + whatever the parasite is doing
                  ## hazard-scale calculation: (1-gamma0) is probability of "no baseline recovery", 1-alpha is probability
                  ##  of "no strain-specific recovery"
                , prob =  1 - (1 - c(state$neg_trait)) * (1 - gamma0)
                , nrow = nrow(state$Imat)
                , ncol = ncol(state$Imat))

                ## [Step 3: host mortality] See "funs_SIR.R"
                
                ## [Step 4]: Update Svec with infections and recoveries (mutations are set asside)
                state$Svec <- state$Svec + rowSums(recover) - rowSums(newinf) 
                
                ## [Step 5.1]: Mutation of new infections. Fraction of new infections -> mutation
               if (sum(newinf) > 0) {
                   mutated  <- rbinom_mat(n = newinf, size = newinf, prob = mu, nrow = nrow(newinf), ncol = ncol(newinf))
               } else {
                   mutated <- matrix(0, nrow = nrow(newinf), ncol = 1)
               }

                ## debug to force mutation to check if that is working
                if (debug3) { 
                  mutated[1, 1] <- 2
                }
                
                mut_counter <- mut_counter + sum(mutated)

                ## sanity checks
                stopifnot(length(recover) == length(mutated))
                stopifnot(length(newinf)  == length(state$Imat))

                ## [Step 5.2]: Update Infecteds
                state$Imat <- state$Imat - recover + newinf - mutated

                dfun("before mutation")
                if (debug) print(mutated)

                ## [Step 5.3]: Find new phenotypes for mutated parasites and hosts.
                tot_mut <- sum(mutated)

                ## Original parasite traits that mutants arise from. Need this for later
                orig_trait <- list(
                   pos_trait = rep(state$lpostrait, colSums(mutated))
                 , neg_trait = rep(state$lnegtrait, colSums(mutated))
                  )

                ## [Step 6] Mut parasite, and update state
                if (tot_mut > 0) {
                    state <- do_mut(
                      state
                    , no_tradeoff          = no_tradeoff
                    , nt_mut_var_pos_trait = nt_mut_var_pos_trait
                    , orig_trait           = orig_trait  ## ^^Just care about intrinsic nature of a strain
                    , mut_mean             = mut_mean
                    , mut_sd               = mut_sd
                    , mut_type             = mut_type
                    , power_c              = power_c
                    , power_exp            = power_exp
                    , mut_link_p           = mut_link_p
                    , agg_eff_adjust       = agg_eff_adjust
                    , eff_hit              = eff_hit
                    , parasite_tuning      = parasite_tuning
                    , tradeoff_only        = tradeoff_only
                    , mutated              = mutated
                    , eff_scale            = eff_scale
                      )
                }

                if (sum(state$Imat)==0) {
                    message(sprintf("system went extinct prematurely (t=%d)", i))
                    break
                }

                ## Look over columns of I for extinct parasites and over rows of I and entries of S for extinct hosts
                extinct_p   <- which(colSums(state$Imat) == 0)
                num_extinct <- num_extinct + length(extinct_p)
#               extinct_h   <- which(rowSums(state$Imat) == 0 & state$Svec == 0)
                if (length(extinct_p) > 0) {
                    state   <- do_extinct(state, extinct = extinct_p, parasite = TRUE)
                }
                dfun("after mutation")

            }  ## rptfreq time steps
      
    if (debug4 && i == debug4_val) {
      browser()
    }
      
        ## summary statistics
        I_tot        <- ncol(state$Imat)
        num_I        <- sum(state$Imat)
        ltrait_mean  <- sum(colSums(state$Imat)*state$lpostrait)/num_I
        ## Not quite sure about this formula. Taken from BMB code
        ltrait_sd    <- sqrt(sum(colSums(state$Imat)*(state$lpostrait-ltrait_mean)^2)/num_I)
        lalpha_mean  <- sum(colSums(state$Imat)*state$lnegtrait)/num_I
        lalpha_sd    <- sqrt(sum(colSums(state$Imat)*(state$lnegtrait-lalpha_mean)^2)/num_I)
        lalpha_q     <- quantile(rep(state$lnegtrait, colSums(state$Imat)), c(0.025, 0.50, 0.975))
        lalpha_est   <- lalpha_q[2]
        lalpha_lwr   <- lalpha_q[1]
        lalpha_upr   <- lalpha_q[3]
        num_S        <- sum(state$Svec)
        S_tot        <- length(state$Svec)
#       lhres_mean   <- sum(state$Svec*state$hrtraitvec)/num_S
#       lhres_sd     <- sqrt(sum(state$Svec*(state$hrtraitvec-lhres_mean)^2)/num_S)
#       lhtol_mean   <- sum(state$Svec*state$httraitvec)/num_S
#       lhtol_sd     <- sqrt(sum(state$Svec*(state$httraitvec-lhtol_mean)^2)/num_S)
        pop_size     <- num_I + num_S
        shann        <- vegan::diversity(state$Imat, index = "shannon", MARGIN = 1, base = exp(1))

        ## actual alpha and beta of all parasites in all host classes
        avg_alpha    <- mean(state$neg_trait)
        sd_alpha     <- sd(state$neg_trait)
        avg_beta     <- mean(state$pos_trait)
        sd_beta      <- sd(state$pos_trait)

        if (progress == "bar") {
          cat(".")
        } else if (progress == "text") {
          if (((i / 50) %% 1) == 0) {
          print(paste(round(i/nrpt, 2)*100, "% Complete"))            
          }
        } else {
          
        }
        res[i,] <- c(
          i*rptfreq
        , num_S
        , S_tot
        , num_I
        , I_tot
        , pop_size
#       , lhres_mean
#       , lhres_sd
#       , lhtol_mean
#       , lhtol_sd
#       , ltrait_mean
#       , ltrait_sd
        , lalpha_mean
        , lalpha_sd
        , avg_alpha
        , sd_alpha
        , lalpha_est
        , lalpha_lwr
        , lalpha_upr
        , avg_beta
        , sd_beta
#       , ifelse(mut_var == "beta", state$gamma[1,1], state$gamma[1,1])
        , mut_counter
        , num_extinct
        , shann
          )

        ## DRY ...
        if (sum(state$Imat) == 0) {
            message(sprintf("system went extinct prematurely (t=%d)", i))
            break
        }

    } ## loop over reporting frequencies
        if (progress=="bar") {
            cat("\n")
        }
    class(res) <- c("hpevosim","data.frame")   
    return(res)

    ## Run the deterministic model
    } else {
      
      ## convert state into a state list with only S and I and everything else put into params
      state_determ <- c(
        Svec = state$Svec
      , Imat = state$Imat
      )
      
      ## Setup in what traits mutation occurs
      if (!no_tradeoff) {
        if (nt_mut_var_pos_trait) {
        mutpos <- 1 
        mutneg <- 0
        } else {
        mutpos <- 0
        mutneg <- 1        
        }
      } else {
        mutpos <- 1 
        mutneg <- 1      
      }
      
      params_determ <- list(
        N         = N
      , beta      = state$pos_trait
      , gamma     = state$neg_trait
      , gamma0    = gamma0
      , mutlev    = mu
      , with_mut  = TRUE
      , mutpos    = mutpos
      , mutneg    = mutneg
      )
      
    hpevosim_determ.out <- as.data.frame(
      ode(
    y        = state_determ
  , times    = seq(0, determ_length, by = determ_timestep)
  , func     = hpevosim_determ
  , parms    = params_determ
  , maxsteps = 5000
  , hini     = lsoda_hini
  , method   = 'rk4'
    ))
  
        return(hpevosim_determ)

    }
}

