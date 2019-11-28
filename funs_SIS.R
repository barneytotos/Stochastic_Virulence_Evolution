#######
## SIS model with no host demographics. Tradeoff involves just transmission and recovery. No mortality.
## Top level simulation function as well as all helper functions are contained in this script. 
## For extensions on the simplest model and other code originally written by BMB email Ben Bolker for access
#######

#######
## Cautions 
#######

## (0) Residual code lingers for resistance and tolerance evolution. Uncommenting won't allow it to run as some has been
 ## removed. See funs_SIR.R for working resistance and tolerance evolution
## (1) Just alpha and gamma evolving
## (2) If you choose tradeoff_only == TRUE, ltrait doesn't mean anything, beta is what is important.
 ## However the summary stuff at the each trait matrix storing step isn't currently set up to handle this
  ## correctly. Need to calculate instrinsic beta from the tradeoff curve (especially if host evolution
   ## is happening and/or if non 0 resistance and tolerance are used)
## (3) Don't give more than one host or pathogen genotype as starting values, this is hypothetically coded
 ## but will currently run into problems establishing starting conditions. 
## (4) Host resistance and tolerance avialable as parameters, and can also be set to evolve, but are best set to
 ## some constant value (all results in the thesis chapter in this repo and the poster are with 0 resistance
 ## and 0 tolerance) as they currently evolve without constraint, leading to uninteresting results and possibly error
## (5) Stochastic vs Determinisitc option should porbably be separated. Works, but the top level sim
 ## function is a bit too long for my liking at present 

######
## Accessory Functions
######
get_mut_p        <- function (orig_trait, mut_var, power_c, power_exp, mut_link_p, mut_mean
  , mut_sd, mut_type = "shift", agg_eff_adjust, eff_hit, parasite_tuning, tradeoff_only) {

   if (mut_type == "shift") {

      ## Assume a neutrally evolving aggressiveness, regardless of method of calculating efficiency
       ## recall I call virulence the parasite's "negative" trait (negative trait could also be host
       ## recover, but that isn't prepped yet)
     if (parasite_tuning == TRUE) {
       neg_trait_adj <- rnorm(length(orig_trait$neg_trait), 0, mut_sd)
       new_trait_neg <- orig_trait$neg_trait + neg_trait_adj
     } else {
       ## Update Nov 1, 2019: Seems that a given parasite trait evolving should push towards faster host killing.
        ## Also, this mut_mean is a pretty massive mutation size... Parameter vals are a bit hard to have for both models
         ## because mut_sd will mean something pretty different for the different models (e.g. for just tradeoff curve, 
          ## evolving from underneath optimal alpha will very quickly get a parasite to opt alpha if mu has any real size)
           ## As a placeholder for now just to get the plots to make sure code is working is just to drop the size a bit
       neg_trait_adj <- rnorm(length(orig_trait$neg_trait), mut_mean*-1, mut_sd)
       new_trait_neg <- orig_trait$neg_trait + neg_trait_adj       
     }

## Three options for the evolution of the parasite's "positive" trait
## *NOTE*: Positive trait means different things depending on the model.
 ## In the tradeoff curve model and the tuning == FALSE model it is beta, in tuning == TRUE model it is tuning
  ## (1) Evolve efficiency directly. Here assume biased efficiency mutations (requires parasite_tuning == FALSE)
   ## (1.1) agg_eff_adjust == TRUE: With a direct adjustment to efficiency due to the evolution of aggressiveness above
   ## (1.2) agg_eff_adjust == FALSE: Independent evolution of efficiency
  ## (2) Evolve in tuning, which is used to calculate efficiency. In this formulation change in aggressiveness always 
    ## lowers parasite efficiency indirectly, regardless of direction, because of a poorer match to tuning.
    ## For now assume independetly evolving tuning

    if (parasite_tuning == FALSE) {

    ## tradeoff_only controls whether parasites evolve only with a tradeoff curve and none of the efficiency stuff
      if (tradeoff_only == FALSE) {
           
        if (agg_eff_adjust == TRUE) {

      new_trait_pos <- orig_trait$pos_trait - abs(neg_trait_adj) * eff_hit

          ## Assume a negatively evolving efficiency
        new_trait_pos <- new_trait_pos +
            rnorm(length(new_trait_pos),
              ifelse(mut_var == "beta",  ## Beta and gamma in opposite directions
               mut_mean
            ,  mut_mean*-1)
            , mut_sd
              )

       } else {

         new_trait_pos <- orig_trait$pos_trait +
            rnorm(length(orig_trait$pos_trait),
              ifelse(mut_var == "beta",  ## Beta and gamma in opposite directions
               mut_mean
            ,  mut_mean*-1)
            , mut_sd)

       }

    ## Only a tradeoff curve. No evolution directly in beta. Beta is given by the tradeoff curve.
         } else {
           
     new_trait_pos <- orig_trait$pos_trait  
     
         }
    ## Parasite Tuning == TRUE
       } else {

         ## Assume neutrally evolving tuning
        new_trait_pos <- orig_trait$pos_trait +
            rnorm(length(orig_trait$pos_trait),
              0
            , mut_sd
              )
       }

        if (any(is.na(c(new_trait_pos, new_trait_neg)))) stop("??")
        return(list(new_trait_pos = new_trait_pos, new_trait_neg = new_trait_neg))
    } else stop("unknown mut_type")
}
do_mut           <- function (state, mut_var, orig_trait, ...) {
    new_trait        <- get_mut_p(orig_trait, mut_var, ...)
    state$ltraitvec  <- c(state$ltraitvec, new_trait$new_trait_pos)
    state$palphavec  <- c(state$palphavec, new_trait$new_trait_neg)
    return(state)
}
## Update mutant strain's trait values using power-law tradeoff. 
update_mut_pt    <- function (state, orig_trait, power_c, power_exp, mut_link_p, mutated, mutated_host, mut_var, ...) {

   ## Scale beta according to tradeoff curve
   new_par_beta <- scale_beta_alpha(state, power_c, power_exp, mut_link_p, ...)

   ## Resistance will act to decrease parasite transmission and virulence following the shape of the tradeoff curve
     ## Need to think criticall about what scale this should be conducted on. Both logit and probability scale feel like
      ## they each have problems
   ## First calculate a tradeoff curve with the same curvature that passes through the parasite's (alpha, beta)
   cvec         <- pt_calc_c(beta = new_par_beta, alpha = mut_link_p$linkinv(state$palphavec), curv = power_exp)

   ## For each of these tradeoff curves, calculate a new alpha and beta for each host that is infected
   new_alphas   <- t(mut_link_p$linkinv(state$palphavec))
   new_betas    <- matrix(power_tradeoff(
     alpha = c(new_alphas)
  ,  c     = rep(cvec, each = nrow(new_alphas)) 
  ,  curv  = power_exp)
    , nrow = nrow(new_alphas), ncol = ncol(new_alphas))

   state$alpha  <- new_alphas
   state$beta   <- new_betas

   ## Update Infected matrix with the new strain, maintaining which S class received that mutation.
    ## For each mutated strain, Imat gets a new column with a single 1, in the row in which the mutation occurred
   if (sum(mutated) > 0) {
     ## First make a matrix of 0s, then add a single one in each column corresponding to the row of the host that the parasite mutated in
     new_mutes            <- matrix(data = 0, nrow = nrow(state$Imat), ncol = sum(mutated))
     num_mutes            <- matrix(c(rep(which(rowSums(mutated) > 0), rowSums(mutated)[rowSums(mutated) != 0])
       , 1:sum(mutated)), ncol = 2, nrow = sum(mutated))
     new_mutes[num_mutes] <- 1
     state$Imat           <- cbind(state$Imat, new_mutes)
   }

   return(state)
}
## Remove extinct strains
do_extinct       <- function (state, mut_var, extinct, parasite) {

    state[[mut_var]] <- state[[mut_var]][,-extinct, drop = FALSE]
    state$alpha      <- state$alpha[,-extinct, drop = FALSE]
    state$ltraitvec  <- state$ltraitvec[-extinct, drop = FALSE]
    state$palphavec  <- state$palphavec[-extinct, drop = FALSE]
    state$Imat       <- state$Imat[,-extinct, drop = FALSE]

    return(state)
}
## This function is borderline hideous but it does what I want it to do (mutlinomial draws and return proper
 ## matrix structure). All ears for a better way to write this
get_inf          <- function (Svec, uninf, Imat, beta) {

matrix(
data =
    apply(
  Svec - uninf
, 2
, function (x) colSums(
matrix(
data = rmultinom(
   1
 , size = x
 , prob = beta * Imat
   )
, ncol = ncol(beta)
, nrow = nrow(beta)
      )
    )
  )
, ncol = ncol(beta)
, nrow = nrow(beta)
, byrow = TRUE
)

}
## calculate c for a power law that goes through the current strain | efficiency
pt_calc_c        <- function (alpha, beta, curv) {
 beta / ( alpha ^ (1 / curv) )
}
## Scale parasite beta and alpha evolution by the tradeoff
scale_beta_alpha <- function (state, power_c, power_exp, mut_link_p, parasite_tuning, tradeoff_only, eff_scale, ...) {

## Determine the beta of a given pathogen strain | that pathogen's current alpha value
 ## Maximum possible beta
max_beta      <- power_tradeoff(c = power_c, alpha = mut_link_p$linkinv(state$palphavec), curv = power_exp)

## realized beta. If efficiency is directly the trait evolving use the "positive" trait, otherwise calculate efficiency from tuning
if (parasite_tuning == FALSE & tradeoff_only == TRUE) {
realized_beta <- max_beta
} else {
  if (parasite_tuning == TRUE) {
realized_beta <- exp(-(state$ltraitvec - state$palphavec)^2 / eff_scale)  * max_beta
  } else {
realized_beta <- mut_link_p$linkinv(state$ltraitvec) * max_beta  
}
}

return(realized_beta)

}
## Calculate starting trait values for parasite and host | desired starting R0
calc_startvals        <- function (alpha0, tuning, N, power_c, power_exp, mut_link_p, eff_scale, parasite_tuning, tradeoff_only) {

## No resistance in SIS, but leaving structure for now...
alpha_r     <- mut_link_p$linkinv(mut_link_p$linkfun(alpha0)) #- mut_link_h$linkfun(res0))

## Also no tolerance in SIS, but leaving structure for now...
alpha_rt    <- mut_link_p$linkinv(mut_link_p$linkfun(alpha_r))# - mut_link_h$linkfun(tol0))

## Symmetrical matrix of efficiencies associated with combination of each parasite trait
if (parasite_tuning == TRUE) {
effic       <- outer(mut_link_p$linkfun(alpha0), mut_link_p$linkfun(tuning), FUN = eff_calc, eff_scale = eff_scale)
} else {
  if (tradeoff_only == FALSE) {
    effic <- matrix(data = tuning, nrow = 1, ncol = 1)
  } else {
    ## Should refer to number of starting strains
    effic <- matrix(data = 1, nrow = 1, ncol = 1)
  }
}


## From these calculate beta of each of the strains
joint_beta  <- sweep(effic, 2, power_tradeoff(alpha = alpha_r, c = power_c, curv = power_exp), FUN = "*")

return(list(
  intrinsic_beta = effic
, tuning         = tuning
, joint_beta     = joint_beta
, joint_alpha    = alpha_rt
  ))

}
## Determinisitc version for calculating start vals
calc_startvals_determ <- function (alpha0, tuning, N, power_c, power_exp, mut_link_p, eff_scale) {

## No resistance in SIS, but leaving structure for now...
alpha_r     <- mut_link_p$linkinv(alpha0)# - mut_link_h$linkfun(res0))

## Also no tolerance in SIS, but leaving structure for now...
alpha_rt    <- mut_link_p$linkinv(mut_link_p$linkfun(alpha_r))# - mut_link_h$linkfun(tol0))

## Symmetrical matrix of efficiencies associated with combination of each parasite trait
effic       <- outer(alpha0, tuning, FUN = eff_calc, eff_scale = eff_scale)

## From these calculate beta of each of the strains
joint_beta  <- sweep(effic, 2, power_tradeoff(alpha = alpha_r, c = power_c, curv = power_exp), FUN = "*")

return(list(
  intrinsic_beta = effic
, tuning         = tuning
, joint_beta     = joint_beta
, joint_alpha    = alpha_rt
  ))

}
## For plotting, want to calculate the surface on a linear scale, which is easier to look at
calc_startvals_determ_for_plotting <- function (alpha0, tuning, N, power_c, power_exp, mut_link_p, eff_scale) {

## No resistance in SIS, but leaving structure for now...
alpha_r     <- mut_link_p$linkinv(mut_link_p$linkfun(alpha0)) #- mut_link_h$linkfun(res0))

## Also no tolerance in SIS, but leaving structure for now...
alpha_rt    <- mut_link_p$linkinv(mut_link_p$linkfun(alpha_r)) #- mut_link_h$linkfun(tol0))

## calculate efficiency based on the matching of tuning and aggressiveness
eff_calc    <- function (x, y, eff_scale) {
  exp(-(x - y)^2 / eff_scale)
}
## Symmetrical matrix of efficiencies associated with combination of each parasite trait
effic       <- outer(mut_link_p$linkfun(alpha0), mut_link_p$linkfun(tuning), FUN = eff_calc, eff_scale = eff_scale)

## From these calculate beta of each of the strains
joint_beta  <- sweep(effic, 2, power_tradeoff(alpha = alpha_r, c = power_c, curv = power_exp), FUN = "*")

return(list(
  intrinsic_beta = effic
, tuning         = tuning
, joint_beta     = joint_beta
, joint_alpha    = alpha_rt
  ))

}
## Wrappers for selecting random numbers from matrices and returning a matrix
rbinom_mat  <- function (n, size, prob, nrow, ncol) {
  matrix(rbinom(n, size, prob), nrow = nrow, ncol = ncol)
}
## calculate efficiency based on the matching of tuning and aggressiveness
eff_calc    <- function (x, y, eff_scale) {
  exp(-(x - y)^2 / eff_scale)
}

######
## RD functions
######
## The same as the RD functions in funs_SIR.R but with no death here
## SIS no mutation
epi.fun.nm   <- function (t, y, parms) {
  
    with(c(as.list(parms))  , {
  
  S <- y[1] 
  I <- y[2:(N_strains + 1)]
  
    ## dS / dt (Susceptible hosts)
   dS <- - sum(beta * I) * S + sum(I * (alpha + gamma))
    
    ## dI / dt (Infected hosts of one of many types)
   dI <- beta * I * S - I * (alpha + gamma)
    
  list(c(dS, dI))
  
    }
      
      )
    
}
## SIS with mutation
epi.fun.wm   <- function (t, y, parms) {
  
    with(c(as.list(parms))  , {
  
  S <- y[1] 
  I <- y[2:(N_strains + 1)]
  
    ## dS / dt (Susceptible hosts)
   dS <- - sum(beta * I) * S + sum(I * (alpha + gamma))
    
    ## dI / dt (Infected hosts of one of many types)
   dI <- beta * I * S - I * (alpha + gamma) + 
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

######
## power tradeoff and R0 functions
######
power_tradeoff <- function (alpha, c, curv) {
  c * alpha ^ (1 / curv)
}
power_R0       <- function (alpha, c, curv, gamma, N) {
  ( N * c * alpha ^ (1 / curv) )  / ( alpha + gamma )
}
power_R0_slope <- function (alpha, c, curv, gamma, N) {
  ( N * c * alpha ^ (1 / (curv - 1)) * (gamma - curv * alpha + alpha) )  / ( curv * (alpha + gamma)^2 )
}
power_R0_grad  <- function (alpha, c, curv, gamma, N, eps) {
 (power_R0(alpha, c, curv, gamma, N) + power_R0_slope(alpha + eps, c, curv, gamma, N)) / 
   power_R0(alpha, c, curv, gamma, N) 
}

######
## AD functions to come later possibly
######

######
## Top level sim function (wrapper)
######
run_sim <- function(
    ## These first parameters are all about what type of simulation to run
     ## Note: hosts are always set to evolve. Crudly can force them never to evolve by just setting the probability to effectively 0
   deterministic        = FALSE   ## Run an advection diffusion version?
 , parasite_tuning      = TRUE    ## Reformulation where parasite efficiency is defined as matching tuning and aggressiveness
 , tradeoff_only        = FALSE   ## Stepping back to ignore tuning. Parasite just evolving according to the tradeoff curve
 , agg_eff_adjust       = FALSE   ## For efficiency model but not tuning model. Does an increase in parasite aggressiveness decrease efficiency (as a cost)
    ## These next parametes all control how the simulation is run
 , R0_init              = 2       ## >1, not actually used for tuning model
 , alpha0               = 0.03    ## Intrinsic parasite mortality probability (without influence of hosts). Default to main thesis result start
 , tune0                = 0.97    ## Strating tuning. Give a value, but only used if parasite_tuning == TRUE. Default to main thesis result start
 , gamma0               = 0.2     ## Background host recovery (immune pressure or however you want to think of it)
 , N                    = 200     ## Host population size, integer > 0
 , mu                   = 0.01    ## Mutation probability, > 0
 , mut_type             = "shift" ## Type of mutation supported. Only shift viable in this version of the code
 , mut_mean             = -1      ## Mutation mean, < 0 (for sensibility). Ignored for parasite_tuning == TRUE
 , mut_sd               = 0.5     ## Mutation sd, > 0
 , mut_var              = "beta"  ## Trait for parasite evolution. Only beta allowed for this code iteration
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
    if (mut_var=="beta") {
        if (is.null(mut_link_p)) mut_link_p <- make.link("cloglog")
  #     if (is.null(mut_link_h)) mut_link_h <- make.link("log")
      } else {
      ## Positive with log link can push gamma overboard when mutation in gamma is on average disadvantageous.
       ## Made to logit link for now, not completely convinced...
        if (is.null(mut_link_p)) mut_link_p <- make.link("cloglog")
  #     if (is.null(mut_link_h)) mut_link_h <- make.link("log")
      }

    ## If deterministic == TRUE, set up the full array of possible strains from the start (for diffusion via mutation)
    if (deterministic == FALSE) {

   ## Based on my new setup it seems better to directly calculate efficiency and starting beta as a function of defined tuning
    ## and aggressiveness starting values instead of going backwards from R0
    startvals  <- calc_startvals(alpha0, tune0, N, power_c, power_exp, mut_link_p, eff_scale, parasite_tuning, tradeoff_only)
    beta0      <- startvals$joint_beta

     if (is.null(Imat)) {
   ## start at equilibrium I ...
      Imat <- max(1, round(N*(1 - 1/R0_init)))
      Imat <- as.matrix(Imat)
    }

    Svec  <- N-Imat
    
    
    } else {

      alpha0     <- c(seq(
        mut_link_p$linkfun(0.01), mut_link_p$linkfun(0.99)
        , length = 100))
      tuning     <- c(seq(
        mut_link_p$linkfun(0.01), mut_link_p$linkfun(0.99)
        , length = 100))
      startvals  <- calc_startvals_determ(
         alpha0, tuning, res0, tol0, N, power_c
       , power_exp, mut_link_p, eff_scale
        )
      beta0      <- startvals$joint_beta

      ## Also adjust Imat to capture how many total strains there are
      Imat       <- matrix(data = 0, ncol = length(alpha0), nrow = length(tuning))
      ## Later can set a parameter for which strains are the starting strains
       # Imat[nrow(Imat), 1] <- N - sum(Svec)
       # Imat[1, 1] <- N - sum(Svec)
      ## With the deterministic rate model need to start at a high enough beta to overcome gamma0
       Imat[Imat_seed[1], Imat_seed[2]] <- N - sum(Svec)

      # Convert to proportions
       Imat <- Imat / N
       Svec <- Svec / N
       N <- 1
    }

    ## The trait that defines parasite efficiency that is evolving will differ depending on the formulation of
     ## how a parasite is evolving. 
    ## Set up this way so that the rest of the code relies upon the same structure, and all the knobs can use the
     ## same code with fewer modifications (does increase confusion at places I suppose)
    if (parasite_tuning == TRUE | (parasite_tuning == FALSE & tradeoff_only == FALSE)) {
      if (deterministic == FALSE) {
      ltraitvec  <- mut_link_p$linkfun(startvals$tuning)     
      } else {
      ltraitvec  <- startvals$tuning 
      }
    } else {
     ltraitvec  <- mut_link_p$linkfun(startvals$joint_beta)
    }

    ## Alpha (intrinsic parasite mortality pressure)
    if (deterministic == FALSE) {
    palphavec  <- mut_link_p$linkfun(alpha0)
    } else {
    palphavec  <- alpha0 
    }

    ## Initial trait vectors for the host genotypes. Assumes all hosts start with identical trait value (for now)
#    hrtraitvec <- rep(mut_link_h$linkfun(res0), length(res0))
#    httraitvec <- rep(mut_link_h$linkfun(tol0), length(tol0))

    ## Without host resistance and tolerance evolution alpha0 is just what was defined
    alpha0     <-  startvals$joint_alpha

    ## Parameters structure (parallel vectors), so these can
     ## Be modified via function and passed back ...
    state      <- list(
      beta       = as.matrix(beta0)
    , alpha      = as.matrix(alpha0)
    , ltraitvec  = ltraitvec
    , palphavec  = palphavec
#    , hrtraitvec = hrtraitvec
#    , httraitvec = httraitvec
    , Imat       = Imat
    , Svec       = Svec)


    dfun("init")

    nrpt <- nt %/% rptfreq

     ## Added tracking host responses
    res <- as.data.frame(matrix(
      NA, nrow = nrpt, ncol = 17
    , dimnames = list(
      NULL
    , c("time"
      , "num_S"
      , "num_S_strains"
      , "num_I"
      , "num_I_strains"
      , "pop_size"
#      , paste0(c("mean_hl", "sd_hl"), "res")
#      , paste0(c("mean_hl", "sd_hl"), "tol")
      , paste0(c("mean_pl", "sd_pl"), mut_var)
      , paste0(c("mean_pl", "sd_pl"), "alpha")
      , "mean_alpha"
      , "sd_alpha"
      , "median_alpha"
      , "lower_alpha"
      , "upper_alpha"
      , "mean_beta"
      , "sd_beta"
      ))))

    ## Run the stochastic version
    if (deterministic == FALSE) {

    for (i in 1:nrpt) {

            for (j in 1:rptfreq) {

               ## Debug code chunck that can be cut and paste to wherever there is a problem
                ## Fill i and j in manually after checking output after an error is encountered
                if (debug2 == TRUE) {
                  print(paste(i, j, sep = "  -  "))
                  print(str(state$Svec))
                  assign("state_check", state, .GlobalEnv)
                  if (i == 20 & j == 20) browser()
                }

                ## [Step 1]: Infection.
                ## Prob of escaping infection completely
                uninf      <- rbinom_mat(
                  n    = length(state$Svec)
                , size = state$Svec
                , prob = prod((1-state$beta)^state$Imat)
                , nrow = nrow(state$Svec)
                , ncol = ncol(state$Svec))
                
                ## Division of new infectives among strains
                 ## 'prob' is internally normalized
                newinf     <- get_inf(state$Svec, uninf, state$Imat, beta = state$beta)
                stopifnot(sum(rowSums(newinf) + uninf) == sum(state$Svec))

                ## [Step 2]: Recovery of I.
                recover    <- rbinom_mat(
                  n    = length(c(state$Imat))
                , size = c(state$Imat)
                  ## total host recovery based on some background probability + whatever the parasite is doing
                , prob =  1 - (1 - c(state$alpha)) * (1 - gamma0)
                , nrow = nrow(state$Imat)
                , ncol = ncol(state$Imat))

                ## [Step 4.1]: Mutation of new infections. Fraction of new infections -> mutation
               if (sum(newinf) != 0) {
                  mutated  <- rbinom_mat(n = newinf, size = newinf, prob = mu, nrow = nrow(newinf), ncol = ncol(newinf))
               } else {
                  mutated  <- rbinom_mat(n = nrow(newinf), size = newinf, prob = mu, nrow = nrow(newinf), ncol = ncol(newinf))
               }

                ## debug to force mutation to check if that is working
                if (debug3 == TRUE) { 
                  mutated[1, 1] <- 2
                  }

                stopifnot(length(recover) == length(mutated))
                stopifnot(length(newinf)  == length(state$Imat))

                ## [Step 4.3]: Update Infecteds
                state$Imat <- state$Imat - recover + newinf - mutated

                dfun("before mutation")
                if (debug) print(mutated)

                ## [Step 4.4]: Find new phenotypes for mutated parasites and hosts.
                tot_mut <- sum(mutated)

                ## Original parasite taits that mutants arise from. Need this for later
                orig_trait <- list(
                   pos_trait = rep(state$ltraitvec, colSums(mutated))
                 , neg_trait = rep(state$palphavec, colSums(mutated))
                  )

                ## Mut parasite first. Another choice of order that may/may not matter
                if (tot_mut > 0) {
                    state <- do_mut(
                      state
                    , mut_var         = mut_var
                    , orig_trait      = orig_trait  ## ^^Just care about intrinsic nature of a strain
                    , mut_mean        = mut_mean
                    , mut_sd          = mut_sd
                    , mut_type        = mut_type
                    , power_c         = power_c
                    , power_exp       = power_exp
                    , mut_link_p      = mut_link_p
                    , agg_eff_adjust  = agg_eff_adjust
                    , eff_hit         = eff_hit
                    , parasite_tuning = parasite_tuning
                    , tradeoff_only   = tradeoff_only
                      )
                }

                ## [Step 5]: Update Svec with infections and recoveries prior to the mutations
                state$Svec <- state$Svec + rowSums(recover) - rowSums(newinf) 

                ## [Step 6]: Update state (a big component beign parasite alpha and beta) with new parasite and host evolution.
                if (tot_mut > 0) {
                state <- update_mut_pt(
                    state           = state
                  , orig_trait      = orig_trait
                  , power_c         = power_c
                  , power_exp       = power_exp
                  , mut_link_p      = mut_link_p
                  , mutated         = mutated
                  , mut_var         = mut_var
                  , parasite_tuning = parasite_tuning
                  , tradeoff_only   = tradeoff_only
                  , eff_scale       = eff_scale)
                }

                if (sum(state$Imat)==0) {
                    message(sprintf("system went extinct prematurely (t=%d)", i))
                    break
                }

                ## Look over columns of I for extinct parasites and over rows of I and entries of S for extinct hosts
                extinct_p <- which(colSums(state$Imat) == 0)
#                extinct_h <- which(rowSums(state$Imat) == 0 & state$Svec == 0)
                if (length(extinct_p) > 0) {
                    state <- do_extinct(state, mut_var, extinct = extinct_p, parasite = TRUE)
                }
                dfun("after mutation")

            }  ## rptfreq time steps
      
    if (debug4 == TRUE & i == debug4_val) {
      browser()
    }
      
        ## summary statistics
        I_tot        <- ncol(state$Imat)
        num_I        <- sum(state$Imat)
        ltrait_mean  <- sum(colSums(state$Imat)*state$ltraitvec)/num_I
        ## Not quite sure about this formula. Taken from BMB code
        ltrait_sd    <- sqrt(sum(colSums(state$Imat)*(state$ltraitvec-ltrait_mean)^2)/num_I)
        lalpha_mean  <- sum(colSums(state$Imat)*state$palphavec)/num_I
        lalpha_sd    <- sqrt(sum(colSums(state$Imat)*(state$palphavec-lalpha_mean)^2)/num_I)
        lalpha_q     <- quantile(rep(state$palphavec, colSums(state$Imat)), c(0.025, 0.50, 0.975))
        lalpha_est   <- lalpha_q[2]
        lalpha_lwr   <- lalpha_q[1]
        lalpha_upr   <- lalpha_q[3]
        num_S        <- sum(state$Svec)
        S_tot        <- length(state$Svec)
#        lhres_mean   <- sum(state$Svec*state$hrtraitvec)/num_S
#        lhres_sd     <- sqrt(sum(state$Svec*(state$hrtraitvec-lhres_mean)^2)/num_S)
#        lhtol_mean   <- sum(state$Svec*state$httraitvec)/num_S
#        lhtol_sd     <- sqrt(sum(state$Svec*(state$httraitvec-lhtol_mean)^2)/num_S)
        pop_size     <- num_I + num_S

        ## actual alpha and beta of all parasites in all host classes
        avg_alpha    <- mean(state$alpha)
        sd_alpha     <- sd(state$alpha)
        avg_beta     <- mean(state$beta)
        sd_beta      <- sd(state$beta)

        if (progress == "bar") {
          cat(".")
        } else if (progress == "text") {
          if (((i / 50) %% 1) == 0) {
          print(paste(round(i/nrpt, 2), "% Complete"))            
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
#        , lhres_mean
#        , lhres_sd
#        , lhtol_mean
#        , lhtol_sd
        , ltrait_mean
        , ltrait_sd
        , lalpha_mean
        , lalpha_sd
        , avg_alpha
        , sd_alpha
        , lalpha_est
        , lalpha_lwr
        , lalpha_upr
        , avg_beta
        , sd_beta
     #   , ifelse(mut_var == "beta", state$gamma[1,1], state$gamma[1,1])
          )

        ## DRY ...
        if (sum(state$Imat) == 0) {
            message(sprintf("system went extinct prematurely (t=%d)", i))
            break
        }

    } ## loop over reporting frequencies
   # if (progress) cat("\n")
    return(res)

    ## Run the deterministic model
    } else {
      
      ## parameter values and setup. *!* Will /need to move lots of this up into the parameter values outside of the sim, 
       ## but for now just get it working.
      
      N_strains     <- 500
      ## Set up diffusion. This gives the proportion of each class that moves up and down.
      mut_rate      <- 0.25
      ## Can manipulate bias or not by adding advective velocity
      biased_mut    <- 0.40
      ## For tran.1D need to adjust mut rate so that the total loss takes into account the biased mut
      mut_rate      <- mut_rate - biased_mut/2
      
      epi.strains   <- setup.grid.1D(x.up = -4, x.down = 4, N = N_strains)
      strains.beta  <- power_tradeoff(alpha = strains.alpha, c = power_c, curv = power_exp)
      epi.length    <- 400
      
      ## Rest of the parameters
      epi.params <- c(
        mut_rate   = mut_rate
        , biased_mut = biased_mut
        , alpha      = list(strains.alpha)
        , beta       = list(strains.beta)
        , gamma      = 0.2
        , d          = 0.01
        , b          = 0.1
        , N_strains  = N_strains
        )

      epi.start <- c(
        S = c(999)
      , I = c(1, rep(0, N_strains - 1))
        )
      
      ## Run deterministic model here
      epi.out <- ode.1D(
        y          = epi.start
        , times      = seq(1, epi.length)
        , func       = epi.fun.wm
        , parms      = epi.params
        , hini       = 0.01
        , hmax       = 0.01
        , method     = "lsoda"
        , nspec      = 1
        )

  epi.out

    }
}
