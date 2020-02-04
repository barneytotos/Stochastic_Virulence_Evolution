#######
## SIS model with no host demographics. Tradeoff involves just transmission and recovery/clearance. No mortality.
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
  ## correctly. Need to calculate intrinsic beta from the tradeoff curve (especially if host evolution
   ## is happening and/or if non 0 resistance and tolerance are used)
## (3) Don't give more than one host or pathogen genotype as starting values, this is hypothetically coded
 ## but will currently run into problems establishing starting conditions. 
## (4) Host resistance and tolerance available as parameters, and can also be set to evolve, but are best set to
 ## some constant value (all results in the thesis chapter in this repo and the poster are with 0 resistance
 ## and 0 tolerance) as they currently evolve without constraint, leading to uninteresting results and possibly error
## (5) Stochastic vs Deterministic option should probably be separated. Works, but the top level sim
 ## function is a bit too long for my liking at present 

######
## Accessory Functions
######

get_mut_p        <- function (orig_trait, no_tradeoff, nt_mut_var_pos_trait, power_c, power_exp, mut_link_p, mut_mean
  , mut_sd, mut_type = "shift", agg_eff_adjust, eff_hit, parasite_tuning, tradeoff_only) {

  if (no_tradeoff) {
    
    if (nt_mut_var_pos_trait) {
       pos_trait_adj <- rnorm(length(orig_trait$pos_trait), mut_mean, mut_sd)
       new_trait_pos <- orig_trait$pos_trait + pos_trait_adj   
       ## Function always returns both new_trait_pos and new_trait_neg so need the other one
       new_trait_neg <- orig_trait$neg_trait
    } else {
       neg_trait_adj <- rnorm(length(orig_trait$neg_trait), mut_mean, mut_sd)
       new_trait_neg <- orig_trait$neg_trait + neg_trait_adj
       new_trait_pos <- orig_trait$pos_trait
    }
    
    
    
  } else {
  
   if (mut_type == "shift") {

      ## Assume a neutrally evolving aggressiveness, regardless of method of calculating efficiency
       ## recall I call virulence the parasite's "negative" trait (negative trait could also be host
       ## recovery, but that isn't prepped yet)
     if (parasite_tuning) {
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
       ## For now assume independently evolving tuning

       if (!parasite_tuning) {

           ## tradeoff_only controls whether parasites evolve only with a tradeoff curve and none of the efficiency stuff
           if (!tradeoff_only) {
               if (agg_eff_adjust) {

                   new_trait_pos <- orig_trait$pos_trait - abs(neg_trait_adj) * eff_hit
                   ## Assume a negatively evolving efficiency
                   new_trait_pos <- new_trait_pos +
                       rnorm(length(new_trait_pos),
                             mut_mean
                           , mut_sd
                             )
               } else {
                   new_trait_pos <- orig_trait$pos_trait +
                       rnorm(length(orig_trait$pos_trait),
                             mut_mean
                           , mut_sd)
          }
    ## Only a tradeoff curve. No evolution directly in beta. Beta is given by the tradeoff curve. Store old trait for now, update in the next function
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
     
   } else stop("unknown mut_type")
  }

       if (any(is.na(c(new_trait_pos, new_trait_neg)))) stop("NAs detected in trait vectors??")
       return(list(new_trait_pos = new_trait_pos, new_trait_neg = new_trait_neg))
}

##' draw new mutation values
##' @param state current system state
##' @param mut_var variable under mutation ("beta" or otherwise)
##' @param orig_trait vector of original trait values
do_mut           <- function (state, no_tradeoff, nt_mut_var_pos_trait, orig_trait, power_c, power_exp, mut_link_p, mutated, parasite_tuning, tradeoff_only, eff_scale, ...) {
  
    ## Goal here will be to always run these first lines, then have the output finished if(neutral)
  
    ## !! I clearly don't understand something about these nested function calls, why tradeoff_only now isn't found, without = tradeoff_only or why my , ... started to break
    new_trait        <- get_mut_p(orig_trait, no_tradeoff, nt_mut_var_pos_trait, power_c, power_exp, mut_link_p, tradeoff_only = tradeoff_only, parasite_tuning = parasite_tuning, ...)
    state$lpostrait  <- c(state$lpostrait, new_trait$new_trait_pos)
    state$lnegtrait  <- c(state$lnegtrait, new_trait$new_trait_neg)
    
    ## if using tradeoff curve, translate the new mutant link-scale values into the probability scale by using the tradeoff curve in 
     ## the appropriate way, determined by the model chosen (tradeoff only, efficiency, or tuning). Otherwise, for no_tradeoff this
      ## function will just use the link inverse 

    state <- update_mut_pt(
                    state           = state
                  , orig_trait      = orig_trait
                  , power_c         = power_c
                  , power_exp       = power_exp
                  , mut_link_p      = mut_link_p
                  , mutated         = mutated
                  , no_tradeoff     = no_tradeoff
                  , parasite_tuning = parasite_tuning
                  , tradeoff_only   = tradeoff_only
                  , eff_scale       = eff_scale)

    
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

##' Update mutant strain's trait values using power-law tradeoff. 
update_mut_pt    <- function (state, orig_trait, power_c, power_exp, mut_link_p, mutated, mutated_host, no_tradeoff, tradeoff_only, ...) {

  if (!no_tradeoff) {
    ## Scale beta according to tradeoff curve
    new_par_beta <- scale_beta_alpha(state, power_c, power_exp, mut_link_p, tradeoff_only, ...)

    ## Resistance will act to decrease parasite transmission and virulence following the shape of the tradeoff curve
    ## Need to think critically about what scale this should be conducted on. Both logit and probability scale feel like
    ## they have problems
    ## First calculate a tradeoff curve with the same curvature that passes through the parasite's (alpha, beta)
    cvec         <- pt_calc_c(beta = new_par_beta, alpha = mut_link_p$linkinv(state$lnegtrait), curv = power_exp)

    ## For each of these tradeoff curves, calculate a new alpha and beta for each host that is infected
    new_negtraits   <- t(mut_link_p$linkinv(state$lnegtrait))
    new_postraits    <- matrix(power_tradeoff(
        alpha = c(new_negtraits)
     ,  c     = rep(cvec, each = nrow(new_negtraits)) 
    ,  curv   = power_exp)
    , nrow    = nrow(new_negtraits), ncol = ncol(new_negtraits))

    state$neg_trait  <- new_negtraits
    state$pos_trait  <- new_postraits
    
    ## !! Placeholder for now because I realized $lpostrait never gets updated if tradeoff_only == TRUE (it isn't used)
     ## in any way, but it should be updated because it is returned in the output without warning. Will want to put this
      ## in a better spot, but want it down for now
    if (tradeoff_only) {
    state$lpostrait <- mut_link_p$linkfun(state$pos_trait)
    }
    
  } else {
    
    ## placeholder for now
    state$neg_trait <- matrix(mut_link_p$linkinv(state$lnegtrait), nrow = 1) 
    state$pos_trait <- matrix(mut_link_p$linkinv(state$lpostrait), nrow = 1)
  }

   return(state)
}

##' Remove extinct strains
do_extinct       <- function (state, extinct, parasite) {

 #  state[[mut_var]] <- state[[mut_var]][,-extinct, drop = FALSE]
    state$pos_trait  <- state$pos_trait[,-extinct, drop = FALSE]
    state$neg_trait  <- state$neg_trait[,-extinct, drop = FALSE] 
    state$lpostrait  <- state$lpostrait[-extinct, drop = FALSE]
    state$lnegtrait  <- state$lnegtrait[-extinct, drop = FALSE]
    state$Imat       <- state$Imat[,-extinct, drop = FALSE]

    return(state)
}
## This function is borderline hideous but it does what I want it to do (multinomial draws and return proper
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
    max_beta      <- power_tradeoff(c = power_c, alpha = mut_link_p$linkinv(state$lnegtrait), curv = power_exp)

    ## realized beta. If efficiency is directly the trait evolving use the "positive" trait, otherwise calculate efficiency from tuning
    if (!parasite_tuning && tradeoff_only) {
        realized_beta <- max_beta
    } else {
        if (parasite_tuning) {
            ## Gaussian fitness as a function of tuning
            realized_beta <- exp(-(state$lpostrait - state$lnegtrait)^2 / eff_scale)  * max_beta
        } else {
            ## fraction of max beta
            realized_beta <- mut_link_p$linkinv(state$lpostrait) * max_beta  
        }
    }
    return(realized_beta)
}

## Calculate starting trait values for parasite and host | desired starting R0
calc_startvals        <- function (neg_trait0, tuning, N, power_c, power_exp, mut_link_p, eff_scale, parasite_tuning, tradeoff_only, no_tradeoff, pos_trait0) {

## No resistance in SIS, but leaving structure for now...
negtrait_r     <- mut_link_p$linkinv(mut_link_p$linkfun(neg_trait0)) #- mut_link_h$linkfun(res0))

## Also no tolerance in SIS, but leaving structure for now...
negtrait_rt    <- mut_link_p$linkinv(mut_link_p$linkfun(negtrait_r)) # - mut_link_h$linkfun(tol0))

## Symmetrical matrix of efficiencies associated with combination of each parasite trait
if (!no_tradeoff) {
if (parasite_tuning) {
effic       <- outer(mut_link_p$linkfun(neg_trait0), mut_link_p$linkfun(tuning), FUN = eff_calc, eff_scale = eff_scale)
} else {
  if (!tradeoff_only) {
    effic <- matrix(data = tuning, nrow = 1, ncol = 1)
  } else {
    ## Should refer to number of starting strains
    effic <- matrix(data = 1, nrow = 1, ncol = 1)
  }
}
  
## From these calculate beta of each of the strains
joint_beta  <- sweep(effic, 2, power_tradeoff(alpha = negtrait_r, c = power_c, curv = power_exp), FUN = "*")
  
} else {
  
effic          <- 1
joint_beta     <- pos_trait0
  
}

return(list(
  intrinsic_beta = effic
, tuning         = tuning
, joint_postrait = joint_beta
, joint_negtrait = negtrait_rt
  ))

}

##' Select binomial random deviates, return a matrix with appropriate dimensions
##' @param n total number of deviates
##' @param size binomial N
##' @param prob binomial probability
##' @param nrow number of rows
##' @param ncol number of cols
## FIXME set default to ncol=1 ? eliminate nrow since determined by ncol? can we use dim(n) if it exists?
rbinom_mat  <- function (n, size, prob, nrow, ncol) {
  matrix(rbinom(n, size, prob), nrow = nrow, ncol = ncol)
}
## calculate efficiency based on the matching of tuning and aggressiveness
eff_calc    <- function (x, y, eff_scale) {
  exp(-(x - y)^2 / eff_scale)
}

######
## For ODE models see other script "hpevosim/R/funs_SIS_determ.R"
######

######
## power tradeoff and R0 functions
######
power_tradeoff <- function (alpha, c, curv) {
    c * alpha ^ (1 / curv)
}
power_R0       <- function (alpha, c, curv, gamma, N) {
    N * power_tradeoff(alpha, c, curv)  / ( alpha + gamma )
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

##' Top level sim function (wrapper)
##' @param deterministic FIXME
##' @param mut_type FIXME
##' @param mut_var  FIXME
##' @param Imat  FIXME
##' @param eff_scale  FIXME
##' @param progress FIXME
##' @export 
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

        return(epi.out)

    }
}

