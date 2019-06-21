#######
## SIS model with host death as well as host reproduction. Top level simulation function as well as all helper
## functions are contained in this script. 
## For extensions on the simplest model and other code originally written by BMB email Ben Bolker for access
#######

#######
## Cautions 
#######

## (1) Don't choose gamma as the mutational parameter right now, it won't work. Just beta.
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
## Get host and parasite mutations
get_mut_h        <- function (state, orig_trait, mut_var, mut_mean, mut_sd
  , res_mut, tol_mut, mut_host_sd_shift, mut_host_mean_shift, mut_type = "shift") {

    if (mut_type == "shift") {

      ## Direction of change different for host and pathogen and different for the two parameters (beta and gamam)

      ## First get new resistant mutatnts
        new_trait_r <- orig_trait$res_mut +
            rnorm(length(orig_trait$res_mut)
            , ifelse(mut_var == "beta"
              , mut_mean/mut_host_mean_shift
              , mut_mean*-1/mut_host_mean_shift)
            , mut_sd/mut_host_sd_shift)  ## allow the sd for the host to be a multiple of the sd for the pathogen (1 for equal)

      ## Second, retrieve resistance values for the new tolerant mutatnts
        repeated_trait_r <- rep(state$hrtraitvec, tol_mut)
        new_trait_r      <- c(new_trait_r, repeated_trait_r)

         new_trait_t <- orig_trait$tol_mut +
            rnorm(length(orig_trait$tol_mut)
            , ifelse(mut_var == "beta"
#              , mut_mean/mut_host_mean_shift
              , 0
              , mut_mean*-1/mut_host_mean_shift)
            , mut_sd/mut_host_sd_shift)  ## allow the sd for the host to be a multiple of the sd for the pathogen (1 for equal)

        repeated_trait_t <- rep(state$httraitvec, res_mut)
        new_trait_t      <- c(repeated_trait_t, new_trait_t)

        if (any(is.na(c(new_trait_r, new_trait_t)))) stop("??")
        return(list(res_trait = new_trait_r, tol_trait = new_trait_t))
    } else stop("unknown mut_type")
}
get_mut_p        <- function (orig_trait, mut_var, power_c, power_exp, mut_link_p, mut_mean
  , mut_sd, mut_type = "shift", agg_eff_adjust, parasite_tuning, tradeoff_only) {

   if (mut_type == "shift") {

      ## Assume a neutrally evolving aggressiveness, regardless of method of calculating efficiency
       ## recall I call virulence the parasite's "negative" trait (negative trait could also be host
       ## recover, but that isn't prepped yet)
     if (parasite_tuning == TRUE) {
       neg_trait_adj <- rnorm(length(orig_trait$neg_trait), 0, mut_sd)
       new_trait_neg <- orig_trait$neg_trait + neg_trait_adj
     } else {
       neg_trait_adj <- rnorm(length(orig_trait$neg_trait), mut_mean, mut_sd)
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

      new_trait_pos <- orig_trait$pos_trait - neg_trait_adj

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
do_mut           <- function (state, mut_var, mut_host, orig_trait, ...) {
    ## If the mutation is in the parasite
    if (mut_host == FALSE) {
    new_trait        <- get_mut_p(orig_trait, mut_var, ...)
    ## Whether or not efficiency is directly evolving or not (tuning in this case), or we are talking
     ## just about beta (tradeoff curve only) I just refer to it as the "positive" trait here
    state$ltraitvec  <- c(state$ltraitvec, new_trait$new_trait_pos)
    state$palphavec  <- c(state$palphavec, new_trait$new_trait_neg)
    } else {
    ## Else if the mutation is in the host (some adjustments to how state gets updated when the host is evolving)
    new_trait        <- get_mut_h(state, orig_trait, mut_var,...)
    state$hrtraitvec <- c(state$hrtraitvec, new_trait$res_trait)
    state$httraitvec <- c(state$httraitvec, new_trait$tol_trait)
    }
    return(state)
}
## Update mutant strain's trait values using power-law tradeoff. 
update_mut_pt    <- function (state, orig_trait, power_c, power_exp, mut_link_p, mut_link_h, mutated, mutated_host, mut_var, ...) {

   ## Scale beta according to tradeoff curve
   new_par_beta <- scale_beta_alpha(state, power_c, power_exp, mut_link_p, ...)

   ## Resistance will act to decrease parasite transmission and virulence following the shape of the tradeoff curve
     ## Need to think criticall about what scale this should be conducted on. Both logit and probability scale feel like
      ## they each have problems
   ## First calculate a tradeoff curve with the same curvature that passes through the parasite's (alpha, beta)
   cvec         <- pt_calc_c(beta = new_par_beta, alpha = mut_link_p$linkinv(state$palphavec), curv = power_exp)

   ## For each of these tradeoff curves, calculate a new alpha and beta for each host that is infected
  # new_alphas  <- t(outer(mut_link_p$linkinv(state$palphavec), mut_link_h$linkinv(state$hrtraitvec), "*"))
   new_alphas   <- t(mut_link_p$linkinv(outer(state$palphavec, state$hrtraitvec, "-")))
   new_betas    <- matrix(power_tradeoff(rep(cvec, each = nrow(new_alphas)), alpha = c(new_alphas)
     , curv = power_exp), nrow = nrow(new_alphas), ncol = ncol(new_alphas))

   state$alpha  <- new_alphas
   state$beta   <- new_betas

   ## ** Important order problem??? Because of the nonlinear scaling, applying resistance first and then tolerance second,
    ## tolerance will have a smaller effect [I think]. But tolerance doesn't affect beta, so unclear how to fix this.
    ## Not important for now until resistance and tolerance get constraints put on them anyway though

   ## Further adjust alpha via tolerance, which will act as a multiple to parasite intrinsic mortality rate
   state$alpha  <- get_alpha_tol(state, numcol = ncol(state$alpha), numrow = nrow(state$alpha), mut_link_h, mut_link_p)

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

   if (sum(mutated_host) > 0) {
   ## Update Svec (mutated_host is a vector that tracks which host mutated)
   state$Svec       <- cbind(state$Svec, matrix(rep(1, sum(mutated_host)), nrow = 1))
   ## Add new rows with 0s for the new S genotypes
   state$Imat       <- rbind(state$Imat, matrix(data = rep(0, sum(mutated_host)*ncol(state$Imat))
     , ncol = ncol(state$Imat), nrow = sum(mutated_host)))
   }

   return(state)
}
## Mutation updating old version (no tradeoff -- no relationship between beta and alpha)
update_mut       <- function (state, mut_link, mutated, mutated_host, mut_var) {

   state[[mut_var]] <- t(mut_link$linkinv(outer(state$ltraitvec, state$hrtraitvec,  "/")))

   ## Parasite intrinsic mortality rate following some relationship yet to be determined
   state$palphavec  <- rep(state$palphavec[1], length(state$ltraitvec))

   ## Tolerance will act as a multiple to parasite intrinsic mortality rate
   state$alpha      <- t(outer(state$palphavec, state$hrtraitvec,  "*"))

   ## Update Infected matrix with the new strain, maintaining which S class received that mutation.
    ## For each mutated strain, Imat gets a new column with a single 1, in the row in which the mutation occurred
   if (sum(mutated) > 0) {
     ## First make a matrix of 0s, then add a single one in each column corresponding to the row of the host that the parasite mutated in
     new_mutes            <- matrix(data = 0, nrow = nrow(state$Imat), ncol = sum(mutated))
     num_mutes            <- matrix(data = c(rep(c(mutated), c(mutated)), 1:sum(mutated)), ncol = 2, nrow = sum(mutated))
     new_mutes[num_mutes] <- 1
     state$Imat           <- cbind(state$Imat, new_mutes)
   }

   if (sum(mutated_host) > 0) {
   ## Update Svec (mutated_host is a vector that tracks which host mutated)
   state$Svec       <- cbind(state$Svec, matrix(rep(1, sum(mutated_host)), nrow = 1))
   ## Add new rows with 0s for the new S genotypes
   state$Imat       <- rbind(state$Imat, matrix(data = rep(0, sum(mutated_host)*ncol(state$Imat)), ncol = ncol(state$Imat), nrow = sum(mutated_host)))
   }

   return(state)
}
## Remove extinct strains
do_extinct       <- function (state, mut_var, extinct, parasite) {
  ## If a parasite strain has gone extinct remove a column
  if (parasite == TRUE) {
    state[[mut_var]] <- state[[mut_var]][,-extinct, drop = FALSE]
  ## Again, will eventually need to clean this up
    state$alpha      <- state$alpha[,-extinct, drop = FALSE]
    state$ltraitvec  <- state$ltraitvec[-extinct, drop = FALSE]
    state$palphavec  <- state$palphavec[-extinct, drop = FALSE]
    state$Imat       <- state$Imat[,-extinct, drop = FALSE]
 ## If a host strain has gone extinct remove a row
  } else {
    state[[mut_var]] <- state[[mut_var]][-extinct, , drop = FALSE]
    state$alpha      <- state$alpha[-extinct, , drop = FALSE]
    state$hrtraitvec <- state$hrtraitvec[-extinct, drop = FALSE]
    state$httraitvec <- state$httraitvec[-extinct, drop = FALSE]
    state$Imat       <- state$Imat[-extinct, , drop = FALSE]
    state$Svec       <- state$Svec[, -extinct, drop = FALSE]
  }
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
## Calculate death, either using a classic density dependence function, or a constant rate
get_death_con    <- function (state, d, S) {
  if (S == TRUE) {
  deathS     <- rbinom_mat(n = length(state$Svec), size = state$Svec, prob = d, nrow = nrow(state$Svec), ncol = ncol(state$Svec))
  state$Svec <- state$Svec - deathS
  return(list(state, deathS))
  } else {
  deathI1     <- rbinom_mat(n = length(c(state$Imat)), size = c(state$Imat), prob = d, nrow = nrow(state$Imat), ncol = ncol(state$Imat))
  state$Imat  <- state$Imat - deathI1
  deathI2     <- rbinom_mat(n = length(c(state$Imat)), size = c(state$Imat), prob = c(state$alpha), nrow = nrow(state$Imat), ncol = ncol(state$Imat))
  state$Imat  <- state$Imat - deathI2
  return(list(state, deathI1, deathI2))
  }
}
get_death_dd     <- function (S, I, N0, d0, decay) {
  d0 * exp(-(S + I) / (N0 / decay))
}
## Calculate birth, either using a classic density dependence function, or a balancing function
 ## dd = density dependence
get_birth_dd     <- function (state, N0, b0, decay) {
  birth_rate <- b0 * exp(-(sum(state$Svec) + sum(state$Imat)) / (N0 / decay))
  birth      <- rbinom_mat(n = length(state$Svec), size = state$Svec + rowSums(state$Imat) ## Assume S and I reproduce
    , prob = birth_rate, nrow = nrow(state$Svec), ncol = ncol(state$Svec))
  birth
}
 ## bal = balance with stochasticity
get_birth_bal    <- function (state, d) {

  ## Set birth rate equal to the total death rate among all infected individuals and susceptible individuals from parasite virulence and background death
 birth_rate <- (sum(state$Imat - (state$Imat * (1 - d) * (1 - state$alpha))) + (sum(state$Svec) * d)) / sum(state$Svec)
  if (birth_rate != 0 & sum(state$Svec > 0)) {
    if (birth_rate < 1) {
    birth <- rbinom_mat(n = length(state$Svec), size = state$Svec, prob = birth_rate, nrow = nrow(state$Svec), ncol = ncol(state$Svec))
    } else {
  ## If the number of deaths is bigger than the number of S, keep birth marginally stochastic by adding births equal to the difference
   ## between death and birth, but at least double the number of S by allowing each S to reproduce
    birth <- state$Svec * floor(birth_rate) + rbinom_mat(n = length(state$Svec), size = state$Svec, prob = (birth_rate - floor(birth_rate)), nrow = nrow(state$Svec), ncol = ncol(state$Svec))
    }
  } else {
    birth <- matrix(rep(0, length(state$Svec)), nrow = 1)
  }
 birth

}
 ## det = deterministic (exact birth = death)
get_birth_det    <- function (deathS, deathI) {
  ## Total birth per S class is the sum of the number of deaths in S and death in the I classes
  rowSums(deathI[[3]]) + rowSums(deathI[[2]]) + deathS[[2]]
}
 ## fill = slots, max pop stable
get_birth_fill   <- function (state, N0) {
## sum the number of each host type and multiply by investment in resistance and tolerance
 ## to get a weighted number of each host type (placeholder for now)
placeholder_res <- 1
placeholder_tol <- 1
weighted_hosts  <- (rowSums(state$Imat) + state$Svec) * placeholder_res * placeholder_tol
## sequence of 1s and 0s for occupied "slots" and non-occupied "slots"
pop_slots       <- c(rep(1, sum(weighted_hosts)), rep(0, N0 - sum(weighted_hosts)))
birth           <- apply(weighted_hosts, 2, function (x) sum(1 - sample(pop_slots, x)))
birth
}
 ## fill = slots, new version with some uncertainty but still a *very strong* form
 ## of birth to stabalize the population 
get_birth_fill   <- function (state, N0) {
## sum the number of each host type and multiply by investment in resistance and tolerance
 ## to get a weighted number of each host type (placeholder for now)
placeholder_res <- 1
placeholder_tol <- 1
weighted_hosts  <- (rowSums(state$Imat) + state$Svec) * placeholder_res * placeholder_tol
## sequence of 1s and 0s for occupied "slots" and non-occupied "slots"
#pop_slots       <- c(rep(1, sum(weighted_hosts)), rep(0, N0 - sum(weighted_hosts)))
#birth           <- apply(weighted_hosts, 2, function (x) sum(1 - sample(pop_slots, x)))
if (N0 > weighted_hosts) {
birth <- rbinom(
      n = 1
    , size = N0
    , prob = (N0 - weighted_hosts) / N0)
} else {
  birth <- 0
}
birth
}
## Function to wrap apply into returning a matrix (problems with always needing a matrix and reducing to a vector sometimes)
get_alpha_tol    <- function (state, numcol, numrow, mut_link_h, mut_link_p) {

  mut_link_p$linkinv(
  matrix(
      apply(mut_link_p$linkfun(state$alpha), 2, function (x) x - matrix(state$httraitvec, ncol = 1))
    , nrow = numrow
    , ncol = numcol
  )
  )

}
## Power law relationship between alpha and beta
power_tradeoff   <- function (alpha, c, curv) {
  c * alpha ^ (1 / curv)
}
## calculate c for a power law that goes through the current strain | efficiency
pt_calc_c        <- function (alpha, beta, curv) {
 beta / ( alpha ^ (1 / curv) )
}
## Scale parasite beta and alpha evolution by the tradeoff
scale_beta_alpha <- function (state, power_c, power_exp, mut_link_p, parasite_tuning, eff_scale, ...) {

## Determine the beta of a given pathogen strain | that pathogen's current alpha value
 ## Maximum possible beta
max_beta      <- power_tradeoff(c = power_c, alpha = mut_link_p$linkinv(state$palphavec), curv = power_exp)

## realized beta. If efficiency is directly the trait evolving use the "positive" trait, otherwise calculate efficiency from tuning
if (parasite_tuning == FALSE) {
realized_beta <- max_beta
} else {
realized_beta <- exp(-(state$ltraitvec - state$palphavec)^2 / eff_scale)  * max_beta
}

return(realized_beta)

}
## Calculate starting trait values for parasite and host | desired starting R0
calc_startvals        <- function (alpha0, tuning, res0, tol0, gamma0, d, N, power_c, power_exp, mut_link_h, mut_link_p, eff_scale, parasite_tuning) {

## alpha only considering host resistance
alpha_r     <- mut_link_p$linkinv(mut_link_p$linkfun(alpha0) - mut_link_h$linkfun(res0))

## alpha given host resistance and tolerance
alpha_rt    <- mut_link_p$linkinv(mut_link_p$linkfun(alpha_r) - mut_link_h$linkfun(tol0))

## calculate efficiency based on the matching of tuning and aggressiveness
eff_calc    <- function (x, y, eff_scale) {
  exp(-(x - y)^2 / eff_scale)
}
## Symmetrical matrix of efficiencies associated with combination of each parasite trait
effic       <- outer(mut_link_p$linkfun(alpha0), mut_link_p$linkfun(tuning), FUN = eff_calc, eff_scale = eff_scale)
if (parasite_tuning == FALSE) {
  effic <- matrix(data = 1, nrow = 1, ncol = 1)
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
calc_startvals_determ <- function (alpha0, tuning, res0, tol0, gamma0, d, N, power_c, power_exp, mut_link_h, mut_link_p, eff_scale) {

## alpha only considering host resistance
alpha_r     <- mut_link_p$linkinv(alpha0 - mut_link_h$linkfun(res0))

## alpha given host resistance and tolerance
alpha_rt    <- mut_link_p$linkinv(mut_link_p$linkfun(alpha_r) - mut_link_h$linkfun(tol0))

## calculate efficiency based on the matching of tuning and aggressiveness
eff_calc    <- function (x, y, eff_scale) {
  exp(-(x - y)^2 / eff_scale)
}
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
calc_startvals_determ_for_plotting <- function (alpha0, tuning, res0, tol0, gamma0, d, N, power_c, power_exp, mut_link_h, mut_link_p, eff_scale) {

## alpha only considering host resistance
alpha_r     <- mut_link_p$linkinv(mut_link_p$linkfun(alpha0) - mut_link_h$linkfun(res0))

## alpha given host resistance and tolerance
alpha_rt    <- mut_link_p$linkinv(mut_link_p$linkfun(alpha_r) - mut_link_h$linkfun(tol0))

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
rbinom_mat       <- function (n, size, prob, nrow, ncol) {
  matrix(rbinom(n, size, prob), nrow = nrow, ncol = ncol)
}

######
## Functions for exploration and such...
######
 ## slope of the tradeoff curve
power_tradeoff_slope     <- function (agg, c, curv) {
  c * plogis(agg) ^ (1/curv - 1) / curv
}
## Additive (rates)
power_tradeoff_R0        <- function (agg, c, curv, d, gamma, res, tol) {
  c * plogis(agg - res) ^ (1 / curv) / (plogis(agg - res - tol) + d + gamma)
}
## multiplicative (probabilities)
power_tradeoff_R0        <- function (agg, c, curv, d, gamma, res, tol) {
  c * plogis(agg - res) ^ (1 / curv) / (1-((1-plogis(agg - res - tol))*(1-d)*(1-gamma)))
}
agg_eff_eps_res_tol      <- function (eff , agg, eps, c, curv, d, gamma, adj, inc_alpha, res, tol, eff_mut, prop_change) {

if (prop_change == FALSE) {

  if (inc_alpha == TRUE) {
    ## Same as the previous function but also with resistance and tolerance, where resistance opperates to
     ## decrease aggressiveness with a correlated decrease in transmission, and where tolerance opperates to
      ## just decrease the size of the effect of aggressiveness on host mortality rate

 ( ( plogis(eff - eps*adj + eff_mut) * c * plogis(agg + eps - res) ^ (1/curv) ) / ( plogis(agg + eps - res - tol) + d + gamma) ) -
 ( ( plogis(eff                    ) * c * plogis(agg       - res) ^ (1/curv) ) / ( plogis(agg       - res - tol) + d + gamma) )


  } else {


 ( ( plogis(eff + eps*adj) * c * plogis(agg - eps - res) ^ (1/curv) ) / ( plogis(agg - eps - res - tol) + d + gamma) ) -
 ( ( plogis(eff          ) * c * plogis(agg       - res) ^ (1/curv) ) / ( plogis(agg       - res - tol) + d + gamma) )


  }

} else {

  if (inc_alpha == TRUE) {
    ## Same as the previous function but also with resistance and tolerance, where resistance opperates to
     ## decrease aggressiveness with a correlated decrease in transmission, and where tolerance opperates to
      ## just decrease the size of the effect of aggressiveness on host mortality rate

 ( ( plogis(eff - eps*adj + eff_mut) * c * plogis(agg + eps - res) ^ (1/curv) ) / ( plogis(agg + eps - res - tol) + d + gamma) ) /
 ( ( plogis(eff                    ) * c * plogis(agg       - res) ^ (1/curv) ) / ( plogis(agg       - res - tol) + d + gamma) )


  } else {


 ( ( plogis(eff + eps*adj) * c * plogis(agg - eps - res) ^ (1/curv) ) / ( plogis(agg - eps - res - tol) + d + gamma) ) /
 ( ( plogis(eff          ) * c * plogis(agg       - res) ^ (1/curv) ) / ( plogis(agg       - res - tol) + d + gamma) )


  }


}

}
agg_eff_eps_res_tol_tune <- function (tune, agg, eps, c, curv, d, gamma,      inc_alpha, res, tol, tune_mut, prop_change, eff_scaling) {

if (prop_change == FALSE) {

  if (inc_alpha == TRUE) {
    ## Same as the previous function but also with resistance and tolerance, where resistance opperates to
     ## decrease aggressiveness with a correlated decrease in transmission, and where tolerance opperates to
      ## just decrease the size of the effect of aggressiveness on host mortality rate

 ( ( (exp(-((tune + tune_mut) - (agg + eps))^2 / eff_scaling)) * c * plogis(agg + eps - res) ^ (1/curv) ) / ( plogis(agg + eps - res - tol) + d + gamma) ) -
 ( ( (exp(-((tune           ) -  agg       )^2 / eff_scaling)) * c * plogis(agg       - res) ^ (1/curv) ) / ( plogis(agg       - res - tol) + d + gamma) )


  } else {


 ( ( (exp(-((tune + tune_mut) - (agg - eps))^2 / eff_scaling)) * c * plogis(agg - eps - res) ^ (1/curv) ) / ( plogis(agg - eps - res - tol) + d + gamma) ) -
 ( ( (exp(-((tune           ) -  agg       )^2 / eff_scaling)) * c * plogis(agg       - res) ^ (1/curv) ) / ( plogis(agg       - res - tol) + d + gamma) )


  }

} else {

  if (inc_alpha == TRUE) {
    ## Same as the previous function but also with resistance and tolerance, where resistance opperates to
     ## decrease aggressiveness with a correlated decrease in transmission, and where tolerance opperates to
      ## just decrease the size of the effect of aggressiveness on host mortality rate

 ( ( (exp(-((tune + tune_mut) - (agg + eps))^2 / eff_scaling)) * c * plogis(agg + eps - res) ^ (1/curv) ) / ( plogis(agg + eps - res - tol) + d + gamma) ) /
 ( ( (exp(-((tune           ) -  agg       )^2 / eff_scaling)) * c * plogis(agg       - res) ^ (1/curv) ) / ( plogis(agg       - res - tol) + d + gamma) )


  } else {


 ( ( (exp(-((tune + tune_mut) - (agg - eps))^2 / eff_scaling)) * c * plogis(agg - eps - res) ^ (1/curv) ) / ( plogis(agg - eps - res - tol) + d + gamma) ) /
 ( ( (exp(-((tune           ) -  agg       )^2 / eff_scaling)) * c * plogis(agg       - res) ^ (1/curv) ) / ( plogis(agg       - res - tol) + d + gamma) )


  }

}

}

######
## RD functions
######
## Deterministic differential equation model with diffusion for mutation, either just parasite evo or also res tol evo
par_evo_determ          <- function (t, y, params, ...) {

  Y <- list(
   Svec = y[1]
 , Imat = matrix(data = y[-1], nrow = nrow(params$beta), ncol = nrow(params$beta)))

  n <- nrow(Y$Imat)
  
## Note that bins are equally spaced from mut_link_p$linkfun(0.001):mut_link_p$linkfun(0.999), this makes diffusion
 ## happen faster or slower depending on where you are in the logit, which is likely ok, but something
  ## to definitely keep in mind
 
  with(as.list(c(params, Y)), {

      dS <- - (Svec * d) -  ## background death
        sum(beta * Imat * sum(Svec)) + ## loss to infection
        sum(Imat * d) +
        Svec * d +
        sum(Imat * alpha) +
        sum(gamma * Imat) ## gain due to recovery

      dI <- - (Imat * d) - 
        Imat * alpha +      ## background death and virulence death
        (beta * Imat * sum(Svec)) -  ## gain due to infection
        gamma * Imat +  ## loss due to recovery
        tran.2D(C = Imat, D.x = mutlev, dx = 1, dy = 1)$dC

    list(c(dS, dI))
  }
  )

}
## Not completed and very slow because bins > # by quite a lot
par_evo_determ_res_tol  <- function (t, y, params, ...) {

  Y <- list(
    Smat = matrix(data = y[1:nrow(params$res)*2], nrow = nrow(params$res), ncol = nrow(params$res))
  , Imat = matrix(data = y[(nrow(params$res)*2 + 1):(nrow(params$res)*2 + 1 + nrow(params$beta)*2)]
    , nrow = nrow(params$beta), ncol = nrow(params$beta))
   )

  n <- nrow(Y$Imat)

  with(as.list(c(params, Y)), {

      dS <- - (Smat * d) -  ## background death
        sum(beta * Imat * sum(Svec)) + ## loss to infection
        sum(Imat * d) +
        Smat * d +
        sum(sweep(Imat, 2, alpha, FUN = "*")) +
        sum(gamma * Imat) ## gain due to recovery

      dI <- - (Imat * d) - sweep(Imat, 2, alpha, FUN = "*") +  ## background death and virulence death
        (beta * Imat * sum(Svec)) -  ## gain due to infection
        gamma * Imat +  ## loss due to recovery
        tran.2D(C = Imat, D.x = mutlev, dx = 1, dy = 1)$dC

    list(c(dS, dI))
  }
  )

}

######
## AD functions
######
## Two types of Adaptive Dynamics: 
 ## (1) A crude form that relies on R0 theory, but doesn't actually solve any epidemiological model
  ## like the other Reaction Diffusion model, movement in bins in mut_link_p$linkfun(0.001):mut_link_p$linkfun(0.999) space 
par_evo_AD       <- function (c, curv, eff_scale, mut_link, numbins, Iseed, simul_mut, max_range = FALSE, debug1 = FALSE, debug1_val = NULL) {
  alpha_range  <- c(seq(
        mut_link_p$linkfun(0.01), mut_link_p$linkfun(0.99)
        , length = numbins))
  tuning_range <- c(seq(
        mut_link_p$linkfun(0.01), mut_link_p$linkfun(0.99)
        , length = numbins))
  
  eff_calc     <- function (x, y, eff_scale) {
   exp(-(x - y)^2 / eff_scale)
  }
  
  ## Calculate efficiency, and use efficiency to calculate beta surface and R0 surface
  effic        <- outer(alpha_range, tuning_range, FUN = eff_calc, eff_scale = eff_scale)
  
  beta_surface <- sweep(effic, 2, power_tradeoff(
    alpha = mut_link$linkinv(alpha_range)
  , c     = c
  , curv  = curv), FUN = "*")
  
  R0_surface   <- sweep(beta_surface, 2, (mut_link$linkinv(alpha_range) + 0.2 + 0.01), FUN = "/")
  
  ## seed the trait-space matrix
  I_matrix     <- matrix(data = 0, nrow = dim(R0_surface)[1], ncol = dim(R0_surface)[2])
  I_matrix[Iseed[1], Iseed[2]] <- 1
  
## move up the surface in an AD sort of way. Check which direction leads to the largest difference in R0 and move
 ## in that direction
  
  ## counter
  step_count <- 1
  I_coord    <- data.frame(
    alpha     = alpha_range[Iseed[2]]
  , alpha_loc = Iseed[2]
  , tune      = tuning_range[Iseed[1]]
  , tune_loc  = Iseed[1]
  , time      = step_count
  , beta      = beta_surface[Iseed[1], Iseed[2]]
  , R0        = R0_surface[Iseed[1], Iseed[2]])

  ## placeholder to initiate the loop
  trait_move <- c(1, 1)
  
  ## Store gradient at each step for debugging purposes
  grad_check <- data.frame(
    left   = numeric(1)
  , down   = numeric(1)
  , diagDL = numeric(1)
  , diagDR = numeric(1)
  , diagUL = numeric(1)
  )
  
  ## Do gradient ascent (FALSE) or find the maximum possible window in which R0 can icnrease (TRUE)?
  if (max_range == FALSE) {
  
  ## step 2: loop until there is no longer a larger R0
  while (trait_move[1] != 0 | trait_move[2] != 0) {
    
  step_count <- step_count + 1
    
  ## step 1: check in what direction there is a larger R0, only in the immediate vicinity (small mutation step)
  trait_row <- which(rowSums(I_matrix) == 1)
  trait_col <- which(colSums(I_matrix) == 1)
  R0_vicinity <- R0_surface[
    seq(trait_row - 1,  trait_row + 1, by = 1)
  , seq(trait_col - 1,  trait_col + 1, by = 1)
    ]
  ## Is there a mutation in both directions allowed?
  if (simul_mut == FALSE) {
  R0_vicinity[1, 1] <- 0; R0_vicinity[1, 3] <- 0; R0_vicinity[3, 1] <- 0; R0_vicinity[3, 3] <- 0
  }
  R0_near_max <- R0_vicinity / R0_surface[trait_row, trait_col]
  trait_move  <- which(R0_near_max == max(R0_near_max), arr.ind = TRUE)
  
  grad_check <- rbind(grad_check
   , data.frame(
     left = R0_near_max[2, 1]
   , down = R0_near_max[3, 2]
   , diagDL = R0_near_max[3, 1]
   , diagDR = R0_near_max[3, 3]
   , diagUL = R0_near_max[1, 1]
   ))
  
  trait_move  <- trait_move - 2
  
  ## update movement. In AD new mutant completely takes over
  I_matrix[trait_row, trait_col] <- 0
  I_matrix[trait_row + trait_move[1], trait_col + trait_move[2]] <- 1
  
  ## store where the system moved
  I_coord <- rbind(
    I_coord
  , data.frame(
     alpha     = alpha_range[trait_col + trait_move[2]]
   , alpha_loc = trait_col + trait_move[2]
   , tune      = tuning_range[trait_row + trait_move[1]]
   , tune_loc  = trait_row + trait_move[1]
   , time      = step_count
   , beta      = beta_surface[trait_row + trait_move[1], trait_col + trait_move[2]]
   , R0        = R0_surface[trait_row + trait_move[1], trait_col + trait_move[2]]
  )
  )
  
  if (step_count %% 20 == 0) { print(step_count) }
  if (debug1 == TRUE) {
  if (step_count == debug1_val) { browser() }
  }
  
  }
    
  } else {
    
  ## check maximum movement in one direction followed by maximum movement in the other direction
    ## Note: Matrix of values and plot are flipped, so be careful.
    ## Matrix is: agg is rows, tuning is cols, so asking which row of column 2 is the largest is a move in agg space,
     ## despite it being movement in the plot in row space because of how I am choosing to plot
  I_coord_agg <- transform(I_coord, priority_move = "Agg")  
  while (trait_move[1] != 2 | trait_move[2] != 2) {
    
  step_count <- step_count + 1
    
  ## step 1: check in what direction there is a larger R0, only in the immediate vicinity (small mutation step)
  trait_row <- which(rowSums(I_matrix) == 1)
  trait_col <- which(colSums(I_matrix) == 1)
  R0_vicinity <- R0_surface[
    seq(trait_row - 1,  trait_row + 1, by = 1)
  , seq(trait_col - 1,  trait_col + 1, by = 1)
    ]
  ## Is there a mutation in both directions allowed?
  if (simul_mut == FALSE) {
  R0_vicinity[1, 1] <- 0; R0_vicinity[1, 3] <- 0; R0_vicinity[3, 1] <- 0; R0_vicinity[3, 3] <- 0
  }
  R0_near_max <- R0_vicinity / R0_surface[trait_row, trait_col]
  
  ## Force a movement in agg if it is available, and only then move in tuning space
  if (R0_near_max[1, 2] > 1 | R0_near_max[3, 2] > 1) {
  trait_move  <- c(which(R0_near_max[, 2] == max(R0_near_max[, 2]), arr.ind = TRUE), 2)
  trait_move  <- trait_move - 2    
  } else if (R0_near_max[2, 1] > 1 | R0_near_max[2, 3] > 1) {
  trait_move  <- c(2, which(R0_near_max[2, ] == max(R0_near_max[2, ]), arr.ind = TRUE))
  trait_move  <- trait_move - 2      
  } else {
  trait_move <- c(2, 2)
  }
  
  ## update movement. In AD new mutant completely takes over
  I_matrix[trait_row, trait_col] <- 0
  I_matrix[trait_row + trait_move[1], trait_col + trait_move[2]] <- 1
  
  ## store where the system moved
  I_coord_agg <- rbind(
    I_coord_agg
  , data.frame(
     alpha     = alpha_range[trait_row + trait_move[2]]
   , alpha_loc = trait_row + trait_move[2]
   , tune      = tuning_range[trait_col + trait_move[1]]
   , tune_loc  = trait_col + trait_move[1]
   , time      = step_count
   , beta      = beta_surface[trait_row + trait_move[1], trait_col + trait_move[2]]
   , R0        = R0_surface[trait_row + trait_move[1], trait_col + trait_move[2]]
   , priority_move = "Agg"
  )
  )
  
    if (step_count %% 20 == 0) { 
      print(step_count) 
 #     assign("I_coord_agg", I_coord_agg, envir = .GlobalEnv)
      }
  
#  if (step_count == 1000) { browser() }
    
  }
  
  ## reseed the trait-space matrix
  I_matrix     <- matrix(data = 0, nrow = dim(R0_surface)[1], ncol = dim(R0_surface)[2])
  I_matrix[Iseed[1], Iseed[2]] <- 1
  
  I_coord_tune  <- transform(I_coord, priority_move = "Tune") 
  trait_move <- c(1, 1)
  step_count <- 1
  while (trait_move[1] != 2 | trait_move[2] != 2) {
    
   step_count <- step_count + 1
    
  ## step 1: check in what direction there is a larger R0, only in the immediate vicinity (small mutation step)
  trait_row <- which(rowSums(I_matrix) == 1)
  trait_col <- which(colSums(I_matrix) == 1)
  R0_vicinity <- R0_surface[
    seq(trait_row - 1,  trait_row + 1, by = 1)
  , seq(trait_col - 1,  trait_col + 1, by = 1)
    ]
  ## Is there a mutation in both directions allowed?
  if (simul_mut == FALSE) {
  R0_vicinity[1, 1] <- 0; R0_vicinity[1, 3] <- 0; R0_vicinity[3, 1] <- 0; R0_vicinity[3, 3] <- 0
  }
  R0_near_max <- R0_vicinity / R0_surface[trait_row, trait_col]
  
  ## Force a movement in tuning if it is available, and only then move in agg space
  if (R0_near_max[2, 1] > 1 | R0_near_max[2, 3] > 1) {
  trait_move  <- c(2, which(R0_near_max[2, ] == max(R0_near_max[2, ]), arr.ind = TRUE))
  trait_move  <- trait_move - 2   
  } else if (R0_near_max[1, 2] > 1 | R0_near_max[3, 2] > 1) {
  trait_move  <- c(which(R0_near_max[, 2] == max(R0_near_max[, 2]), arr.ind = TRUE), 2)
  trait_move  <- trait_move - 2     
  } else {
  trait_move <- c(2, 2)
  }
  
  ## update movement. In AD new mutant completely takes over
  I_matrix[trait_row, trait_col] <- 0
  I_matrix[trait_row + trait_move[1], trait_col + trait_move[2]] <- 1
  
  ## store where the system moved
  I_coord_tune <- rbind(
    I_coord_tune
  , data.frame(
     alpha     = alpha_range[trait_row + trait_move[2]]
   , alpha_loc = trait_row + trait_move[2]
   , tune      = tuning_range[trait_col + trait_move[1]]
   , tune_loc  = trait_col + trait_move[1]
   , time      = step_count
   , beta      = beta_surface[trait_row + trait_move[1], trait_col + trait_move[2]]
   , R0        = R0_surface[trait_row + trait_move[1], trait_col + trait_move[2]]
   , priority_move = "Tune"
  )
  ) 
  
    if (step_count %% 20 == 0) { print(step_count) }
    
  }
    
  I_coord <- rbind(I_coord_tune, I_coord_agg)
  
  }
  
  if (max_range == FALSE) {
  return(list(I_coord, grad_check))
  } else {
  return(I_coord)    
  }
  
}
 ## (2) Another crude form that relies on R0 theory, but doesn't actually solve any epidemiological model
  ## this differs in that there is randomness, any mutant strain with eps can invade (e.g + tune, - agg or + tune, + agg if
   ## each of these options leads to higher R0)
par_evo_AD_rand  <- function (c, curv, eff_scale, mut_link, numbins, Iseed, simul_mut, max_range = FALSE, debug1 = FALSE, debug1_val = NULL) {
  alpha_range  <- c(seq(
        mut_link_p$linkfun(0.01), mut_link_p$linkfun(0.99)
        , length = numbins))
  tuning_range <- c(seq(
        mut_link_p$linkfun(0.01), mut_link_p$linkfun(0.99)
        , length = numbins))
  
  eff_calc     <- function (x, y, eff_scale) {
   exp(-(x - y)^2 / eff_scale)
  }
  
  ## Calculate efficiency, and use efficiency to calculate beta surface and R0 surface
  effic        <- outer(alpha_range, tuning_range, FUN = eff_calc, eff_scale = eff_scale)
  
  beta_surface <- sweep(effic, 2, power_tradeoff(
    alpha = mut_link$linkinv(alpha_range)
  , c     = c
  , curv  = curv), FUN = "*")
  
  R0_surface   <- sweep(beta_surface, 2, (mut_link$linkinv(alpha_range) + 0.2 + 0.01), FUN = "/")
  
  ## seed the trait-space matrix
  I_matrix     <- matrix(data = 0, nrow = dim(R0_surface)[1], ncol = dim(R0_surface)[2])
  I_matrix[Iseed[1], Iseed[2]] <- 1
  
## move up the surface in an AD sort of way. Check which direction leads to the largest difference in R0 and move
 ## in that direction
  
  ## counter
  step_count <- 1
  I_coord    <- data.frame(
    alpha     = alpha_range[Iseed[2]]
  , alpha_loc = Iseed[2]
  , tune      = tuning_range[Iseed[1]]
  , tune_loc  = Iseed[1]
  , time      = step_count
  , beta      = beta_surface[Iseed[1], Iseed[2]]
  , R0        = R0_surface[Iseed[1], Iseed[2]])

  ## placeholder to initiate the loop
  trait_move <- c(1, 1)
  
  ## Store gradient at each step for debugging purposes
  grad_check <- data.frame(
    left   = numeric(1)
  , down   = numeric(1)
  , right  = numeric(1)
  , up     = numeric(1)
  , diagDL = numeric(1)
  , diagDR = numeric(1)
  , diagUL = numeric(1)
  )
  
  ## Do gradient ascent (FALSE) or find the maximum possible window in which R0 can icnrease (TRUE)?
  if (max_range == FALSE) {
  
  end_run <- FALSE
  ## step 2: loop until there is no longer a larger R0
  while (trait_move[1] != 0 | trait_move[2] != 0) {
    
  step_count <- step_count + 1
    
  ## step 1: check in what direction there is a larger R0, only in the immediate vicinity (small mutation step)
  trait_row <- which(rowSums(I_matrix) == 1)
  trait_col <- which(colSums(I_matrix) == 1)
  R0_vicinity <- R0_surface[
    seq(trait_row - 1,  trait_row + 1, by = 1)
  , seq(trait_col - 1,  trait_col + 1, by = 1)
    ]
  ## Is there a mutation in both directions allowed?
  if (simul_mut == FALSE) {
  R0_vicinity[1, 1] <- 0; R0_vicinity[1, 3] <- 0; R0_vicinity[3, 1] <- 0; R0_vicinity[3, 3] <- 0
  }
  R0_near_max <- R0_vicinity / R0_surface[trait_row, trait_col]
  
  ## I think this is the only difference needed in the short term to very casually test
    ## randomness in AD for thesis defense
  which_move_g1 <- which(R0_near_max > 1, arr.ind = TRUE)
  if (length(which_move_g1) == 0) {
    trait_move <- c(0,0)
    end_run <- TRUE
  } else {
    trait_move <- which_move_g1[sample(nrow(which_move_g1))[1], ] 
  }

  
  grad_check <- rbind(grad_check
   , data.frame(
     left   = R0_near_max[2, 1]
   , down   = R0_near_max[3, 2]
   , right  = R0_near_max[2, 3]
   , up     = R0_near_max[1, 2]
   , diagDL = R0_near_max[3, 1]
   , diagDR = R0_near_max[3, 3]
   , diagUL = R0_near_max[1, 1]
   ))
  
  trait_move  <- trait_move - 2
  
  ## update movement. In AD new mutant completely takes over
  I_matrix[trait_row, trait_col] <- 0
  I_matrix[trait_row + trait_move[1], trait_col + trait_move[2]] <- 1
  
  ## store where the system moved
  I_coord <- rbind(
    I_coord
  , data.frame(
     alpha     = alpha_range[trait_col + trait_move[2]]
   , alpha_loc = trait_col + trait_move[2]
   , tune      = tuning_range[trait_row + trait_move[1]]
   , tune_loc  = trait_row + trait_move[1]
   , time      = step_count
   , beta      = beta_surface[trait_row + trait_move[1], trait_col + trait_move[2]]
   , R0        = R0_surface[trait_row + trait_move[1], trait_col + trait_move[2]]
  )
  )
  
  if (step_count %% 20 == 0) { print(step_count) }
  if (debug1 == TRUE) {
  if (step_count == debug1_val) { browser() }
  }
  
  if (end_run == TRUE) {
    trait_move[1] <- 0; trait_move[2] <- 0
  }
  
  }
    
  } else {
    
  ## check maximum movement in one direction followed by maximum movement in the other direction
    ## Note: Matrix of values and plot are flipped, so be careful.
    ## Matrix is: agg is rows, tuning is cols, so asking which row of column 2 is the largest is a move in agg space,
     ## despite it being movement in the plot in row space because of how I am choosing to plot
  I_coord_agg <- transform(I_coord, priority_move = "Agg")  
  while (trait_move[1] != 2 | trait_move[2] != 2) {
    
  step_count <- step_count + 1
    
  ## step 1: check in what direction there is a larger R0, only in the immediate vicinity (small mutation step)
  trait_row <- which(rowSums(I_matrix) == 1)
  trait_col <- which(colSums(I_matrix) == 1)
  R0_vicinity <- R0_surface[
    seq(trait_row - 1,  trait_row + 1, by = 1)
  , seq(trait_col - 1,  trait_col + 1, by = 1)
    ]
  ## Is there a mutation in both directions allowed?
  if (simul_mut == FALSE) {
  R0_vicinity[1, 1] <- 0; R0_vicinity[1, 3] <- 0; R0_vicinity[3, 1] <- 0; R0_vicinity[3, 3] <- 0
  }
  R0_near_max <- R0_vicinity / R0_surface[trait_row, trait_col]
  
  ## Force a movement in agg if it is available, and only then move in tuning space
  if (R0_near_max[1, 2] > 1 | R0_near_max[3, 2] > 1) {
  trait_move  <- c(which(R0_near_max[, 2] == max(R0_near_max[, 2]), arr.ind = TRUE), 2)
  trait_move  <- trait_move - 2    
  } else if (R0_near_max[2, 1] > 1 | R0_near_max[2, 3] > 1) {
  trait_move  <- c(2, which(R0_near_max[2, ] == max(R0_near_max[2, ]), arr.ind = TRUE))
  trait_move  <- trait_move - 2      
  } else {
  trait_move <- c(2, 2)
  }
  
  ## update movement. In AD new mutant completely takes over
  I_matrix[trait_row, trait_col] <- 0
  I_matrix[trait_row + trait_move[1], trait_col + trait_move[2]] <- 1
  
  ## store where the system moved
  I_coord_agg <- rbind(
    I_coord_agg
  , data.frame(
     alpha     = alpha_range[trait_row + trait_move[2]]
   , alpha_loc = trait_row + trait_move[2]
   , tune      = tuning_range[trait_col + trait_move[1]]
   , tune_loc  = trait_col + trait_move[1]
   , time      = step_count
   , beta      = beta_surface[trait_row + trait_move[1], trait_col + trait_move[2]]
   , R0        = R0_surface[trait_row + trait_move[1], trait_col + trait_move[2]]
   , priority_move = "Agg"
  )
  )
  
    if (step_count %% 20 == 0) { 
      print(step_count) 
 #     assign("I_coord_agg", I_coord_agg, envir = .GlobalEnv)
      }
  
#  if (step_count == 1000) { browser() }
    
  }
  
  ## reseed the trait-space matrix
  I_matrix     <- matrix(data = 0, nrow = dim(R0_surface)[1], ncol = dim(R0_surface)[2])
  I_matrix[Iseed[1], Iseed[2]] <- 1
  
  I_coord_tune  <- transform(I_coord, priority_move = "Tune") 
  trait_move <- c(1, 1)
  step_count <- 1
  while (trait_move[1] != 2 | trait_move[2] != 2) {
    
   step_count <- step_count + 1
    
  ## step 1: check in what direction there is a larger R0, only in the immediate vicinity (small mutation step)
  trait_row <- which(rowSums(I_matrix) == 1)
  trait_col <- which(colSums(I_matrix) == 1)
  R0_vicinity <- R0_surface[
    seq(trait_row - 1,  trait_row + 1, by = 1)
  , seq(trait_col - 1,  trait_col + 1, by = 1)
    ]
  ## Is there a mutation in both directions allowed?
  if (simul_mut == FALSE) {
  R0_vicinity[1, 1] <- 0; R0_vicinity[1, 3] <- 0; R0_vicinity[3, 1] <- 0; R0_vicinity[3, 3] <- 0
  }
  R0_near_max <- R0_vicinity / R0_surface[trait_row, trait_col]
  
  ## Force a movement in tuning if it is available, and only then move in agg space
  if (R0_near_max[2, 1] > 1 | R0_near_max[2, 3] > 1) {
  trait_move  <- c(2, which(R0_near_max[2, ] == max(R0_near_max[2, ]), arr.ind = TRUE))
  trait_move  <- trait_move - 2   
  } else if (R0_near_max[1, 2] > 1 | R0_near_max[3, 2] > 1) {
  trait_move  <- c(which(R0_near_max[, 2] == max(R0_near_max[, 2]), arr.ind = TRUE), 2)
  trait_move  <- trait_move - 2     
  } else {
  trait_move <- c(2, 2)
  }
  
  ## update movement. In AD new mutant completely takes over
  I_matrix[trait_row, trait_col] <- 0
  I_matrix[trait_row + trait_move[1], trait_col + trait_move[2]] <- 1
  
  ## store where the system moved
  I_coord_tune <- rbind(
    I_coord_tune
  , data.frame(
     alpha     = alpha_range[trait_row + trait_move[2]]
   , alpha_loc = trait_row + trait_move[2]
   , tune      = tuning_range[trait_col + trait_move[1]]
   , tune_loc  = trait_col + trait_move[1]
   , time      = step_count
   , beta      = beta_surface[trait_row + trait_move[1], trait_col + trait_move[2]]
   , R0        = R0_surface[trait_row + trait_move[1], trait_col + trait_move[2]]
   , priority_move = "Tune"
  )
  ) 
  
    if (step_count %% 20 == 0) { print(step_count) }
    
  }
    
  I_coord <- rbind(I_coord_tune, I_coord_agg)
  
  }
  
  if (max_range == FALSE) {
  return(list(I_coord, grad_check))
  } else {
  return(I_coord)    
  }
  
}

######
## Top level sim function (wrapper)
######
run_sim <- function(
    ## These first parameters are all about what type of simulation to run
     ## Note: hosts are always set to evolve. Crudly can force them never to evolve by just setting the probability to effectively 0
   deterministic       = FALSE   ## Run an advection diffusion version?
 , parasite_tuning     = TRUE    ## Reformulation where parasite efficiency is defined as matching tuning and aggressiveness
 , tradeoff_only       = FALSE   ## Stepping back to ignore tuning. Parasite just evolving according to the tradeoff curve
 , agg_eff_adjust      = FALSE   ## For efficiency model but not tuning model. Does an increase in parasite aggressiveness decrease efficiency (as a cost)
 , host_dyn_only       = FALSE   ## TRUE means no parasites, just to check for host dynamics with birth and death
 , host_evo_delay      = FALSE   ## Should hosts evolve with some delay?
 , host_evo_delay_start= NULL    ## Parameter to control the window of this delay. See function comments
 , host_evo_delay_stop = NULL    ## Parameter to control the window of this delay. See function comments
    ## These next parametes all control how the simulation is run
 , R0_init             = 2       ## >1, not actually used for tuning model
 , gamma0              = 1/5     ## >0, constant, not evolving in this version
 , gamma_max           = Inf     ## Not actually used here, plan to use when merging with Ben's code
 , alpha0              = 0.03    ## Intrinsic parasite mortality probability (without influence of hosts). Default to main thesis result start
 , tune0               = 0.97    ## Strating tuning. Give a value, but only used if parasite_tuning == TRUE. Default to main thesis result start
 , d                   = 0.01    ## Background death
 , b                   = 0.1     ## Birth rate, only used in dd birth model (give a value regardless, could be NULL)
 , b_decay             = 2.3     ## Rate of birth rate decay, only used in dd birth model (give a value regardless, could be NULL)
 , N                   = 50      ## Host population size, integer > 0
 , mu                  = 0.01    ## Mutation probability, > 0
 , mut_type            = "shift" ## Type of mutation supported. Only shift viable in this version of the code
 , mut_mean            = -1      ## Mutation mean, < 0 (for sensibility). Ignored for parasite_tuning == TRUE
 , mut_sd              = 0.5     ## Mutation sd, > 0
 , mut_var             = "beta"  ## Trait for parasite evolution. Only beta allowed for this code iteration
 , mut_link_p          = NULL    ## Default of logit scale setup in the function
 , mut_link_h          = NULL    ## Default of logit scale setup in the function
 , mut_host_sd_shift   = 1       ## **1 for identical sd to the parasite
 , mut_host_mean_shift = 1       ## **1 for identical mean to the parasite
 , mut_host_mu_shift   = 2       ## Proportion less frequent host mutation is than parasite mutation, only used if hosts are evolving 
 , mut_host_res_bias   = 0.5     ## 0.5 means equal probability of evolving resistance or tolerance, larger proportion favors more resistance evolution
 , res0                = 1       ## Starting host mean resistance value
 , res0_sd             = 0       ## Variation in resistance among starting host strains (if > 1 host strain). Not yet used, but plan to
 , tol0                = 1       ## Starting host mean tolerance value
 , tol0_sd             = 0       ## Variation in tolerance among starting host strains (if > 1 host strain). Not yet used, but plan to
 , power_c             = 0.75    ## Power law tradeoff scaling
 , power_exp           = 2       ## Power law tradeoff exponent
 , Imat                = NULL    ## Setup within function in this version, don't adjust
 , birth_type          = "dd"    ## Type of host birth, include dd (density dependent); bal (balance); det (deterministic balance); fill (strong matching to death with some stochasticity)
                                  ## See helper functions for details on these options
 , eff_scale           = 50      ## Weighting of the matching of parasite tuning and aggressiveness
 , nt                  = 100000  ## Length of simulation (time steps)
 , rptfreq             = max(nt / 500, 1) ## How often the state of the system is saved
 , seed                = NULL    ## Can set seed if desired
 , progress            = FALSE   ## Progress bar?
 , debug               = FALSE   ## A debug option for stopping at various points in the sim. See code to check where.
 , debug2              = FALSE   ## A debug option for stopping at various points in the sim. See code to check where.
 , debug3              = FALSE   ## A debug option for stopping at various points in the sim. See code to check where.
    ## A few deterministic parameters
 , determ_length       = 200     ## length of time to run the deterministic model
 , determ_timestep     = 5       ## Lsoda parameter for RD model
 , lsoda_hini          = NULL    ## Lsoda parameter for RD model
 , Imat_seed           = NULL    ## Set up inside function, keep as null
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
    , gamma0  > 0
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

    if (is.null(Imat)) {
        ## start at equilibrium I ...
        Imat <- max(1, round(N*(1 - 1/R0_init)))
        Imat <- as.matrix(Imat)
    }

    Svec  <- N-Imat

   ## Set up mutation link
    if (mut_var=="beta") {
        if (is.null(mut_link_p)) mut_link_p <- make.link("logit")
        if (is.null(mut_link_h)) mut_link_h <- make.link("log")
      } else {
      ## Positive with log link can push gamma overboard when mutation in gamma is on average disadvantageous.
       ## Made to logit link for now, not completely convinced...
        if (is.null(mut_link_p)) mut_link_p <- make.link("logit")
        if (is.null(mut_link_h)) mut_link_h <- make.link("log")
      }

    ## If deterministic == TRUE, set up the full array of possible strains from the start (for diffusion via mutation)
    if (deterministic == FALSE) {

    ## Determine intial alpha from parasite and host initial traits
     ## Note: Naming conventions are a bit odd here. alpha0 as a parameter in the function refers to the desired
      ## starting parasite mortality rate, irrespective of the host trait. That parameter is back transformed
       ## to the intrinsic parasite scale and used to calculate a true starting alpha (also called alpha0 here
        ## that is based on both parasite and host triats)

   ## Based on my new setup it seems better to directly calculate efficiency and starting beta as a function of defined tuning
    ## and aggressiveness starting values instead of going backwards from R0
    startvals  <- calc_startvals(         
      alpha0, tune0, res0, tol0, gamma0, d, N, power_c
    , power_exp, mut_link_h, mut_link_p, eff_scale, parasite_tuning)
    beta0      <- startvals$joint_beta

    } else {

      alpha0     <- c(seq(
        mut_link_p$linkfun(0.01), mut_link_p$linkfun(0.99)
        , length = 100))
      tuning     <- c(seq(
        mut_link_p$linkfun(0.01), mut_link_p$linkfun(0.99)
        , length = 100))
      startvals  <- calc_startvals_determ(
         alpha0, tuning, res0, tol0, gamma0, d, N, power_c
       , power_exp, mut_link_h, mut_link_p, eff_scale
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
    if (parasite_tuning == FALSE) {
    ltraitvec  <- mut_link_p$linkfun(startvals$joint_beta)
    } else {
      if (deterministic == FALSE) {
    ltraitvec  <- mut_link_p$linkfun(startvals$tuning)
      } else {
    ltraitvec  <- startvals$tuning        
      }
    }

    ## Alpha (intrinsic parasite mortality pressure)
    if (deterministic == FALSE) {
    palphavec  <- mut_link_p$linkfun(alpha0)
    } else {
    palphavec  <- alpha0 
    }

    ## Initial trait vectors for the host genotypes. Assumes all hosts start with identical trait value (for now)
    hrtraitvec <- rep(mut_link_h$linkfun(res0), length(res0))
    httraitvec <- rep(mut_link_h$linkfun(tol0), length(tol0))

    ## Alpha0 now becomes starting host mortality rate as a function of both host and parasite gentoype
    alpha0     <-  startvals$joint_alpha

    ## Parameters structure (parallel vectors), so these can
     ## Be modified via function and passed back ...
    state      <- list(
      beta       = as.matrix(beta0)
    , gamma      = as.matrix(gamma0) ## ...as will gamma (make matrix after, because ltraitvec needs to remain a vector)
    , alpha      = as.matrix(alpha0)
    , ltraitvec  = ltraitvec
    , palphavec  = palphavec
    , hrtraitvec = hrtraitvec
    , httraitvec = httraitvec
    , Imat       = Imat
    , Svec       = Svec)

    ## For no infection (to check host dynamics in the absence of infection)
     ## Maybe a _very_ inefficient/hacky place/way to do this, but I think makes the least clutter
    if (host_dyn_only == TRUE) {
      state$Imat[1, 1] <- 0
      state$Svec[1, 1] <- N
    }

    dfun("init")

    nrpt <- nt %/% rptfreq

     ## Added tracking host responses
    res <- as.data.frame(matrix(
      NA, nrow = nrpt, ncol = 19
    , dimnames = list(
      NULL
    , c("time"
      , "num_S"
      , "num_S_strains"
      , "num_I"
      , "num_I_strains"
      , "pop_size"
      , paste0(c("mean_hl", "sd_hl"), "res")
      , paste0(c("mean_hl", "sd_hl"), "tol")
      , paste0(c("mean_pl", "sd_pl"), mut_var)
      , paste0(c("mean_pl", "sd_pl"), "alpha")
      , "mean_alpha"
      , "sd_alpha"
      , "mean_beta"
      , "sd_beta"
      ## just list the non-evolving param
      , ifelse(mut_var == "beta", "gamma", "beta") 
      ))))

    ## Run the stochastic version
    if (deterministic == FALSE) {

    for (i in 1:nrpt) {

      ## First few if statements are: Code to check my sanity of patterns by manually starting and stoping host
       ## evolution and defined points in time. Non-dynamic and not really part of the model, mostly debugging stuff
      if (host_evo_delay == TRUE & i == host_evo_delay_start) {
      mut_host_mu_shift <- 10
      }

      if (host_evo_delay == TRUE & i >= host_evo_delay_stop) {
      state$httraitvec <- state$httraitvec / 1.05
      mut_host_mu_shift <- 1000000000
      }

            for (j in 1:rptfreq) {

               ## Debug code chunck that can be cut and paste to wherever there is a problem
                ## Fill i and j in manually after checking output after an error is encountered
                if (debug2 == TRUE) {
                  print(paste(i, j, sep = "  -  "))
                  print(str(state$Svec))
                  assign("state_check", state, .GlobalEnv)
                  if (i == 20 & j == 20) browser()
                }

                ## [Step 1.1]: Birth. Depending on type of birth birth occurs here or after a few steps (either way
                 ## new births are not accessible to death or infection until the next time step)
               if (birth_type == "dd") {
                birth <- get_birth_dd(state, N0 = N, b0 = b, decay = b_decay)
               } else if (birth_type == "bal") {
                birth <- get_birth_bal(state, d) 
               }
              
                ## [Step 2]: Death of S. Returning a list of updated state and death, so death can be used to calculate deterministic birth
                deathS     <- get_death_con(state, d, S = TRUE)
                state      <- deathS[[1]]

                ## [Step 3]: Infection.
                ## Prob of escaping infection completely
                uninf      <- rbinom_mat(n = length(state$Svec), size = state$Svec, prob = prod((1-state$beta)^state$Imat), nrow = nrow(state$Svec), ncol = ncol(state$Svec))
                ## Division of new infectives among strains
                 ## 'prob' is internally normalized
                newinf     <- get_inf(state$Svec, uninf, state$Imat, beta = state$beta)
                stopifnot(sum(rowSums(newinf) + uninf) == sum(state$Svec))

                ## [Step 4]: Death of I.
                ## Natural death + parasite induced death
                deathI     <- get_death_con(state, d, S = FALSE)
                state      <- deathI[[1]]

                ## [Step 1.2]: Other birth types
                if (birth_type == "det") {
                birth      <- get_birth_det(deathS, deathI)
                } else if (birth_type == "fill") {
                birth <- get_birth_fill(state, N0 = N) 
                }

                ## [Step 5]: Recovery of I.
                recover    <- rbinom_mat(n = length(c(state$Imat)), size = c(state$Imat), prob = c(state$gamma), nrow = nrow(state$Imat), ncol = ncol(state$Imat))

                ## [Step 6]: Mutation of new infections. Fraction of new infections -> mutation
               if (host_dyn_only == FALSE & sum(newinf) != 0) {
                  mutated  <- rbinom_mat(n = newinf, size = newinf, prob = mu, nrow = nrow(newinf), ncol = ncol(newinf))
               } else {
                  mutated  <- rbinom_mat(n = nrow(newinf), size = newinf, prob = mu, nrow = nrow(newinf), ncol = ncol(newinf))
               }

                ## [Step 7]: Mutation of new hosts (during birth).
                mutated_host <- rbinom(length(birth), size = birth, prob = mu/mut_host_mu_shift)
                ## Of the mutated hosts sort into resistance and tolerance mutants
                mutated_host_r <- rbinom(length(mutated_host), size = mutated_host, prob = mut_host_res_bias)
                mutated_host_t <- mutated_host - mutated_host_r

                ## debug to force mutation to check if that is working
                if (debug3 == TRUE) { 
                  mutated[1,1] <- 2
                  mutated_host[1] <- 2
                  mutated_host_r[1] <- 1
                  mutated_host_t[1] <- 1
                  }

                stopifnot(length(recover) == length(mutated))
                stopifnot(length(newinf)  == length(state$Imat))

                ## Updated birth from mutated hosts (these are the genetically identical birthed hosts)
                birth      <- birth - mutated_host

                ## [Step 8]: Update Infecteds
                state$Imat <- state$Imat - recover + newinf - mutated

                dfun("before mutation")
                if (debug) print(mutated)

                ## [Step 9]: Find new phenotypes for mutated parasites and hosts.
                tot_mut <- sum(mutated)

                ## Original parasite taits that mutants arise from. Need this for later
                orig_trait <- list(
                   pos_trait = rep(state$ltraitvec, colSums(mutated))
                 , neg_trait = rep(state$palphavec, colSums(mutated)))

                ## Mut parasite first. Another choice of order that may/may not matter
                if (tot_mut > 0) {
                    state <- do_mut(
                      state
                    , mut_var         = mut_var
                    , orig_trait      = orig_trait  ## ^^Just care about intrinsic nature of a strain
                    , mut_mean        = mut_mean
                    , mut_sd          = mut_sd
                    , mut_host        = FALSE
                    , mut_type        = mut_type
                    , power_c         = power_c
                    , power_exp       = power_exp
                    , mut_link_p      = mut_link_p
                    , agg_eff_adjust  = agg_eff_adjust
                    , parasite_tuning = parasite_tuning
                    , tradeoff_only   = tradeoff_only)
                }

                ## Host mutations
                if (length(which(is.na(mutated_host)) > 0)) {
                browser()
                }

                if (sum(mutated_host) > 0) {
                    state <- do_mut(
                      state
                    , mut_var             = mut_var
                    , orig_trait          = list(
                        res_mut           = rep(state$hrtraitvec, mutated_host_r)
                      , tol_mut           = rep(state$httraitvec, mutated_host_t))
                    , res_mut             = mutated_host_r
                    , tol_mut             = mutated_host_t
                    , mut_mean            = mut_mean
                    , mut_sd              = mut_sd
                    , mut_host            = TRUE
                    , mut_host_sd_shift   = mut_host_sd_shift
                    , mut_host_mean_shift = mut_host_mean_shift
                    , mut_type            = mut_type)
                }

                ## [Step 10]: Update Svec with infections and recoveries prior to the mutations (and host birth).
                state$Svec <- state$Svec + rowSums(recover) - rowSums(newinf) + birth

                ## [Step 11]: Update state (a big component beign parasite alpha and beta) with new parasite and host evolution.
                if (tot_mut > 0 | sum(mutated_host) > 0) {
                state <- update_mut_pt(
                    state           = state
                  , orig_trait      = orig_trait
                  , power_c         = power_c
                  , power_exp       = power_exp
                  , mut_link_p      = mut_link_p
                  , mut_link_h      = mut_link_h
                  , mutated         = mutated
                  , mutated_host    = mutated_host
                  , mut_var         = mut_var
                  , parasite_tuning = parasite_tuning
                  , eff_scale       = eff_scale)
                }

                if (sum(state$Imat)==0 & host_dyn_only == FALSE) {
                    message(sprintf("system went extinct prematurely (t=%d)", i))
                    break
                }
                if (sum(state$Svec)==0 & host_dyn_only == TRUE) {
                    message(sprintf("system went extinct prematurely (t=%d)", i))
                    break
                }

                ## Look over columns of I for extinct parasites and over rows of I and entries of S for extinct hosts
                extinct_p <- which(colSums(state$Imat) == 0)
                extinct_h <- which(rowSums(state$Imat) == 0 & state$Svec == 0)
                if (length(extinct_p) > 0 & host_dyn_only == FALSE) {
                    state <- do_extinct(state, mut_var ,extinct = extinct_p, parasite = TRUE)
                }
                if (length(extinct_h) > 0) {
                    state <- do_extinct(state, mut_var, extinct = extinct_h, parasite = FALSE)
                }
                dfun("after mutation")

            }  ## rptfreq time steps
      
        browser()
      
        ## summary statistics
        I_tot        <- ncol(state$Imat)
        num_I        <- sum(state$Imat)
        ltrait_mean  <- sum(colSums(state$Imat)*state$ltraitvec)/num_I
        ## Not quite sure about this formula. Taken from BMB code
        ltrait_sd    <- sqrt(sum(colSums(state$Imat)*(state$ltraitvec-ltrait_mean)^2)/num_I)
        lalpha_mean  <- sum(colSums(state$Imat)*state$palphavec)/num_I
        lalpha_sd    <- sqrt(sum(colSums(state$Imat)*(state$palphavec-lalpha_mean)^2)/num_I)
        num_S        <- sum(state$Svec)
        S_tot        <- length(state$Svec)
        lhres_mean   <- sum(state$Svec*state$hrtraitvec)/num_S
        lhres_sd     <- sqrt(sum(state$Svec*(state$hrtraitvec-lhres_mean)^2)/num_S)
        lhtol_mean   <- sum(state$Svec*state$httraitvec)/num_S
        lhtol_sd     <- sqrt(sum(state$Svec*(state$httraitvec-lhtol_mean)^2)/num_S)
        pop_size     <- num_I + num_S

        ## actual alpha and beta of all parasites in all host classes
        avg_alpha    <- mean(state$alpha)
        sd_alpha     <- mean(state$alpha)
        avg_beta     <- mean(state$beta)
        sd_beta      <- mean(state$beta)

        if (progress) cat(".")
        res[i,] <- c(
          i*rptfreq
        , num_S
        , S_tot
        , num_I
        , I_tot
        , pop_size
        , lhres_mean
        , lhres_sd
        , lhtol_mean
        , lhtol_sd
        , ltrait_mean
        , ltrait_sd
        , lalpha_mean
        , lalpha_sd
        , avg_alpha
        , sd_alpha
        , avg_beta
        , sd_beta
        , ifelse(mut_var == "beta", state$gamma[1,1], state$gamma[1,1]))

        ## DRY ...
        if (sum(state$Imat) == 0 & host_dyn_only == FALSE) {
            message(sprintf("system went extinct prematurely (t=%d)", i))
            break
        }
        if (sum(state$Svec) == 0 & host_dyn_only == TRUE) {
            message(sprintf("system went extinct prematurely (t=%d)", i))
            break
        }

    } ## loop over reporting frequencies
    if (progress) cat("\n")
    return(res)

    ## Run the deterministic model
    } else {

      ## convert state into a state list with only S and I and everything else put into params
      state_determ <- c(
        Svec = state$Svec
      , Imat = state$Imat
      )
      
      ## Create an alpha matrix instead of using an alpha vector
      alpha_mat <- mut_link_p$linkinv(matrix(data = rep(palphavec, each = length(tuning)), nrow = length(tuning), ncol = length(alpha0)))
      
      params_determ <- list(
        N      = N
      , d      = d
      , beta   = state$beta
      , alpha  = alpha_mat
      , gamma  = gamma0
      , mutlev = mu
      )

      ## Run deterministic model here
  par_evo.out <- as.data.frame(
    ode(
    y        = state_determ
  , times    = seq(0, determ_length, by = determ_timestep)
  , func     = par_evo_determ
  , parms    = params_determ
  , maxsteps = 5000
  , hini     = lsoda_hini
  , method   = 'rk4'
    ))

  par_evo.out

    }
}
