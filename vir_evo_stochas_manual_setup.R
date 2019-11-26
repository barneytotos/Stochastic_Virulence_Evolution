#############################################################################################
## Because I need to run such a sparse response surface design, set up parameters manually ##
#############################################################################################

######
## Set up parameter values to explore:
## 1) Time to equilibrium
## 2) Movement away from equilibrium / sd around equilibrium
## 3) Feedback of population size, mutation probability and sd on eco-evo feedbacks and 1 and 2 above
######

######
## Set up parameter values to range over:
## 1) 
## 2) 
## 3) 
## 4) 
##### 

######
## We will want to explore variation in these parameter values between:
## 1) Tradeoff curve with no second trait
## 2) Tradeoff curve with proportional beta
##### 

nt            <- 5e5
num_points    <- 1500
rptfreq       <- max(nt / num_points, 1) 
nrpt          <- nt %/% rptfreq
num_runs      <- 250
deterministic <- FALSE

