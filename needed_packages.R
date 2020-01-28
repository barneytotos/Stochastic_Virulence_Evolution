#########################
### Required packages ###
#########################

needed_packages <- c(
  "arm"
, "ggplot2"
, "gridExtra"
, "dplyr"
, "ReacTran"
, "deSolve"
, "abind"
, "gganimate"
, "reshape2"
, "lhs"
, "vegan"
<<<<<<< HEAD
, "viridis"
=======
>>>>>>> 48da966643a47dda183e2eb0002c0fe0d670645e
  )

## Check if the packages are installed. If they are not install them, then load them
if (length(setdiff(needed_packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(needed_packages, rownames(installed.packages())))  
}

lapply(needed_packages, require, character.only = TRUE)
