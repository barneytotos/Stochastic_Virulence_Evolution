#######################################
## Load and clean up stochastic runs ## 
#######################################

## Some place you previously saved the file
# stochas.res <- readRDS("res_out/res_out_DTS/res_1000_all_stochas_all.Rds")

## Maybe not a complete run, so droplevels first
stochas.res <- droplevels(stochas.res)

## See what is covered in this set of runs
with(stochas.res, list(unique(mu), unique(mut_sd), unique(alpha0), unique(tune0), unique(N)))

## For this first pass using res_1000_all_stochas, examine some results from only the starting 
 ## location of c(0.03, 0.97)
stochas.res_subset <- stochas.res %>% 
  filter(
    mu        == 0.005 &
    mut_sd    == 0.05 & 
    N         == 200  & 
    eff_scale == 30   &
    alpha0    == 0.03 &
    tune0     == 0.97
    )

## Check that the subset data frame has only the parameters that I want
with(stochas.res_subset, list(mu = unique(mu), mut_sd = unique(mut_sd), alpha0 = unique(alpha0), tune0 = unique(tune0), N = unique(N)))
paste("Num runs = ", nrow(stochas.res_subset) / length(unique(stochas.res_subset$time)), sep = " ")

stochas.res_mu0.001 <- stochas.res %>% filter(mu == 0.001 & N == 600)
stochas.res_mu0.005 <- stochas.res %>% filter(mu == 0.005 & N == 600)
stochas.res_mu0.025 <- stochas.res %>% filter(mu == 0.025 & N == 600)

stochas.res_mutsd0.05 <- stochas.res %>% filter(mut_sd == 0.05 & N == 600 & eff_scale == 30)
stochas.res_mutsd0.01 <- stochas.res %>% filter(mut_sd == 0.01 & N == 600)
stochas.res_mutsd0.25 <- stochas.res %>% filter(mut_sd == 0.25 & N == 600)

stochas.res_Ncheck  <- stochas.res %>% filter(mu == 0.005 & mut_sd == 0.05)

stochas.res_N200  <- stochas.res %>% filter(N == 200 & mu == 0.005 & mut_sd == 0.05)
stochas.res_N600  <- stochas.res %>% filter(N == 600 & mu == 0.005 & mut_sd == 0.05)
stochas.res_N1800 <- stochas.res %>% filter(N == 1800 & mu == 0.005 & mut_sd == 0.05)

## remove the full results data frame to save space
rm(stochas.res)
