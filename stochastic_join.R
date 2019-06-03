####################################################
## Join the results of the stochastic simulations ## 
####################################################

stochastic_piece <- readRDS("res_out/res_out_DTS/res_1000_all_stochas_all.Rds")

stochastic_total  <- rbind(
  stochastic_piece
, stochastic_piece2
, stochastic_piece3
, stochastic_piece4
, stochastic_piece5
, stochastic_piece6
, stochastic_piece7)

stochastic_total <- rbind(
   stochastic_piece
,  stochastic_piece8
)

which_duplicated  <- which(duplicated(stochastic_total[, c(1, 27, 29, 30, 31, 43, 49, 52)]) == TRUE)
stochastic_total  <- stochastic_total[-which_duplicated, ]

saveRDS(stochastic_total, "res_1000_all_stochas_all.Rds")
