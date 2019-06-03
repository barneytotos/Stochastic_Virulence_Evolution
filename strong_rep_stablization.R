pop_slots <- c(rep(0, 220), rep(1, 380))
N0 <- 33

birth_1000 <- data.frame(
  time    = seq(1, 1001)
, N_start = numeric(1001)
, N_end   = numeric(1001)
)

birth_1000[1, 2] <- N0

for (i in 1:1000) {
  
  ## death (background)
birth_1000[i, 3] <- birth_1000[i, 2] - rbinom(1, birth_1000[i, 2], 0.01)
  ## death (due to virulence)
birth_1000[i, 3] <- birth_1000[i, 3] - rbinom(1, birth_1000[i, 3], 0.21)
  ## birth
if (birth_1000[i, 3] > N0) {
  birth_1000[(i+1), 2] <- birth_1000[i, 3]
} else {
  birth_1000[(i+1), 2] <- birth_1000[i, 3] + 
    rbinom(
      n = 1
    , size = N0
    , prob = (N0 - birth_1000[i, 3]) / N0)
}

}

ggplot(birth_1000, aes(time, N_start)) + geom_line()

 
hist(
rbinom(n = 10000, size = sum(pop_slots), prob = (N0 - sum(pop_slots)) / N0)
  , breaks = 200
)

## calculate probability of birth

birth_prob <- data.frame(
  N0   = rep(600, 600)
, N    = seq(1, 600, by = 1)
, prob = numeric(600)
)

birth_prob <- transform(birth_prob
  , prob = 1 - (N / N0)
  )

ggplot(birth_prob, aes(N, prob)) + geom_line()

