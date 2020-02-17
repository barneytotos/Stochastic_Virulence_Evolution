## library(hpevosim)
library(testthat)

## FIXME: add context() statements
## FIXME: shorten runs to speed up tests?

sum_vars <- c("num_S","num_I","num_I_strains","mean_negtrait","sd_negtrait",
              "mean_postrait","sd_postrait","shannon")

run_sim_params <- stoch_params0[intersect(names(stoch_params0),names(formals(run_sim)))]
res <- do.call(run_sim,run_sim_params)

expect_equal(nrow(res),50)
expect_equal(ncol(res),18)
expect_true(all(res$pop_size==400))

m <- colMeans(tail(res[sum_vars],10))
## dput(round(m,8))
expect_equal(m,
             c(num_S = 89.8, num_I = 310.2, num_I_strains = 15.9,
               mean_negtrait = 0.1216931, sd_negtrait = 0.01837436,
               mean_postrait = 0.0493278, sd_postrait = 0.0025132, 
               shannon = 1.79367083)
)

plot(res)  ## goes to look for plot.evosim function

## re-running without tradeoff
run_sim_params_nt <- transform(run_sim_params,
                                    no_tradeoff = TRUE,
                                    mut_mean = -0.05)
res_nt <- do.call(run_sim,run_sim_params_nt)

plot(res_nt)

m_nt <- colMeans(tail(res_nt[sum_vars],10))

## dput(round(m_nt,8))
expect_equal(m_nt,
             c(num_S = 69.3, num_I = 330.7, num_I_strains = 17.2,
               ## fixed at orig values
               mean_negtrait = 0.01, sd_negtrait = 0,
               ## arbitrarily close to 1
               mean_postrait = 0.98985252, sd_postrait = 0.00909884, 
               shannon = 1.99599423)
             )

## running again, this time with a deterministic model.
## evolve only in beta
## imatseed
run_sim_params_determ <- transform(run_sim_params_nt,
                                    deterministic = TRUE,
                                    mut_mean = -0.02,
                                    mu = 0.1,
                                    gamma0 = 0.02)

res_det <- do.call(run_sim,run_sim_params_determ)
## FIXME: progress bar?

## FIXME: Need to improve the way these results are returned and the clarity of these checks.
m_det <- colSums(res_det)
nt <- nrow(res_det)
nI_cat <- sqrt(ncol(res_det)-2)
names(m_det) <- NULL
expect_equal(which(m_det != 0), c(1,2,seq(5403, 5502, by = 1)))

m_det <- m_det[which(m_det != 0)]

expect_equal(m_det[1:5],
        c(4100.0000, 9.8869994595, 0.5694942941, 0.6238255844, 0.7240922323))


image(aa[,,55])
res_det_fp <- matrix(data = unlist(res_det[nt, -c(1, 2)]),
                     nrow = nI_cat, ncol = nI_cat, byrow = TRUE)

## row 55 is initial infected genotype
init_gamma_cat <- eval(formals(run_sim)$Imat_seed)[2]
beta_prob <- rev(seq(0.00, 0.99, by = 0.01))
plot(beta_prob,
     res_det_fp[init_gamma_cat, ],
     xlab = "Beta", ylab = "Proportion") 

if (require(Matrix)) {
    image(Matrix::Matrix(res_det[,-(1:2)]
}
## base-R plotting
aa <- array(unlist(res_det[,-(1:2)]),c(41,100,100))
dimnames(aa) <- list(time=0:40, ## not really
                     beta=beta_prob,
                     gamma=beta_prob) ## the same in this case
bmat <- aa[,,init_gamma_cat]
bmat_flip <- bmat[,rev(seq(ncol(bmat)))]

par(las=1)
image(y=as.numeric(colnames(bmat_flip)),
      z=bmat_flip, xlab="time",ylab="beta prob")

## normalize; if you normalize by column (i.e. distribution across time)
## it looks weird and interesting but maybe meaningless)
bmat_flip_norm <- sweep(bmat_flip,MARGIN=1,rowSums(bmat_flip),"/")
image(y=as.numeric(colnames(bmat_flip)),
      z=bmat_flip_norm, xlab="time",ylab="beta prob",
      zlim=c(0,0.4) ## restrict scale so not overwhelmed by initial
                    ## concentrations (could use log10(0.01+x))
      )

if (FALSE) {
    library(rgl)
    persp3d(y=as.numeric(colnames(bmat_flip)),
            z=bmat_flip_norm, xlab="time",ylab="beta prob",
            zlim=c(0,0.4), ## restrict scale so not overwhelmed by initial
            ## concentrations (could use log10(0.01+x))
            col="gray"
            )
}

## could do this in tidyverse instead ...

