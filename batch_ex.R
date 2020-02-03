## 
## Run via Slurm job array: https://slurm.schedmd.com/job_array.html
## https://github.com/bbolker/betararef/tree/master/inst/batchfiles
## also see: https://cran.r-project.org/web/packages/future.batchtools/index.html

w <- as.integer(Sys.get("SLURM_ARRAY_TASK_ID"))
outfn <- sprintf("myjob_%d",w)
##
## figure out which rows of the data matrix to run
## run them
## send the output
