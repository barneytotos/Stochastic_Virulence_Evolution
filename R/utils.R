## from lme4; construct a list with names equal to deparse() by default
namedList <- function (...) 
{
    L <- list(...)
    snm <- sapply(substitute(list(...)), deparse)[-1]
    if (is.null(nm <- names(L))) 
        nm <- snm
    if (any(nonames <- nm == "")) 
        nm[nonames] <- snm[nonames]
    setNames(L, nm)
}

plot.evosim <- function(x, ..., type="postrait") {
    require("ggplot2")
    switch(type,
           postrait= (ggplot(x, aes(time, y=mean_postrait,
                                    ymin=mean_postrait-sd_postrait,
                                    ymax=mean_postrait+sd_postrait))
               + geom_line()
               + geom_ribbon(colour=NA,alpha=0.2)
           ) ## postrait
           ) ## switch
}
    

