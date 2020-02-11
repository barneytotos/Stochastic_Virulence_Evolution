## from lme4; construct a list with names equal to deparse() by default
#' @importFrom stats setNames
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

##' @export
##' @importFrom ggplot2 ggplot aes geom_line geom_ribbon
plot.hpevosim <- function(x, ..., type="postrait") {
    time <- mean_postrait <- sd_postrait <- NULL ## global variable check
    switch(type,
           postrait= (ggplot(x, aes(time, y=mean_postrait,
                                    ymin=mean_postrait-sd_postrait,
                                    ymax=mean_postrait+sd_postrait))
               + geom_line()
               + geom_ribbon(colour=NA,alpha=0.2)
           ) ## postrait
           ) ## switch
}
    

