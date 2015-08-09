#' Calculates confidence interval for a specified statistic
#'
#' Calculates the confidence interval for a statistic based on a given sample.
#'
#' The details of this routine come here... (s.e. for median is based on asymptotic formula for normal distro/large smaple size, ref http://www.amazon.com/dp/0716786044/?tag=stackoverfl08-20 pg 139)
#'
#' @param x vector of observations
#' @param FUN base statistic for the confidence interval.
#'      Defaults to \code{FUN = "mean"}.
#' @param alpha desired significance level for the interval.
#'      Defaults to \code{alpha = 0.05}.
#' @param method method used to calculate the interval.
#'      Defaults to \code{method = "parametric"}.
#' @param nboot number of bootstrap resamples (if \code{method} = "bootstrap").
#'      Defaults to \code{nboot = 1000}.
#' @param ncpus number of cores to use for bootstrap
#'      (if \code{method} = "bootstrap"). Defaults to \code{ncpus = 1}.
#'
#' @return a list object containing the following items:
#' \itemize{
#'    \item \code{ci} - numeric vector of length 2 containing the confidence
#'      interval
#'    \item \code{se} - standard error of the base statistic
#' }
#'
#' @author Felipe Campelo
#' @export

calc_ci <- function(x,                  # numeric:   vector of observations
                    FUN = "mean",       # character: statistical function to use
                    alpha = 0.05,       # numeric:   significance level for CI
                    method = "parametric", # character: method for obtaining the CI
                    nboot = 1000,       # integer:   number of bootstrap samples
                    ncpus = 1           # integer:   number of cores to use for
                                        #            parallel processing
){
    n <- length(x)

    # Using bootstrap:
    if (method == "bootstrap"){
        # Define the bootstrap function
        bootfun <- switch(FUN,
                          mean = function(x,d){return(mean(x[d]))},
                          median = function(x,d){return(median(x[d]))})

        # Perform bootstrap
        boot.x <- boot::boot(data = x,
                             statistic = bootfun,
                             R = nboot,
                             parallel = "multicore",
                             ncpus = ncpus)

        # get CI
        tmp <- ifelse(length(unique(x)) > 1,
                      ci <- boot::boot.ci(boot.x,
                                          conf = 1 - alpha,
                                          type = "bca")$bca[4:5],
                      ci <-rep(botox$t0, 2))
        se  <- sd(boot.x$t)

    } else if (method == "parametric"){
        if (FUN == "mean"){
            se <- sd(x) / sqrt(n)
            ci <- mean(x) + c(-1, 1) * qt(1 - alpha/2, n - 1) * se
        } else if (FUN == "median"){
            lowCIndx    <- max(1, floor((n - 1.96*sqrt(n))/2))
            uppCIndx    <- min(n, ceiling(1 + (n + 1.96*sqrt(n))/2))
            rankedX     <- sort(x)
            ci          <- c(rankedX[lowCIndx], rankedX[uppCIndx])
            se          <- 1.253*sd(x)/sqrt(n)
        } else stop("Unrecognized statistical function used in calc_ci()")
    }
    CI <- list(ci = ci,
               se = se)
    return(CI)
}
