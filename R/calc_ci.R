#' Calculates confidence interval for a specified statistic
#'
#' Calculates the confidence interval for a statistic, based on a given sample.
#'
#' The details of this routine come here...
#' (s.e. for median is based on asymptotic formula for normal
#' distribution/large sample size, see
#' R.R. Sokal, F.J. Rohlf, "Biometry", W.H. Freeman, 4th ed., 2011, p. 139)
#'
#' @param x vector of observations
#' @param FUN name of statistic function
#'      Defaults to \code{FUN = "mean"}.
#' @param alpha desired significance level for the interval.
#'      Defaults to \code{alpha = 0.05}.
#' @param method method used to calculate the interval.
#'      Defaults to \code{method = "parametric"}.
#' @param nboot number of bootstrap resamples (if \code{method} = "bootstrap").
#'      Defaults to \code{nboot = 1000}.
#' @param ncpus number of cores to use for bootstrap
#'      (if \code{method} = "bootstrap"). Defaults to \code{ncpus = 4}.
#'
#' @return a list object containing the following items:
#' \itemize{
#'    \item \code{ci} - numeric vector of length 2 containing the confidence
#'      interval
#'    \item \code{se} - standard error of the base statistic
#' }
#'
#' @author Felipe Campelo (\email{fcampelo@@ufmg.br}), Fernanda Takahashi (\email{fernandact@@ufmg.br})

calc_ci <- function(x,                    # numeric:   vector of observations
                    FUN = "mean",         # character: statistic function
                    alpha = 0.05,         # numeric:   significance level
                    method = "parametric",# character: method for calculating CI
                    nboot = 1000,         # integer:   bootstrap samples
                    ncpus = 4             # integer:   number of cores
){

  # ========== Error catching ========== #
  assertthat::assert_that(
    is.numeric(x) && is.vector(x) && length(x) > 1,
    is.character(FUN) && length(FUN) == 1,
    FUN %in% c("mean", "median"),
    is.numeric(alpha) && alpha > 0 && alpha < 1,
    is.character(method) && length(method) == 1,
    method %in% c("parametric", "bootstrap"),
    assertthat::is.count(nboot) && assertthat::is.count(ncpus))
  # ==================================== #
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
                      ci <-rep(boot.x$t0, 2))
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
