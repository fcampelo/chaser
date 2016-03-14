#' Calculates confidence interval for a specified statistic
#'
#' Calculates the confidence interval for a population parameter, based on a
#' given sample.
#'
#' This routine is used to calculate confidence intervals on the mean or median,
#' based either in the parametric formulas or in bootstrap resampling.
#'
#' @section References:
#' R.R. Sokal, F.J. Rohlf.
#' "Biometry".
#' W.H. Freeman, 4th ed., 2011, p. 139
#'
#' D.C. Montgomery, G.C. Runger.
#' "Applied Statistics and Probability for Engineers"
#' Wiley, 6th ed., 2013
#'
#' A.C. Davison, D.V. Hinkley.
#' "Bootstrap Methods and Their Application".
#' Cambridge University Press, 1st ed., 1997
#'
#' @param x vector of observations
#' @param stat name of statistic for which the CI is desired.
#' @param alpha desired significance level for the interval.
#' @param method method used to calculate the interval. Accepts "param" (the
#'          default), "boot" (for bootstrap) or "binom" (for binomial median CI)
#' @param ... further parameters to be passed on to \code{boot}
#'          (if method == "boot")
#'
#' @return a list object containing the following items:
#' \itemize{
#'    \item \code{x.est} - estimated value
#'    \item \code{ci} - numeric vector of length 2 containing the confidence
#'      interval
#'    \item \code{se} - standard error
#' }
#'
#' @importFrom boot boot
#' @importFrom boot boot.ci
#'
#' @author Felipe Campelo (\email{fcampelo@@ufmg.br}),
#'         Fernanda Takahashi (\email{fernandact@@ufmg.br})
#'
#' @examples
#' \dontrun{
#' x <- rnorm(n = 30, mean = 5, sd = 3)
#' calc_ci(x) # parametric 95% CI on the mean
#' calc_ci(x, stat = "median", alpha = .01) # parametric 99% CI on the median
#'
#' x <- rexp(n = 25, rate = 0.5)
#' # bootstrap 95% CI on the mean
#' calc_ci(x, method = "boot",
#'         R = 1000, parallel = "multicore", ncpus = 4)}


calc_ci <- function(x,                  # numeric:   vector of observations
                    stat   = "mean",    # character: statistic function
                    method = "param",   # character: method for calculating CI
                    alpha  = 0.05,      # numeric:   significance level
                    ...)                # other arguments for boot (ncpus, etc.)
{

    # ========== Error catching ========== #
    assertthat::assert_that(
        is.numeric(x), is.vector(x), length(x) > 1,
        is.character(stat), length(stat) == 1,
        stat %in% c("mean", "median"),
        is.numeric(alpha), alpha > 0, alpha < 1,
        is.character(method), length(method) == 1,
        method %in% c("param", "boot") || (method == "binom" && stat == "median"))
    # ==================================== #
    n <- length(x)

    # Prepare bootstrap parameters
    if (method == "boot"){
        # define bootstrap parameters
        boot.R     <- 1000
        boot.ncpus <- 1
        boot.par   <- "multicore"
        myarg <- list(...)
        if ("R" %in% names(myarg))        boot.R     <- myarg$R
        if ("parallel" %in% names(myarg)) boot.par   <- myarg$parallel
        if ("ncpus" %in% names(myarg))    boot.ncpus <- myarg$ncpus

        # Define the bootstrap function
        bootfun <- switch(stat,
                          mean   = ,
                          median = function(x, d){return(median(x[d]))})
    }

    # ==================================== #

    # CI on the mean
    if (stat == "mean"){
        if (method == "param"){
            est    <- mean(x)
            se     <- sd(x) / sqrt(n)
            ci     <- est + c(-1, 1) * qt(1 - alpha/2, n - 1) * se
        } else if (method == "boot"){
            boot.x <- boot::boot(data      = x,
                                 statistic = function(x, d){return(mean(x[d]))},
                                 R         = boot.R,
                                 parallel  = boot.par,
                                 ncpus     = boot.ncpus)
            est <- boot.x$t0
            se  <- sd(boot.x$t)
            tmp <- ifelse(length(unique(x)) > 1,
                          ci <- boot::boot.ci(boot.out = boot.x,
                                              conf     = 1 - alpha,
                                              type     = "bca")$bca[4:5],
                          ci <- rep(boot.x$t0, 2))
        } else {stop("method = ", method, " not available for for stat = 'mean'")}
    } else if (stat == "median"){
        if (method == "param"){
            lCIndx <- max(1, floor((n + qnorm(alpha/2) * sqrt(n)) / 2))
            uCIndx <- min(n, ceiling(1 + (n + qnorm(1-alpha/2) * sqrt(n)) / 2))
            est    <- median(x)
            ci     <- sort(x)[c(lCIndx, uCIndx)]
            se     <- 1.253 * sd(x) / sqrt(n)
        } else if (method == "boot"){
            boot.x <- boot::boot(data      = x,
                                 statistic = function(x, d){return(median(x[d]))},
                                 R         = boot.R,
                                 parallel  = boot.par,
                                 ncpus     = boot.ncpus)
            est <- boot.x$t0
            se  <- sd(boot.x$t)
            tmp <- ifelse(length(unique(x)) > 1,
                          ci <- boot::boot.ci(boot.out = boot.x,
                                              conf     = 1 - alpha,
                                              type     = "bca")$bca[4:5],
                          ci <- rep(boot.x$t0, 2))
        } else if (method == "binom"){
            ci  <- sort(x)[qbinom(c(alpha / 2,
                                    1 - alpha / 2),
                                  length(x), 0.5)]
            est <- median(x)
            se  <- NA
        } else {stop("unrecognized method for calc_ci()")}
    } else {stop("method = ", method, " not available for for stat = 'median'")}

    return(list(est = est,
                ci  = ci,
                se  = se))
}
