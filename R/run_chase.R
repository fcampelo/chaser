#' Comparison of Heuristics with Adaptive Sample size Estimation
#'
#' Run experiment using adaptive sample size estimation.
#'
#' The details of this routine come here...
#'
#' @section Problems and algorithms:
#' Parameters \code{problems} and \code{algorithms} must be a list of
#' \code{problem} (\code{algo}) specifications, each one a list itself.
#' See \code{\link{run_nreps}} for details on the specification of each
#' \code{problem} and \code{algo} component.
#'
#'@section Random seed:
#' The \code{seed} argument receives the desired seed for the PRNG (used, for
#' instance, in the bootstrap resampling). This value can be set for
#' reproducibility purposes. The value of this parameter defaults to NULL, in
#' which case the seed is arbitrarily set using \code{as.numeric(Sys.time())}.
#' The value used as seed is always reported in the output structure.
#'
#' @section Warning:
#' For distribution-free techniques the initial number of observations
#' \code{nstart} should be relatively high (> 25 if extreme outliers are not
#' expected, > 50 for the general case) to guarantee good statistical
#' properties (particularly compliance with the nominal type-I error
#' rate \code{alpha}). However, if some distributional assumptions can be
#' made (e.g., low skewness of the population), then \code{nstart} can in
#' principle be as small as 5.
#' In general, higher sample sizes are the price to pay for abandoning
#' distributional assumptions. Use lower values of \code{nstart} with caution.
#'
#'@seealso \code{\link{run_nreps}} for more information on the definition of
#'  problem and algorithm lists
#'
#' @param problems a list object containing lists defining all problem
#'    instances to be used in the experiment.
#'    See \code{\link{Problems and algorithms}} for details.
#' @param algorithms a list object containing lists defining all algorithms to
#'    be used in the experiment.
#'    See \code{\link{Problems and algorithms}} for details.
#' @param dmax desired confidence interval halfwidth for the estimated mean
#'    performance of each algorithm on each instance.
#' @param method method used to calculate the sample size for each pair
#'    algorithm-problem.
#' @param alpha significance level for the confidence intervals on the means of
#'    each algo-problem pair. Defaults to \code{alpha = 0.05}.
#' @param nstart initial number of algorithm runs.
#'      Defaults to \code{nstart = 25}.
#'      See \link{Warning} for details.
#' @param nmax maximum allowed sample size. Defaults to \code{nmax = Inf}.
#' @param seed seed for the random number generator.
#'      Defaults to NULL.
#'      See \link{Random seed} for details.
#' @param nboot number of bootstrap samples (if \code{method = "bootstrap"}).
#'      Defaults to \code{nboot = 1000}.
#' @param ncpus number of cores to use for bootstrap
#'      (if \code{method} = "bootstrap"). Defaults to \code{ncpus = 1}.
#'
#' @return a list object containing the following items:
#' \itemize{
#'    \item \code{rawdata} - data frame containing all observations generated
#'    \item \code{datameans} - data frame containing the means of each algorithm
#'    on each problem.
#'    \item \code{datasterr} - data frame containing the standard errors of the
#'    means for each algorithm on each problem.
#' }
#'
#' @author Felipe Campelo
#' @export

run_chase <- function(problems,
                      algorithms,
                      dmax,
                      method = c("parametric", "bootstrap"),
                      alpha = 0.05,
                      nstart = 25,
                      nmax = Inf,
                      seed = NULL,
                      nboot = 1000,
                      ncpus = 1)
{
    # Recover input information
    nprobs <- length(problems)
    nalgos <- length(algorithms)

    if(is.null(seed)) {seed <- as.numeric(Sys.time())}

    # Initialize output data frames
    rawdata <- data.frame(Algorithm = character(),
                          Instance  = character(),
                          Value     = numeric())
    datameans <- rawdata
    datasterr <- rawdata

    # Iterative cycle
    for (i in 1:nalgos){
        for (j in 1:nprobs){
            data.ij <- run_nreps(problem = problems[[j]],
                                 algo    = algorithms[[i]],
                                 dmax    = dmax,
                                 method  = method,
                                 alpha   = alpha,
                                 nstart  = nstart,
                                 nmax    = nmax,
                                 seed    = seed,
                                 nboot   = nboot,
                                 ncpus   = ncpus)
            rawdata     <- rbind(rawdata,
                                 data.frame(Algorithm = rep(algorithms[[i]]$name,
                                                            times = data.ij$n),
                                            Instance  = rep(problems[[j]]$name,
                                                            times = data.ij$n),
                                            Value     = data.ij$xbar))
            datameans   <- rbind(datameans,
                                 data.frame(Algorithm = algorithms[[i]]$name,
                                            Instance  = problems[[j]]$name,
                                            Value     = data.ij$xbar))
            datasterr   <- rbind(datasterr,
                                 data.frame(Algorithm = algorithms[[i]]$name,
                                            Instance  = problems[[j]]$name,
                                            Value     = data.ij$se))
        }
    }

    output<-list(xbar = mean(x),
                 s = sd(x),
                 n = n,
                 delta = delta,
                 ci = ci,
                 se = se,
                 x = x,
                 seed = seed,
                 alpha = alpha)

    return(output)
}
