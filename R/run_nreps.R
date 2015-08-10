#' Automatically determine sample size for an algorithm on a problem instance
#'
#' Iteratively calculates the required sample size for an algorithm on a
#'    problem, so that the final number of replicates yields an estimation of
#'    mean performance with a fixed uncertainty.
#'
#' The details of this routine come here...
#'
#'
#' @param problem a list object containing the definitions of the problem
#'    instance. See \code{\link{run_chase}} for details.
#' @param algo a list object containing the definitions of the algorithm See
#'    \code{\link{run_chase}} for details.
#' @param dmax desired confidence interval halfwidth for the estimated mean
#'    performance of algorithm \code{algo} on instance \code{problem}.
#' @param method method used to calculate the confidence interval. Currently
#'      supported methods are "parametric" and "bootstrap"
#' @param alpha significance level for the confidence intervals.
#' @param nstart initial number of algorithm runs.
#'      See \code{\link{run_chase}} for details.
#' @param nmax maximum allowed sample size.
#'      Defaults to \code{nmax = Inf}.
#' @param seed seed for the random number generator.
#'      Defaults to NULL. See \code{\link{run_chase}} for details.
#' @param nboot number of bootstrap samples.
#' @param ncpus number of cores to use for bootstrap.
#'
#' @return a list object containing the following items:
#' \itemize{
#'    \item \code{observations} - vector of observed performance values
#'    \item \code{xbar} - sample mean of the performance of \code{algo} on \code{problem}
#'    \item \code{se} - standard error of the mean
#'    \item \code{n} - number of observations
#'    \item \code{seed} - the seed used for the PRNG
#' }
#'
#' @seealso run_chase
#' @author Felipe Campelo
#' @export

run_nreps <- function(problem,        # problem parameters
                      algo,           # algorithm parameters
                      dmax,           # maximum allowed CI halfwidth
                      method = c("parametric", "bootstrap"), # method for obtaining the CI
                      alpha = 0.05,   # significance level for CI
                      nstart = 25,    # initial number of samples
                      nmax = Inf,     # maximum allowed sample size
                      seed = NULL,    # seed for PRNG
                      nboot = 1000,   # number of bootstrap samples
                      ncpus = 1       # number of cores to use for bootstrap
){

    # initial definitions
    n           <- nstart - 1       # initial number of observations
    delta       <- Inf              # initial CI halfwidth

    # set PRNG seed
    if(is.null(seed)) {seed <- as.numeric(Sys.time())}
    set.seed(seed)

    # initialize vector of observations
    x <- get_observations(algo, problem, n)

    # Iterative cycle
    while(!(delta < dmax | n > nmax)){
        # Generate a new observation
        x       <- c(x,
                     get_observations(algo, problem))
        n       <- n + 1

        # Calculate CI halfwidth
        CI      <- calc_ci(x,
                           "mean",
                           alpha,
                           method,
                           nboot)
        delta   <- diff(CI$ci)/2
    }

    output<-list(observations = x,
                 xbar         = mean(x),
                 n            = n,
                 se           = CI$se,
                 seed         = seed)

    return(output)
}
