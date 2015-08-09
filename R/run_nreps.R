#' Automatically determine sample size for an algorithm on a problem instance
#'
#' Iteratively calculates the required sample size for an algorithm on a
#'    problem, so that the final number of replicates yields an estimation of
#'    mean performance with a fixed uncertainty.
#'
#' The details of this routine come here...
#'
#' @section Problem structure:
#' The \code{problem} list must contain all relevant parameters that define the
#' problem instance. \code{problem} is a named list containing at least the
#' following fields:
#' \itemize{
#'    \item \code{$name} - name of the problem instance function, that is, a
#'          routine that calculates y = f(x)
#'    \item \code{$xmin} - vector of lower bounds of each variable
#'    \item \code{$xmax} - vector of upper bounds of each variable
#' }
#'
#' All other required parameters for the problem must be included as fields in
#' this list.
#'
#' @section Algorithm structure:
#' The \code{algo} list must contain all relevant parameters that define the
#' algorithm that will be applied to solve the \code{problem}. \code{algo} is a
#' named list containing the following fields:
#' \itemize{
#'    \item \code{$name} - name of the function that calls the algorithm
#'    \item \code{$pars} - named list of all parameters needed to set up the
#'          algorithm (e.g., population size, stop criteria, operator names and
#'          parameters, stop criteria, etc.).
#' }
#'
#' The function defined by the routine \code{algo$name} must have the following
#' structure:
#'
#'    \preformatted{
#'          myalgo <- function(pars, problem){
#'                ...
#'                return(results)
#'          }
#'    }
#'
#' That is, it must be able to run if called as
#'    \code{myalgo(algo$pars, problem)}.
#'
#' The algorithm routine should return a list containing at least the function
#' value obtained (result$f). More information can also be returned (such as the
#' solution found, number of function evaluations, etc.) but this additional
#' info is not used by \code{get_observation()}.
#'
#'@section Random seed:
#' The\code{seed} argument receives the desired seed for the PRNG (used, for
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
#' @param problem a list object containing the definitions of the problem
#'    instance. See \code{\link{Problem structure}} for details.
#' @param algo a list object containing the definitions of the algorithm See
#'    \code{\link{Algorithm structure}} for details.
#' @param dmax desired confidence interval halfwidth for the estimated mean
#'    performance of algorithm \code{algo} on instance \code{problem}.
#' @param method method used to calculate the confidence interval. Currently
#'      supported methods are "parametric" and "bootstrap"
#' @param alpha significance level for the confidence intervals.
#'      Defaults to \code{alpha = 0.05}.
#' @param nstart initial number of algorithm runs.
#'      Defaults to \code{nstart = 25}. See \code{\link{Warning}} for details.
#' @param nmax maximum allowed sample size.
#'      Defaults to \code{nmax = Inf}.
#' @param seed seed for the random number generator.
#'      Defaults to NULL. See \code{\link{Random seed}} for details.
#' @param nboot number of bootstrap samples (if \code{method = "bootstrap"}).
#'      Defaults to \code{nboot = 1000}.
#' @param ncpus number of cores to use for bootstrap
#'      (if \code{method} = "bootstrap"). Defaults to \code{ncpus = 1}.
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
    while((delta > dmax) & (n < nmax)){
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
