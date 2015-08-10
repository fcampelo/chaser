#' Comparison of Heuristics with Adaptive Sample size Estimation
#'
#' Run experiment using adaptive sample size estimation.
#'
#' The details of this routine come here...
#'
#' @section Problems and Algorithms:
#' Parameters \code{problems} and \code{algorithms} must each be a list of
#' \code{problem} (\code{algorithm}) specifications, defined according to the
#' instructions given below.
#'
#' Each component of \code{problems} is a list containing all relevant
#' parameters that define the problem instance. \code{problems[[j]]} must be
#' a named list containing at least the following fields:
#'
#' \itemize{
#'    \item \code{$name} - name of the problem instance function, that is, a
#'          routine that calculates y = f(x)
#'    \item \code{$xmin} - vector of lower bounds of each variable
#'    \item \code{$xmax} - vector of upper bounds of each variable
#' }
#'
#' If the problem requires additional parameters, these must also be provided
#' as fields.
#'
#' Similarly, the \code{algorithms} parameter must contain a list where each
#' element provides the full specification of an algorithm. More specifically,
#' \code{algorithms[[i]]} is a named list containing all relevant
#' parameters that define an algorithm that will be applied to solve the
#' \code{problems}.
#'
#' \code{algorithms[[i]]} must contain a \code{$name} field (the name of the
#' function that calls the algorithm) and any other elements/parameters that
#' \code{algorithms[[i]]$name} requires (e.g., population size, stop criteria,
#' operator names and parameters, stop criteria etc.).
#'
#' The function defined by the routine \code{algorithms[[i]]$name} must have the
#' following structure: supposing that the list in \code{algorithms[[1]]} has
#' fields \code{$name = myalgo} and \code{$par1, $par2}, then:
#'
#'    \preformatted{
#'          myalgo <- function(par1, par2, probpars, ...){
#'                # do stuff
#'                # ...
#'                return(results)
#'          }
#'    }
#'
#' That is, it must be able to run if called as:
#'
#'    \preformatted{
#'          # remove '$name' field from list of arguments
#'          # and include the problem definition as field 'probpars'
#'          myargs          <- algorithms[[1]][names(algorithms[[1]]) != "name"]
#'          myargs$probpars <- probpars
#'
#'          # call function
#'          do.call(algorithms[[1]]$name,
#'                  args = myargs)
#'    }
#'
#' The \code{algorithms[[i]]$name} routine must return a list containing (at
#' least) the function value of the final solution obtained
#' (\code{result$Fbest)} after a given run.
#'
#' @section Random Seed:
#' The \code{seed} argument receives the desired seed for the PRNG (used, for
#' instance, in the bootstrap resampling). This value can be set for
#' reproducibility purposes. The value of this parameter defaults to NULL, in
#' which case the seed is arbitrarily set using \code{as.numeric(Sys.time())}.
#' The value used as seed is always reported in the output structure.
#'
#' Notice that this \code{seed} is not passed to the algorithm routines in the
#' current version.
#'
#' @section Warning:
#' In the general case the initial number of observations
#' (\code{nstart}) should be relatively high (> 25 if extreme outliers are not
#' expected, > 50 if that assumption can't be made) to guarantee good
#' statistical properties (particularly compliance with nominal type-I error
#' rate \code{alpha}). However, if some distributional assumptions can be
#' made (e.g., low skewness of the population), then \code{nstart} can in
#' principle be as small as 5.
#'
#' In general, higher sample sizes are the price to pay for abandoning
#' distributional assumptions. Use lower values of \code{nstart} with caution.
#'
#'
#' @param problems a list object containing lists defining all problem
#'    instances to be used in the experiment.
#'    See \code{Problems and Algorithms} for details.
#' @param algorithms a list object containing lists defining all algorithms to
#'    be used in the experiment.
#'    See \code{Problems and Algorithms} for details.
#' @param dmax desired confidence interval halfwidth for the estimated mean
#'    performance of each algorithm on each instance.
#' @param method method used to calculate the sample size for each pair
#'    algorithm-problem.
#' @param stat statistic of interest (mean, median, ...)
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
#'
#' @author Felipe Campelo (\email{fcampelo@@ufmg.br})
#' @export

run_chase <- function(problems,
                      algorithms,
                      dmax,
                      method = c("parametric", "bootstrap"),
                      stat   = c("mean", "median"),
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
