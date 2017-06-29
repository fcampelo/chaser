#' Automatically determine sample size for an algorithm on a problem instance
#'
#' Iteratively calculates the required sample size for an algorithm on a
#'    problem instance, so that the final number of replicates yields an
#'    estimation of expected performance with a predefined maximum uncertainty.
#'
#' @section Instances and Algorithms:
#' Parameters \code{instance} and \code{algorithm} must each be a list of
#' instance (algorithm) specifications, defined according to the instructions
#' given below.
#'
#' \code{instance} is a named list containing all relevant parameters that
#' define the problem instance. This list must contain at least the following
#' fields:
#'
#' \itemize{
#'    \item \code{$name} - name of the problem instance function, that is, a
#'          routine that calculates y = f(x)
#'    \item \code{$xmin} - vector of lower bounds of each variable
#'    \item \code{$xmax} - vector of upper bounds of each variable
#' }
#'
#' If the instance requires additional parameters, these must also be provided
#' as named fields.
#'
#' Similarly, \code{algorithm} must be a named list containing all relevant
#' parameters that define the algorithm to be applied for solving the problem
#' instance.
#'
#' \code{algorithm} must contain a \code{$name} field (the name of the
#' function that calls the algorithm) and any other elements/parameters that
#' \code{algorithm$name} requires (e.g., population size, stop criteria,
#' operator names and parameters, etc.).
#'
#' The function defined by the routine \code{algorithm$name} must have the
#' following structure: supposing that the list in \code{algorithm} has
#' fields \code{$name = myalgo} and \code{$par1 = "a", $par2 = 5}, then:
#'
#'    \preformatted{
#'          myalgo <- function(par1, par2, instance, ...){
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
#'          # and include the problem definition as field 'instance'
#'          myargs          <- algorithm[names(algorithm) != "name"]
#'          myargs$instance <- instance
#'
#'          # call function
#'          do.call(algorithm$name,
#'                  args = myargs)
#'    }
#'
#' The \code{algorithm$name} routine must return a list containing (at
#' least) the function value of the final solution obtained
#' (\code{result$Fbest)} after a given run.
#'
#' @section Initial Number of Observations:
#' In the general case the initial number of observations / algorithm / instance
#' (\code{nstart}) should be relatively high (> 20 if outliers are not
#' expected, > 50 (at least) if that assumption can't be made) to guarantee good
#' statistical properties (particularly compliance with nominal type-I error
#' rate \code{alpha}). However, if some distributional assumptions can be
#' made - particularly low skewness of the population of algorithm results on
#' the test instances), then \code{nstart} can in principle be as small as 5.
#'
#' In general, higher sample sizes are the price to pay for abandoning
#' distributional assumptions. Use lower values of \code{nstart} with caution.
#'
#' @inheritParams calc_ci
#' @param instance a list object containing the definitions of the problem
#'    instance. See Section \code{Problems and Algorithms} for details.
#' @param algo a list object containing the definitions of the algorithm.
#' See Section \code{Problems and Algorithms} for details.
#' @param dmax desired confidence interval halfwidth for the estimated
#'          mean/median performance of algorithm \code{algo} on instance
#'          \code{problem}.
#' @param nstart initial number of algorithm runs.
#'      See \code{Initial Number of Observations} for details.
#' @param nmax maximum allowed sample size.
#' @param seed seed for the random number generator
#'          (\code{NULL} for using \code{Sys.time()}).
#' @param ... further parameters to be passed on to \code{boot}
#'          (if method == "boot")
#'
#' @return a list object containing the following items:
#' \itemize{
#'    \item \code{x} - vector of observed performance values
#'    \item \code{x.est} - estimated value for the statistic of interest of the
#'          performance of \code{algo} on \code{problem}
#'    \item \code{se} - standard error of the estimate
#'    \item \code{delta} - the CI halfwidth
#'    \item \code{n} - number of observations generated
#'    \item \code{seed} - the seed used for the PRNG
#' }
#'
#' @author Felipe Campelo (\email{fcampelo@@ufmg.br}),
#'         Fernanda Takahashi (\email{fernandact@@ufmg.br})
#'
#' @section References:
#' P. Mathews.
#' "Sample size calculations: Practical methods for engineers and scientists".
#' Mathews Malnar and Bailey, 2010.
#' J. Botella, C. Ximenez, J. Revuelta, M. Suero.
#' "Optimization of sample size in controlled experiments: the CLAST rule".
#' Behavior Research Methods, 38(1) 65 - 76, 2006
#'
#' @examples
#' instance <- list(name = "dummyinstance", xmax = 1, xmin = 0)
#' algo     <- list(name              = "dummyalgo",
#'                  distribution.fun  = "rexp",
#'                  distribution.pars = list(rate = 0.5))
#' out      <- run_nreps(instance, algo, dmax = 1,
#'                       stat = "mean", method = "boot")
#' @export

run_nreps <- function(instance,                    # instance parameters
                      algo,                        # algorithm parameters
                      dmax,                        # desired (max) CI halfwidth
                      #stat   = c("mean", "median"),# statistic to use
                      stat   = "mean",
                      #method = c("param", "boot", "binom"), # technique to calculate CI
                      method = "param",
                      alpha  = 0.05,               # significance level for CI
                      nstart = 20,                 # initial number of samples
                      nmax   = Inf,                # maximum allowed sample size
                      seed   = NULL,               # seed for PRNG
                      ...)                         # parameters for boot
{

  # ========== Error catching ========== #
    # If the calling function was run_chase(), then skip the testing
    # since it has already been done
    if(!exists("run_chase.errorckeck.performed", envir = parent.frame(1))){
        assertthat::assert_that(
            is.list(instance), is.list(algo),
            all(assertthat::has_name(instance, c("xmax", "xmin", "name"))),
            assertthat::has_name(algo, "name"),
            is.numeric(instance$xmin), is.numeric(instance$xmax),
            length(instance$xmin) == length(instance$xmax),
            all(instance$xmin < instance$xmax),
            is.numeric(alpha), alpha > 0, alpha < 1,
            assertthat::is.count(nstart),
            is.infinite(nmax) || assertthat::is.count(nmax),
            nmax > nstart,
            is.null(seed) || assertthat::is.count(seed))
}
  # ==================================== #

  # set PRNG seed
  if(is.null(seed)) {seed <- as.integer(Sys.time())}
  set.seed(seed)

  # initial definitions
  n     <- nstart - 1       # initial number of observations
  delta <- Inf              # initial CI halfwidth

  # initialize vector of observations
  x <- get_observations(algo, instance, n)

  # Iterative cycle
  while(delta > dmax && n <= nmax){
    # Generate a new observation
    x       <- c(x, get_observations(algo, instance, 1))
    n       <- n + 1

    # Calculate CI halfwidth
    CI      <- calc_ci(x,
                       stat   = stat,
                       method = method,
                       alpha  = alpha,
                       ...)
    delta   <- diff(CI$ci) / 2
  }

  output<-list(x      = x,
               x.est  = CI$est,
               x.CI   = CI$ci,
               n      = n,
               se     = CI$se,
               delta  = delta,
               seed   = seed,
               stat   = stat,
               method = method)

  return(output)
}
