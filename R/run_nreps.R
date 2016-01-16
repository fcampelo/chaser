#' Automatically determine sample size for an algorithm on a problem instance
#'
#' Iteratively calculates the required sample size for an algorithm on a
#'    problem instance, so that the final number of replicates yields an
#'    estimation of mean performance with a predefined maximum uncertainty.
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
#'
#' @param instance a list object containing the definitions of the problem
#'    instance. See Section \code{Problems and Algorithms} for details.
#' @param algo a list object containing the definitions of the algorithm.
#' See Section \code{Problems and Algorithms} for details.
#' @param dmax desired confidence interval halfwidth for the estimated mean
#'    performance of algorithm \code{algo} on instance \code{problem}.
#' @param alpha significance level for the confidence intervals.
#'      Defaults to \code{alpha = 0.05}.
#' @param nstart initial number of algorithm runs.
#'      Defaults to \code{nstart = 10}.
#' @param nmax maximum allowed sample size.
#'      Defaults to \code{nmax = Inf}.
#' @param seed seed for the random number generator.
#'      Defaults to \code{seed = NULL}.
#'
#' @return a list object containing the following items:
#' \itemize{
#'    \item \code{x} - vector of observed performance values
#'    \item \code{xbar} - sample mean of the performance of \code{algo} on \code{problem}
#'    \item \code{se} - standard error of the mean
#'    \item \code{delta} - the CI halfwidth
#'    \item \code{n} - number of observations
#'    \item \code{seed} - the seed used for the PRNG
#' }
#'
#' @author Fernanda Takahashi (\email{fernandact@@ufmg.br}), Felipe Campelo (\email{fcampelo@@ufmg.br})
#'
#' @section References:
#' Paul Mathews. "Sample size calculations: Practical methods for engineers and
#' scientists". Mathews Malnar and Bailey, 2010.
#' Botella, J., Ximenez, C., Revuelta, J. and Suero, M.,
#' "Optimization of sample size in controlled experiments: the CLAST rule.",
#' Behavior Research Methods, 38(1) 65 - 76, 2006
#'
#' @export

run_nreps <- function(instance,       # instance parameters
                      algo,           # algorithm parameters
                      dmax,           # desired (maximum) CI halfwidth
                      alpha = 0.05,   # significance level for CI
                      nstart = 10,    # initial number of samples
                      nmax = Inf,     # maximum allowed sample size
                      seed = NULL     # seed for PRNG
){

  # ========== Error catching ========== #
    # If the calling function was run_chase(), then skip the testing
    # since it has already been done
    if(!exists("run_chase.errorckeck.performed", envir = parent.frame(1))){
        assertthat::assert_that(
            is.list(instance) && is.list(algo),
            all(assertthat::has_name(instance, c("xmax", "xmin", "name"))),
            assertthat::has_name(algo, "name"),
            is.numeric(instance$xmin) && is.numeric(instance$xmax),
            length(instance$xmin) == length(instance$xmax),
            all(instance$xmin < instance$xmax),
            is.numeric(alpha) && alpha > 0 && alpha < 1,
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
                       FUN    = "mean",
                       alpha  = alpha,
                       method = "parametric")
    delta   <- diff(CI$ci) / 2
  }

  output<-list(x    = x,
               xbar = mean(x),
               n    = n,
               se   = CI$se,
               delta = delta,
               seed = seed)

  return(output)
}
