#' Run algorithm on problem
#'
#' Call algorithm routine for the solution of a problem instance
#'
#'
#'
#' @param problem a list object containing the definitions of the problem
#'    instance. See \code{\link{run_chase}} for details.
#' @param algo a list object containing the definitions of the algorithm.
#'    See \code{\link{run_chase}} for details.
#' @param n number of observations to generate.
#'
#' @return a list object containing the following items:
#' \itemize{
#'    \item \code{xbar} - sample mean of the performance of \code{algo} on \code{problem}
#'    \item \code{s} - sample standard deviation of the performance of \code{algo} on \code{problem}
#'    \item \code{n} - number of observations
#'    \item \code{delta} - half-width of the confidence interval obtained
#'    \item \code{ci} - confidence interval
#'    \item \code{x} - vector of observed performance values
#'    \item \code{alpha} - confidence level of the interval
#' }
#'
#' @seealso \link{run_chase}
#' @author Felipe Campelo

get_observations <- function(algo,        # algorithm parameters
                             problem,     # problem parameters
                             n = 1)       # number of observations to generate.
{

    # remove '$name' field from list of arguments
    # and include the problem definition as field 'probpars'
    myargs          <- algo[names(algo) != "name"]
    myargs$probpars <- problem

    # Get observation
    f <- numeric(n)
    for (i in 1:n){
        result <- do.call(algo$name,
                          myargs)
        f[i] <- result$Fbest
    }
    return(f)
}
