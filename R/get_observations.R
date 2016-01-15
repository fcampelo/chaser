#' Run algorithm on problem
#'
#' Call algorithm routine for the solution of a problem instance
#'
#'
#' @param instance a list object containing the definitions of the problem
#'    instance. See \code{\link{run_nreps}} for details.
#' @param algo a list object containing the definitions of the algorithm.
#'    See \code{\link{run_nreps}} for details.
#' @param n number of observations to generate.
#'
#' @return vector of observed performance values
#'
#' @seealso \link{run_nreps}
#' @author Felipe Campelo (fcampelo@@ufmg.br)
#'

get_observations <- function(algo,        # algorithm parameters
                             instance,    # problem parameters
                             n = 1)       # number of observations to generate.
{

  # ========== Error catching ========== #
  # Most of error catching is already performed by the calling routine
  # run_nreps(), so no need to do much here, except:
  assertthat::assert_that(assertthat::is.count(n))
  # ==================================== #


  # remove '$name' field from list of arguments
  # and include the problem definition as field 'probpars'
  myargs          <- algo[names(algo) != "name"]
  myargs$instance <- instance

  # Get observation(s)
  f <- numeric(n)
  for (i in 1:n){
    result <- do.call(algo$name,
                      myargs)
    f[i] <- result$Fbest
  }
  return(f)
}
