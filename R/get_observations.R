#' Run algorithm on problem
#'
#' Call algorithm routine for the solution of a problem instance
#'
#' The details of this routine come here...
#'
#' @section Problem structure:
#' The \code{problem} list must contain all relevant parameters that define the
#' problem instance. \code{problem} is a named list containing at least the following
#' fields:
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
#' @param problem a list object containing the definitions of the problem
#'    instance. See \code{\link{Problem structure}} for details.
#' @param algo a list object containing the definitions of the algorithm See
#'    \code{\link{Algorithm structure}} for details.
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
#' @author Felipe Campelo

get_observations <- function(algo,        # algorithm parameters
                             problem,     # problem parameters
                             n = 1)       # number of observations to generate.
{
      result <- do.call(algo$name,
                        algo$pars,
                        problem)
      return(result$f)
}
