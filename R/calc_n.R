#' Sample size for fixed-width confidence interval for the mean
#'
#' Number of observations for obtaining a CI for the mean with predefined
#'    halfwidth
#'
#' \code{calc_n} calculates the required number of observations needed to
#'    generate a fixed-width (1 - \code{alpha})-CI for the mean based on the t
#'    distribution.
#'
#' @param delta Desired halfwidth of the confidence interval
#' @param sd Estimated standard deviation
#' @param alpha Significance level of the CI.
#' @return The required sample size for obtaining a (1 - \code{alpha})-CI with
#'    halfwidth \code{delta}
#' @author Felipe Campelo
#' @examples
#' calc_n(delta = 3, sd = 8)
#' calc_n(delta = .25, sd = 1, alpha = 0.01)

calc_n <- function(delta,           # Desired CI halfwidth
                   sd,              # Estimated standard deviation
                   alpha = 0.05)    # Significance level
{
      # Error catching
      if (delta <= 0)   {stop("CI halfwidths must be positive")}
      if (sd <= 0)      {stop("Standard deviations must be positive")}

      tq <- qnorm(1 - alpha/2)      # Initial point for the search
      n <- 0                        # Initial sample size
      while (n < (tq*sd/delta)^2){
            n     <- (tq*sd/delta)^2            # Estimate n
            tq    <- qt(1 - alpha/2, n - 1)     # Update the t-quantile
      }
      return(n)
}
