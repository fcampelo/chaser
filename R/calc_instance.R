#' Power calculations for the CHASE procedure
#'
#' Compute the number of instances required for a given comparison of algorithms,
#' or determine the power curve of the comparision.
#'
#' This routine calculates the number of instances necessary for a comparison
#' of \code{K} algorithms, considering a desired power level and a given value
#' for the minimally interesting effect size (either standardized in terms of
#' Cohen's d coefficient or given as an absolute value, in which case an
#' estimate for the residual standard deviation must be provided).
#'
#' If the number of instaces is predefined (e.g., when using standard benchmark
#' sets), this routine can return the power for a given effect size, or the
#' effect size that can be detected with a certain power.
#'
#' @section Usage:
#' This routine can be used in several ways. Five parameters are involved in
#' power calculations for statistical comparisons of algorithms:
#'
#' \itemize{
#'  \item the significance level (\code{alpha})
#'  \item the number of instances (\code{N})
#'  \item the power of the comparison (\code{cpower})
#'  \item the magnitude of the actual difference (\code{delta})
#'  \item the residual standard deviation (\code{sigma})
#' }
#'
#' The two last items in this list can often be combined into a standardized
#' effect size \code{d = delta/sigma}, that is, in a non-dimensional effect size
#' indicator that is given as the number of standard deviations.
#'
#' Depending on which inputs are provided to the routine, it performs different
#' calculations. The most common usages are listed below:
#'
#' \strong{Case 1:} \code{N = NULL}, with given \code{d}: calculates the
#' required number of instances, \code{N}
#'
#' \strong{Case 2:} \code{N = NULL}, with given \code{delta} and \code{sigma}:
#' since \code{d = delta / sigma}, this is equivalent to \strong{Case 1})
#'
#' \strong{Case 3:} \code{d = delta = sigma = NULL}: calculates the smallest
#' standardized effect size \code{d} that can be detected with power
#' \code{cpower}.
#'
#' \strong{Case 4:} \code{cpower = NULL}: calculates the power of the comparison
#' to detect a standardized effect size of \code{d} or greater.
#'
#' \strong{Case 5:} \code{d = cpower = NULL}: calculates the power curves
#' \code{d X cpower} of the comparison.
#'
#' \strong{Case 6:} \code{N = cpower = NULL}: calculates the power curves
#' \code{N X cpower} of the comparison.
#'
#' \strong{Case 7:} \code{N = d = NULL}: calculates the power curves
#' \code{N X d} of the comparison.
#'
#' @section MHT correction:
#' The power formulas used in this routine are based on the pairwise tests
#' (e.g., paired t-tests), and as such... [CONTINUE HERE]
#'
#'
#' @param N the number of instances to be used in the experiment.
#'    See \code{Number of Instances} for details.
#' @param cpower power deisred for the comparision
#'    See \code{Number of Instances} for details.
#' @param d desired effect size normalized by the standart deviation
#'    See \code{Number of Instances} for details.
#' @param algorithms the number of algorithms to be comparred or a list object
#'    containing lists defining all algorithms to
#'    be used in the experiment
#' @param alpha significance level for the confidence intervals on the means of
#'    each algo-problem pair.
#' @param alpha.correction is the type of alpha correction used - default is
#'    \code{holm}
#' @param comparetype defines if the comparition is 'all vc all' or 'one vs all'-
#'    default is \code{one}
#' @param direction is the direction of the distribution 1 for 1-sided and 2 for
#'    2-sided - default is \code{2}
#' @param nmax maximum allowed number of instances
#'
#' @return a list object containing the following items:
#' \itemize{
#'    \item \code{N} - number of instances
#'    \item \code{cpower} - the power of the comparision
#'    \item \code{d} - the effect size
#'    \item \code{corrected.alpha} - the corrected alpha
#' }
#' The output list also includes some of the input parameters, for convenience:
#' algorithms, alpha and nmax as given
#'
#' @author  Felipe Campelo (\email{fcampelo@@ufmg.br}),
#'          Fernanda Takahashi (\email{fernandact@@ufmg.br})
#'
#' @export



calc_instance <- function(N=NULL,           # number of instances
                          cpower=NULL,               # power
                          d=NULL,                    # proportion of standard deviations
                          algorithms,                # list of algorithms
                          alpha = 0.05,              # significance level for CI
                          alpha.correction = 'holm', #
                          comparetype = 'one',       # 'one' for one vs. all and 'all' for all vs. all
                          direction = 2,             # 1 for one sided and 2 for 2 sided
                          nmax = Inf                 # maximum allowed sample size
)
{

    assertthat::assert_that(
        (is.null(N)&& !is.null(cpower)&& !is.null(d)) || (!is.null(N)&& is.null(cpower)|| is.null(d)),
        is.list(algorithms)|| assertthat::is.count(algorithms),
        is.numeric(alpha) && alpha > 0 && alpha < 1,
        is.infinite(nmax) || assertthat::is.count(nmax)
    )
    #run_chase(instances = instance, algorithms = algo, dmax = 1)
    if (is.list(algorithms)){
        nalg <- length(algorithms)
    }else{
        nalg <- algorithms
    }
    if (tolower(comparetype == "one")){
        k.correction <- nalg-1
    }else{
        k.correction <- nalg*(nalg-1)/2
    }
    if (tolower(alpha.correction) == "holm"){
        corrected.alpha <-  (1-(1-alpha)^(k.correction))/direction
    }
    if (is.null(N)){
        t_b <- qnorm(cpower)
        t_a <- qnorm(1-corrected.alpha)
        N <- ceiling(((t_b + t_a)*(1/d))^2)
        while ((t_a <= t_b)&&(N < nmax)){
            N <- N+1
            t_a <- qt(1-corrected.alpha,N-1)
            t_b <- qt(cpower, N-1)
        }
    }
    else{
        if(is.null(cpower)&&is.null(d)){
            t_a <- qt(1-corrected.alpha,N-1)
            cpower <- c(0.1, 0.3, 0.5, 0.7, 0.9)
            d <- sapply(cpower, function(x) ((qt(x, N-1) +t_a)/sqrt(N)))
        }
        else{
            if(is.null(cpower)){
                t_b <- (d) * sqrt(N) - qt(1-corrected.alpha, N-1)
                cpower <- pt(t_b, N-1)
            }
            else{
                t_a <- qt(1-corrected.alpha,N-1)
                t_b <- qt(cpower, N-1)
                d <- ((t_b +t_a)/sqrt(N))
            }
        }
    }
    #print(N)
    #print(cpower)
    #print(d)
    output<-list(N     = N,
                 cpower         = cpower,
                 d              = d,
                 corrected.alpha = corrected.alpha,
                 algorithms   = algorithms,
                 alpha        = alpha,
                 nmax         = nmax)
    return(output)
}
