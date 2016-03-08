#' Power and sample size calculations for the CHASE procedure
#'
#' Compute the number of instances required for a given comparison of
#' algorithms, or determine the power curve of the comparision.
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
#' (e.g., paired t-tests), which mean that in cases where more than 2 algorithms
#' are to be compared (that is, if there will be \emph{multiple hypotheses
#' testing} - MHT) the significance level of each test must be corrected to
#' avoid inflating the familywise error rate.
#'
#' In the current version only the single-step Holm correction is considered in
#' the calculation of the sample size. However, if the analysis is eventually
#' performed using the Bonferroni correction the loss of power is in most cases
#' quite small.
#'
#'
#' @section Types of comparisons:
#' Depending on the type of comparison being planned, a different set of tests
#' can be performed: for instance, for a general comparison of methods (in which
#' there is no algorithm that is of particular interest to the researcher) an
#' "all-vs-all" comparison (\code{comparetype = "all-vs-all"}) should be
#' performed, to obtain a complete ordering of performance among the
#' alternatives. In this case \code{K * (K - 1) / 2} comparisons are performed.
#' On the other hand, if there is a single method in which the researcher is
#' particularly interested in (usually the "proposed algorithm") an "all-vs-one"
#' comparison can be used (\code{comparetype = "one-vs-all"}), since in this
#' case the question of interest is how one particular algorithm compares
#' against the others. In this case the total number of comparisons is
#' \code{K - 1}, which means that the MHT correction of \code{alpha} is smaller,
#' resulting in a smaller sample size required (or a higher statistical power
#' for a predefined sample size).
#'
#' @section Directionality of the alternative hypothesis:
#' In the case of an "all-vs-one" comparison, the researcher may be interested
#' in two kinds of alternative hypotheses: \code{direction = "two-sided"},
#' if he wants to infer whether the "proposed algorithm" is \emph{different}
#' from each of the competing approaches; or \code{direction = "one-sided"}, if
#' he is only interested in determining whether it is \emph{better} than the
#' others. In the former, the statistical question is
#' "\emph{is the method better or worse than the others, or is it not
#' different?}"; in the later, "\emph{is the method better han the others, or is
#' it not?}". Using a one-sided alternative results in a smaller sample size
#' required (or a higher statistical power for a predefined sample size).
#'
#' @section References:
#' S. Holm.
#' "A simple sequentially rejective multiple test procedure".
#' Scandinavian Journal of Statistics (6):65-70, 1979.
#'
#' P. Mathews.
#' "Sample size calculations: Practical methods for engineers and scientists".
#' Mathews Malnar and Bailey, 2010.
#'
#' Y. Benjamini, Y.,Hochberg.
#' "Controlling the false discovery rate: a practical and powerful approach to
#' multiple testing".
#' Journal of the Royal Statistical Society Series B 57:289-300, 1995.
#'
#' Y. Benjamini and D. Yekutieli
#' "The control of the false discovery rate in multiple testing under
#' dependency".
#' Annals of Statistics 29:1165-1188, 2001.
#'
#' @param N number of instances to be used in the experiment.
#' @param cpower desired power for the comparison.
#' @param d minimally interesting effect size - standardized (see \code{Usage}
#'          for details).
#' @param delta minimally interesting effect size - absolute value (see
#'          \code{Usage} for details).
#' @param sigma expected residual standard deviation (see \code{Usage} for
#'          details).
#' @param algorithms either an integer representing the number of algorithms to
#'          be compared, or a list object containing the definitions of the
#'          algorithms to be compared (as defined in \code{\link{run_chase}})
#' @param alpha significance level (familywise)
#' @param mht.correction method used for correcting the value of \code{alpha}
#'          for multiple hypotheses testing (see \code{MHT correction} for
#'          details).
#' @param comparetype type of comparison planned (see \code{Types of
#'          comparisons} for details).
#' @param direction type of alternative hypotheses planned (see
#'          \code{Directionality of the alternative hypothesis} for details).
#' @param max.instances maximum allowed number of instances.
#'
#' @return List object containing the input fields, plus the calculated values
#'          as described in section \code{Usage}.
#'
#' @author  Felipe Campelo (\email{fcampelo@@ufmg.br}),
#'          Fernanda Takahashi (\email{fernandact@@ufmg.br})
#'
#' @export

chase_sample_size <- function(N              = NULL,
                              cpower         = NULL,
                              d              = NULL,
                              delta          = NULL,
                              sigma          = NULL,
                              algorithms,
                              alpha,
                              mht.correction = "holm",
                              comparetype    = c("one-vs-all", "all-vs-all"),
                              direction      = c("one-sided" , "two-sided"),
                              max.instances  = Inf)
{
    # ========== Error catching ========== #
    # If d is NULL and (delta + sigma) are defined then define d = delta/sigma
    if(is.null(d) & !any(sapply(list(delta, sigma), is.null))) d <- delta/sigma

    # if "algorithms" is given as a list use the list length
    if(is.list(algorithms)) nalg <- length(algorithms) else nalg <- algorithms

    # Check consistency of the inputs
    assertthat::assert_that(
        sum(sapply(list(N, cpower, d), is.null) <= 2),
        is.null(N) || assertthat::is.count(N),
        is.null(cpower) || (assertthat::is.number(cpower) && cpower > 0 && cpower < 1),
        is.null(d) || (assertthat::is.number(d) && d > 0),
        assertthat::is.count(nalg) && nalg > 1,
        is.numeric(alpha) && alpha > 0 && alpha < 1,
        is.character(mht.correction) && length(mht.correction) == 1,
        mht.correction %in% c("holm"),
        is.character(comparetype) && length(comparetype) == 1,
        comparetype %in% c("one-vs-all", "all-vs-all"),
        is.character(direction) && length(direction) == 1,
        direction %in% c("one-sided" , "two-sided"),
        assertthat::are_equal(direction, "two-sided") || assertthat::are_equal(comparetype, "one-vs-all"),
        is.infinite(max.instances) || assertthat::is.count(max.instances))
    # ==================================== #

    # ========== Standard defitinions ========== #
    if (comparetype == "one-vs-all"){
        k.correction <- nalg-1
    } else{
        k.correction <- nalg * (nalg - 1) / 2
    }
    if (direction == "two-sided") dir <- 2 else dir <- 1


    if (mht.correction == "holm"){
        corrected.alpha <-  (1 - (1 - alpha) ^ k.correction) / dir
    }

    #==== STOPPED HERE

    if (is.null(N)){
        if (is.null(cpower)){ # Return N x cpower curves


        } else if (is.null(d)){ # Return N x d curves


        } else { # Calculate N
            # [FIXME: CORRECT NONCENTRALITY PARAMETER]
            N   <- 2                           # 2 is the minimal sample size
            t_b <- qt(p = cpower, df = N - 1, ncp = d)
            t_a <- qt(p = 1 - corrected.alpha, df = N - 1)
            while ((t_a > t_b) && (N <= max.instances)){
                N <- N + 1
                t_b <- qt(p = cpower, df = N - 1, ncp = d)
                t_a <- qt(p = 1 - corrected.alpha, df = N - 1)
            }
        }

    }
    else{
        if(is.null(cpower) && is.null(d)){ # Return d x cpower curves
            t_a    <- qt(1 - corrected.alpha, N - 1)
            cpower <- seq(from = 0.1, to = 0.9, by = 0.05)
            d      <- sapply(cpower, function(x) ((qt(x, N - 1) + t_a) / sqrt(N)))
        } else if(is.null(cpower)){ # Calculate power
            t_b    <- (d) * sqrt(N) - qt(1 - corrected.alpha, N - 1)
            cpower <- pt(t_b, N - 1)
        } else if(is.null(d)){ # Calculate minimal detectable effect size at power 'cpower'
            t_a    <- qt(1 - corrected.alpha, N - 1)
            t_b    <- qt(cpower, N - 1)
            d      <- ((t_b + t_a) / sqrt(N))
        } else stop("All parameters provided to chase_sample_size()\nNothing to calculate")
    }

    output<-list(N               = N,
                 cpower          = cpower,
                 d               = d,
                 corrected.alpha = corrected.alpha,
                 algorithms      = algorithms,
                 alpha           = alpha,
                 nmax            = nmax)
    return(output)
}
