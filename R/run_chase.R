#' Comparison of Heuristics with Adaptive Sample size Estimation
#'
#' Run experiment using adaptive sample size estimation.
#'
#' This routine executes a full experiment using the CHASE method. It
#' essentially runs each algorithm (i) on each problem instance (j) until the
#' uncertainty of the estimate of the mean performance is below a given
#' threshold.
#'
#' The result is an unbalanced experiment, in which the within-instance sample
#' size is proportional to the variance of the algorithm on that instance, so
#' that it yields an estimate of mean performance with fixed "measurement
#' error".
#'
#' Based on this "measurement error" (which is defined in advance) it is
#' possible to use standard formulas to calculate the effective sample size
#' (i.e., the number of problem instances required to achieve a certain
#' statistical power for a desired effect size) for a given comparative
#' experiment with metaheuristics (or other algorithms).
#'
#' @section Instances and Algorithms:
#' Parameters \code{instances} and \code{algorithms} must lists of
#' \code{instance} (\code{algorithm}) specifications, each defined according to
#' the following instructions:
#'
#' Each component of \code{instances} is itself a list containing all relevant
#' parameters that define the problem instance. Each element
#' \code{instances[[i]]} must be a named list as defined in the documentation
#' of function \code{\link{run_nreps}}.
#'
#' Similarly, the \code{algorithms} parameter must contain a list in which each
#' element provides the full specification of an algorithm, i.e., a named list
#' as defined in the documentation of function \code{\link{run_nreps}}
#'
#' @section Random Seed:
#' The \code{seed} argument receives the desired seed for the PRNG. This value
#' can be set for reproducibility purposes. The value of this parameter defaults
#' to NULL, in which case the seed is arbitrarily set using
#' \code{as.numeric(Sys.time())}. The value used as seed is always reported in
#' the output structure.
#'
#' Notice that this \code{seed} is not passed to the algorithm routines in the
#' current version.
#'
#' @section Initial Number of Observations:
#' In the general case the initial number of observations
#' (\code{nstart}) should be relatively high (> 15 if outliers are not
#' expected, > 50 if that assumption can't be made) to guarantee good
#' statistical properties (particularly compliance with nominal type-I error
#' rate \code{alpha}). However, if some distributional assumptions can be
#' made (e.g., low skewness of the population of algorithm results on the test
#' instances), then \code{nstart} can in principle be as small as 5.
#'
#' In general, higher sample sizes are the price to pay for abandoning
#' distributional assumptions. Use lower values of \code{nstart} with caution.
#'
#' @param instances a list object containing lists defining all problem
#'    instances to be used in the experiment.
#'    See \code{Instances and Algorithms} for details.
#' @param algorithms a list object containing lists defining all algorithms to
#'    be used in the experiment.
#'    See \code{Instances and Algorithms} for details.
#' @param dmax desired confidence interval halfwidth for the estimated mean
#'    performance of each algorithm on each instance.
#' @param alpha significance level for the confidence intervals on the means of
#'    each algo-problem pair. Defaults to \code{alpha = 0.05}.
#' @param nstart initial number of algorithm runs.
#'      Defaults to \code{nstart = 15}.
#'      See \code{Initial Number of Observations} for details.
#' @param nmax maximum allowed sample size. Defaults to \code{nmax = Inf}.
#' @param seed seed for the random number generator.
#'      Defaults to NULL.
#'      See \code{Random seed} for details.
#'
#' @return a list object containing the following items:
#' \itemize{
#'    \item \code{data.raw} - data frame containing all observations generated
#'    \item \code{data.summary} - data frame containing the means, standard
#'    errors and sample sizes of each algorithm on each problem instance.
#' }
#'
#' @author Felipe Campelo (\email{fcampelo@@ufmg.br}), Fernanda Takahashi (\email{fernandact@@ufmg.br})
#'
#' @export

run_chase <- function(instances,    # list of instances
                      algorithms,   # list fo algorithms
                      dmax,         # desired (maximum) CI halfwidth
                      alpha = 0.05, # significance level for CI
                      nstart = 10,  # initial number of samples
                      nmax = Inf,   # maximum allowed sample size
                      seed = NULL   # seed for PRNG
)
{

    # ========== Error catching ========== #
    assertthat::assert_that(
        is.list(instances) && is.list(algorithms),
        is.numeric(alpha) && alpha > 0 && alpha < 1,
        assertthat::is.count(nstart),
        is.infinite(nmax) || assertthat::is.count(nmax),
        nmax > nstart,
        is.null(seed) || assertthat::is.count(seed))

    # Number of instances and algorithms
    nprobs <- length(instances)
    nalgos <- length(algorithms)

    for (i in 1:nprobs){
        assertthat::assert_that(
            all(assertthat::has_name(instances[[i]], c("xmax", "xmin", "name"))),
            is.numeric(instances[[i]]$xmin) && is.numeric(instances[[i]]$xmax),
            length(instances[[i]]$xmin) == length(instances[[i]]$xmax),
            all(instances[[i]]$xmin < instances[[i]]$xmax))
    }
    for (i in 1:nalgos){
        assertthat::assert_that(assertthat::has_name(algorithms[[i]], "name"))
    }
    # Define an absurdly-named variable to help avoiding duplicate testing
    # in run_nreps()
    run_chase.errorckeck.performed <- TRUE
    # ==================================== #

    # set PRNG seed
    if(is.null(seed)) {seed <- as.integer(Sys.time())}
    set.seed(seed)

    # Initialize data.raw dataframe with its smallest possible size,
    # i.e., nprobs * nalgos * nstart
    algonames <- unlist(lapply(algorithms, function(x) x$name))
    probnames <- unlist(lapply(instances, function(x) x$name))
    data.raw <- data.frame(Algorithm   = factor(x = character(nprobs * nalgos * nstart),
                                                levels = unique(algonames)),
                           Instance    = factor(x = character(nprobs * nalgos * nstart),
                                                levels = unique(probnames)),
                           Observation = numeric(nprobs * nalgos * nstart))

    # Keep an "empty" copy of data.raw in case we need to grow it in the
    # iterative cycle (which is almost guaranteed)
    empty.raw.df <- data.raw

    # Initialize data.summary dataframe with its exact size
    nrows.summary <- nprobs * nalgos
    data.summary <- data.frame(Algorithm   = factor(x = character(nrows.summary),
                                                    levels = unique(algonames)),
                               Instance    = factor(x = character(nrows.summary),
                                                    levels = unique(probnames)),
                               x.mean    = numeric(nrows.summary),
                               x.se      = numeric(nrows.summary),
                               x.n       = numeric(nrows.summary))

    # Iterative cycle
    rawcount <- 0
    sumcount <- 0
    for (i in 1:nalgos){
        for (j in 1:nprobs){
            # Generate the data for algorithm i on instance j
            data.ij <- run_nreps(instance = instances[[j]],
                                 algo     = algorithms[[i]],
                                 dmax     = dmax,
                                 alpha    = alpha,
                                 nstart   = nstart,
                                 nmax     = nmax,
                                 seed     = seed)

            # If needed, grow the data.raw dataframe
            nij <- data.ij$n
            if (rawcount + nij > nrow(data.raw)){
                data.raw <- rbind(data.raw, empty.raw.df)
            }

            # Update data.raw
            data.raw[(rawcount + 1):(rawcount + nij), ] <- data.frame(
                Algorithm   = rep(algorithms[[i]]$name, nij),
                Instance    = rep(instances[[j]]$name, nij),
                Observation = data.ij$x)
            rawcount <- rawcount + nij

            # Update data.summary
            sumcount <- sumcount + 1
            data.summary[sumcount, ] <- data.frame(
                Algorithm = algorithms[[i]]$name,
                Instance  = instances[[j]]$name,
                x.mean    = data.ij$xbar,
                x.se      = data.ij$se,
                x.n       = nij)
        }
    }

    # Remove any extra rows from data.raw
    data.raw <- data.raw[1:rawcount, ]

    output<-list(data.raw     = data.raw,
                 data.summary = data.summary,
                 instances    = instances,
                 algorithms   = algorithms,
                 dmax         = dmax,
                 alpha        = alpha,
                 nstart       = nstart,
                 nmax         = nmax,
                 seed         = seed)

    return(output)
}
