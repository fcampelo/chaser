#' Calcultes the number of instances
#'
#' Return either the number of isntances, the effect size or the power of the
#' comparision considering the other two
#'
#' This routine uses the closed formula of the t-test to calculates the number
#' of instances necessary for a comparition of k algorithms, considering given
#' power and  effect size.
#' In case the number of instaces is already known, it can return the effect
#' size or power, given the other
#'
#' @section Number of Instances:
#'
#'
#' @param ninstances the number of instances to be used in the experiment.
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
#' @param compare.type defines if the comparition is 'all vc all' or 'one vs all'-
#'    default is \code{one}
#' @param direction is the direction of the distribution 1 for 1-sided and 2 for
#'    2-sided - default is \code{2}
#' @param nmax maximum allowed number of instances
#'
#' @return a list object containing the following items:
#' \itemize{
#'    \item \code{ninstances} - number of instances
#'    \item \code{cpower} - the power of the comparision
#'    \item \code{d} - the effect size
#'    \item \code{corrected.aplha} - the corrected alpha
#'    and some input parameters: algorithms, alpha and nmax as given
#' }
#'
#' @author Felipe Campelo (\email{fcampelo@@ufmg.br}), Fernanda Takahashi (\email{fernandact@@ufmg.br})
#'
#' @export



calc_instance <- function(ninstances=NULL,           # number of instances
                          cpower=NULL,               # power
                          d=NULL,                    # normalized standard deviations
                          algorithms,                # list or number of algorithms
                          alpha = 0.05,              # significance level for CI
                          alpha.correction = 'holm', #
                          compare.type = 'one',       # 'one' for one vs. all and 'all' for all vs. all
                          direction = 2,             # 1 for one sided and 2 for 2 sided
                          nmax = Inf                 # maximum allowed sample size
)
{

    assertthat::assert_that(
        (is.null(ninstances)&& !is.null(cpower)&& !is.null(d)) || (!is.null(ninstances)&&assertthat::is.count(ninstances)&& is.null(cpower)|| is.null(d)),
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
    if (tolower(compare.type == "one")){
        k.correction <- nalg-1
    }else{
        k.correction <- nalg*(nalg-1)/2
    }
    if (tolower(alpha.correction) == "holm"){
        corrected.alpha <-  (1-(1-alpha)^(k.correction))/direction
    }
    if (is.null(ninstances)){
        t_b <- qnorm(cpower)
        t_a <- qnorm(1-corrected.alpha)
        ninstances <- ceiling(((t_b + t_a)*(1/d))^2)
        while ((t_a <= t_b)&&(ninstances < nmax)){
            ninstances <- ninstances+1
            t_a <- qt(1-corrected.alpha,ninstances-1)
            t_b <- qt(cpower, ninstances-1)
        }
    }
    else{
        if(is.null(cpower)&&is.null(d)){
            t_a <- qt(1-corrected.alpha,ninstances-1)
            cpower <- c(0.1, 0.3, 0.5, 0.7, 0.9)
            d <- sapply(cpower, function(x) ((qt(x, ninstances-1) +t_a)/sqrt(ninstances)))
        }
        else{
            if(is.null(cpower)){
                t_b <- (d) * sqrt(ninstances) - qt(1-corrected.alpha, ninstances-1)
                cpower <- pt(t_b, ninstances-1)
            }
            else{
                t_a <- qt(1-corrected.alpha,ninstances-1)
                t_b <- qt(cpower, ninstances-1)
                d <- ((t_b +t_a)/sqrt(ninstances))
            }
        }
    }
    #print(ninstances)
    #print(cpower)
    #print(d)
    output<-list(ninstances     = ninstances,
                 cpower         = cpower,
                 d              = d,
                 corrected.alpha = corrected.alpha,
                 algorithms   = algorithms,
                 alpha        = alpha,
                 nmax         = nmax)
    return(output)
}
