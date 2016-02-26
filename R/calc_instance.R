#' Calcultes the number os instances
#'
#'
#'
#' @export



calc_instance <- function(ninstances=NULL,           # number of instances
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
        (is.null(ninstances)&& !is.null(cpower)&& !is.null(d)) || (!is.null(ninstances)&& is.null(cpower)|| is.null(d)),
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
