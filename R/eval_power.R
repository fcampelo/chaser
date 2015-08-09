#' Post-experimental power analysis
#'
#' Power analysis of the comparison of mean performance of two algorithms on N
#' problem instances.
#'
#' The details of this routine come here...
#'
#'
#' @param x
#'
#' @return y
#'
#' @author Felipe Campelo
#' @export

# Example (temporary - incorporate to documentation later)
# data<-read.table("data/algo.txt",header = TRUE,sep=",")
# data <- subset(data,Algorithm%in%c("Alg.1","Alg.2"))
# data$Algorithm<-factor(data$Algorithm)

# TODO: - replace power.t.test by explicit (exact) formula
#       - decide how to outout results (initial version: scalar. Afterwards: vector+ggplot)

eval_power <- function(data,            # data frame containing the experimental data
                       alpha = 0.05,    # significance level
                       delta = NA,      # minimal difference of practical interest
                       power = NA,      # desired power level
                       alternative = "two.sided"  # alternative hypothesis
                       #method = c("param", "bootstrap") # type of test
                       #nboot = 1000,   # number of bootstrap samples
                       #ncpus = 1       # number of cores to use for bootstrap
){
    # Step -1: initial setup
    data.names  <- c("Algorithm", "Instance", "Result")
    names(data) <- data.names
    data[, 1]   <- as.factor(data[, 1])
    data[, 2]   <- as.factor(data[, 2])
    nProbs      <- nlevels(data[, 2])

    # Step 0: error catching
    if (nlevels(data[, 1]) != 2 | nProbs < 2){
        stop("eval_power() requires a data frame with exactly two algorithms
             and more than two problem instances")
    }

    if (!is.na(delta)){
        if (delta <= 0) {stop("eval_power() requires delta to be positive")}
    }


    # Step 1: aggregate data as mean values
    aggdata <- with(data,
                    aggregate(x   = Result,
                              by  = list(Algorithm, Instance),
                              FUN = mean))
    names(aggdata) <- data.names

    # Step 2: Perform significance test
    signif.test <- t.test(Result~Algorithm,
                          data = aggdata,
                          paired = TRUE,
                          alternative = alternative,
                          conf.level = 1 - alpha)

    # Step 3: power analysis
    # 3.1 - paired differences: confidence interval on sd
    pds <- with(aggdata,
                Result[Algorithm == levels(Algorithm)[1]] -
                    Result[Algorithm == levels(Algorithm)[2]])
    sd.est  <- sd(pds)
    sd.ci   <- sqrt((nProbs-1)*sd.est^2 / qchisq(c(1-alpha/2,alpha/2), nProbs-1))

    # 3.2 - power
    power.est <- power.t.test(n           = nProbs,
                              sd          = sd.est,
                              sig.level   = alpha,
                              power       = power,
                              type        = "paired",
                              alternative = alternative)



#     # Step 3: power analysis
#     # 3.1 - paired differences: confidence interval on sd
#     pds <- with(aggdata,
#                 Result[Algorithm == levels(Algorithm)[1]] -
#                     Result[Algorithm == levels(Algorithm)[2]])
#     sd.est  <- sd(pds)
#     sd.ci   <- sqrt((nProbs-1)*sd.est^2 / qchisq(c(1-alpha/2,alpha/2), nProbs-1))
#
#     # 3.2 Generate power curve
#     power.seq <- (1:99)/100
#     delta.fun <- function(x){
#         sapply(X = power.seq,
#                FUN = function(x,...){return(power.t.test(power = x, ...)$delta)},
#                n = nProbs,
#                sd = x,
#                sig.level = alpha,
#                type = "paired",
#                alternative = alternative,
#                simplify = TRUE)
#     }
#
#     delta.est <- delta.fun(sd.est)
#
#     sapply(X = (1:10)/100,
#            FUN = function(x,...){power.t.test(power=x,...)$delta},
#            n = 30, sd = 30, type = "paired",
#            simplify = TRUE)
#
#
#
#     # Plot delta x power curves, highlight specified powers
}
