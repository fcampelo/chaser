#' Sampling Routine
#'
#' Samples the number of executions needed to obtain the desired error
#'
#' @section References:
#' Botella, J., Ximenez, C., Revuelta, J. and Suero, M.,
#' "Optimization of sample size in controlled experiments: the CLAST rule.",
#' Behavior Research Methods, 38(1) 65 - 76, 2006
#'
#' @param A name of the function (algorithm) to be evaluated
#' @param P parameters  of A (problem to be evaluated)
#' @param e desired standard error
#' @param n0 initial sample size
#' @param a desired significance level (alpha)
#'
#' @return vector \code{values} containing the sampled
#'
#' @export

sampling <- function(A, P, e = 0.05, n0 = 10, a = 0.1){

  ## run algorithm the initial amount of times
  values <- rep(NA, n0)
  for (i in 1:n0){
    values[i] <- do.call(A, P);
  }

  #calculate the initial delta interval
  t_a2 = qnorm(1-a/2); #initialy uses the z distribution (why not use the t?)
  smp_sd <- sd(values)
  delta <- (t_a2 * smp_sd)/sqrt(n0);

  #similar to the CLAST
  #executes the algorithm one more time and recalculate delta
  #repeat until the desired delta is achieved
  i <-1;
  while (delta > e){
    new <- do.call(A, P)
    values<-c(values, new)
    smp_sd <- sd(values)
    n <- n0+i
    t_a2 <- qt(1-a/2,n-1)#uses the t distribution considering wich uses the number of samples
    delta <- (t_a2 * smp_sd)/sqrt(n);
    i <- i+1;
  }
  #print(delta)
  #print(n)

  #return the sampled data
  #the mean, standart deviation and delta can be calculated from
  #the given information
  sampling <- values;
}
