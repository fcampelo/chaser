#' Envelope Routine
#'
#' Call the sampling routine for a given set o algorithms and problems
#'
#' @param Alg list of algorithms - list()
#' @param Prob list of problems, as lists of parameters of Alg - list(list())
#' @param e desired standard error - positive numeric
#' @param n0 initial sample size - positive integer
#' @param a desired significance level (alpha) - positive numeric
#'
#' @return list with \code{list(db.full, db.short)} containig a database with all the generated data and a collapsed database with the sample size, mean and standart deviation - list of databases
#'
#' @export

envelope <- function(Alg, Prob, e, n0, a){
  # making list with optional parameters
  param <- list()
  if (!missing(e)){
    param <- c(param, e=e)
  }
  if (!missing(n0)){
    param <- c(param, n0=n0)
  }
  if (!missing(a)){
    param <- c(param, a=a)
  }

  N <- length(Alg); #number of tested algorithms
  K <- length(Prob); #number of tested problems
  #collapsed database with only one entry for every AlgxProb pair
  data.size <- N*K;
  df_short <- data.frame(Alg = numeric(data.size), Prob= numeric(data.size), n= numeric(data.size), x_mean= numeric(data.size), x_sd= numeric(data.size), stringsAsFactors=FALSE);
  #complete database of variable containg all the values generated
  df_full <- data.frame();

  l <- 1; #counter to itarate on the collapsed database
  for (i in 1:N){
    for (j in 1:K){
      #values <- sampling(Alg(i), Prob(j), e=e, n0=n0, a=a)
      input <- c(Alg[i], Prob[j], param)
      values <- do.call(sampling, input)
      #adding entry to collapsed database of fixed size
      df_short$Alg[l] <- i
      df_short$Prob[l] <- j
      df_short$n[l] <- length(values);
      df_short$x_mean[l] <- mean(values);
      df_short$x_sd[l] <- sd(values)
      l <- l+1;
      #adding n entries to complete database, n is variable
      df_full <- rbind(df_full, cbind(i, j, values))
    }
  }
  colnames(df_full) <- c("Alg", "Prob", "value") #for some reason the rbind was removing the names of the database, so added it last

  envelope <- list(df.full = df_full, df.short = df_short)
}
