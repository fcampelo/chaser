Alg <- list(list(name = "dummyalgo",
                 distribution.fun = "rnorm",
                 distribution.pars = list(mean = 0, sd = 1)),
            list(name = "dummyalgo",
                 distribution.fun = "rnorm",
                 distribution.pars = list(mean = 0, sd = 2)),
            list(name = "dummyalgo",
                 distribution.fun = "rnorm",
                 distribution.pars = list(mean = 0, sd = 4)),
            list(name = "dummyalgo",
                 distribution.fun = "rnorm",
                 distribution.pars = list(mean = 0, sd = 8)),
            list(name = "dummyalgo",
                 distribution.fun = "rnorm",
                 distribution.pars = list(mean = 0, sd = 16)))

Prob <- list(list(name = "dummyinstance", xmin = 0, xmax = 1))

#out <-envelope(Alg, Prob, e=1)



out <- run_chase(instances = Prob, algorithms = Alg, dmax = 1)

