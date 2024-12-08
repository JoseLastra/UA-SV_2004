rm(list=ls())
# .libPaths("C:/GH/R-4.0.2/library")  # Gerard when working from home

install.packages("asbio")
library(asbio)

anm.mc.bvn(start = c(-4, -4), mu = c(0, 0), 
           sigma = matrix(2, 2, data = c(1, 0, 0, 1)),
           length = 1000, sim = "M", jump.kernel = 0.2, 
           xlim = c(-4, 4), ylim = c(-4, 4), interval = 0.01,
           show.leg = TRUE)

anm.mc.bvn(start = c(-4, -4), mu = c(0, 0), 
           sigma = matrix(2, 2, data = c(1, 0, 0, 1)),
           length = 1000, sim = "M", jump.kernel = 0.5, 
           xlim = c(-4, 4), ylim = c(-4, 4), interval = 0.01,
           show.leg = TRUE)

anm.mc.bvn(start = c(4, -4), mu = c(0, 0), 
           sigma = matrix(2, 2, data = c(1, 0.9, 0.9, 1)),
           length = 1000, sim = "M", jump.kernel = 0.2, 
           xlim = c(-4, 4), ylim = c(-4, 4), interval = 0.01,
           show.leg = TRUE)

anm.mc.bvn(start = c(0, 0), mu = c(0, 0), 
           sigma = matrix(2, 2, data = c(1, 0, 0, 1)),
           length = 10000, sim = "M", jump.kernel = 0.5, 
           xlim = c(-4, 4), ylim = c(-4, 4), interval = 0.01,
           show.leg = TRUE)
