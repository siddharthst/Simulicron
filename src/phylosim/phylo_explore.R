#!/usr/bin/env Rscript

source("./phylosim.R")


library(parallel)
mc.cores <- detectCores() - 1

maxG <- 200
maxTEpop <- 200
n0 <- 1
maxreps <- 10000
out.file <- "explore.txt"

cat("u\tv\tk.regul\tk.pi\tdeg.rate\tdeg.effect\tsize\tdepth\n", file=out.file)

u.range        <- c(0.5,0.5)
v.range        <- c(0.05,0.05) 
deg.rate.range <- c(0.02, 0.02)
deg.effect.range <- c(0.1, 0.1)
k.regul.range  <- c(0, 0.1)
k.pi.range     <- c(0, 0.05) 

mclapply(1:maxreps, function(r) {
	u          <- runif(1, u.range[1], u.range[2])
    v          <- runif(1, v.range[1], v.range[2])
    k.regul    <- runif(1, k.regul.range[1], k.regul.range[2])
    k.pi       <- runif(1, k.pi.range[1], k.pi.range[2])
	deg.rate   <- runif(1, deg.rate.range[1], deg.rate.range[2])
	deg.effect <- runif(1, deg.effect.range[1], deg.effect.range[2])
	ss <- phylosim(n0=n0, G=maxG, u=u, v=v, k.regul=k.regul, k.pi=k.pi, deg.rate=deg.rate, deg.effect=deg.effect, maxTEpop=maxTEpop)
    cat(file=out.file, u, "\t", v, "\t", k.regul, "\t", k.pi, "\t", deg.rate, "\t", deg.effect, "\t", ncol(ss$dist) - 1, "\t", max(ss$dist[-1,-1], na.rm=TRUE), "\n", append=TRUE)
}, mc.cores=mc.cores)

