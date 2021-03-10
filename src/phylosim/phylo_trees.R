#!/usr/bin/env Rscript

source("./phylosim.R")

library(parallel)
mc.cores <- detectCores() - 1

maxG <- 200
tree.size <- c(40, 60)
tree.depth <- c(50, 150)

maxTEpop <- 200
target.replicates <- 5
max.replicates <- 50
n0 <- 1

sim.sets <- list(
    A=c(u=0.2, v=0.1, k.regul=0.01386294, k.pi=0, deg.rate=0, deg.effect=0),
#~     B=c(u=0.2, v=0.1, k.regul=0, k.pi=0, deg.rate=0.1, deg.effect=0.1),
    C=c(u=0.2, v=0.1, k.regul=0, k.pi=0.03, deg.rate=0, deg.effect=0)
)

sims <- mclapply(sim.sets, function(td) {
    ans <- list()
    repl <- 1
    while(length(ans) < target.replicates && repl <= max.replicates) {
        ss <- phylosim(n0=n0, G=maxG, 
                    u=td["u"], v=td["v"], k.regul=td["k.regul"], k.pi=td["k.pi"], 
                    deg.rate=td["deg.rate"], deg.effect=td["deg.effect"], 
                    maxTEpop=maxTEpop)
        if (ncol(ss$dist) - 1 >= tree.size[1] && 
            ncol(ss$dist) - 1 <= tree.size[2] && 
            max(ss$dist[-1,-1], na.rm=TRUE) >= 2*tree.depth[1] && 
            max(ss$dist[-1,-1], na.rm=TRUE) <= 2*tree.depth[2])
                ans[[length(ans)+1]] <- ss
        repl <- repl + 1
    }
    ans
}, mc.cores=mc.cores)


library(ape)


