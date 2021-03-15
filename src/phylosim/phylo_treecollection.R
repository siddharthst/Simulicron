#!/usr/bin/env Rscript

# Generates a collection of trees generated from various phyloegentic parameters

source("./phylosim.R")

library(parallel)
mc.cores <- detectCores() - 1

sims.file <- "./treecollection.Rdata"

maxG <- 201
tree.size <- c(40, 60)
tree.depth <- c(100, 200)

maxTEpop <- 200
target.replicates <- 50
max.replicates <- 1000
n0 <- 1

sim.sets <- list(
    'Regul'    =c(u=0.5, v=0.1, k.regul=0.03, k.pi=0, deg.rate=0, deg.effect=0),
    'Regul+deg'=c(u=0.5, v=0.1, k.regul=0.03, k.pi=0, deg.rate=0.02, deg.effect=0.1),
    'Pi'       =c(u=0.5, v=0.1, k.regul=0, k.pi=0.007, deg.rate=0, deg.effect=0),
    'Pi+deg'   =c(u=0.5, v=0.1, k.regul=0, k.pi=0.009, deg.rate=0.02, deg.effect=0.1)
)

sims <- mclapply(sim.sets, function(td) {
    ans <- list()
    repl <- 1
    while(length(ans) < target.replicates && repl <= max.replicates) {
        ss <- phylosim(n0=n0, G=maxG, 
                    u=td["u"], v=td["v"], k.regul=td["k.regul"], k.pi=td["k.pi"], 
                    deg.rate=td["deg.rate"], deg.effect=td["deg.effect"], 
                    maxTEpop=maxTEpop)
#~         print(c(ncol(ss$dist)-1, max(ss$dist[-1,-1], na.rm=TRUE)/2))
        if (ncol(ss$dist) - 1 >= tree.size[1] && 
            ncol(ss$dist) - 1 <= tree.size[2] && 
            max(ss$dist[-1,-1], na.rm=TRUE) >= 2*tree.depth[1] && 
            max(ss$dist[-1,-1], na.rm=TRUE) <= 2*tree.depth[2])
                ans[[length(ans)+1]] <- ss
        repl <- repl + 1
    }
    ans
}, mc.cores=mc.cores)

for (simn in names(sims)) {
	if (length(sims[[simn]]) < target.replicates)
		warning("Simulation ", simn, ", ", length(sims[[simn]]), "/", target.replicates, " generated.")
}

saveRDS(sims, file=sims.file)
