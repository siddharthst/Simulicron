#!/usr/bin/env Rscript

library(ape)

source("./phylosim.R")

library(parallel)
mc.cores <- detectCores() - 1
mut.rate <- 1

sims.file <- "./treecollection.Rdata"

sims <- try(readRDS(sims.file))
if (class(sims) == "try-error") stop("Impossible to load the dataset. Run phylo_treecollection.R first.")


trees <- mclapply(sims, function(sim) lapply(sim, function(ss) {
	dist2tree(ss$dist)
}), mc.cores=mc.cores)

trees.mut <- mclapply(sims, function(sim) lapply(sim, function(ss) {
	dist2tree(ss$dist, mut.rate=mut.rate)
}), mc.cores=mc.cores)

pdf("fig-trees-examples.pdf", width=5*length(sims), height=5)
	layout(t(1:length(sims)))
	for (si in seq_along(sims)) {
		plot(trees[[si]][[1]], main=names(sim.sets)[si])
		axis(1)
	}
dev.off()



#~ coff <- seq(0.5, 0.8, length.out=101)
#~ pdf("fig-trees-pastfuture.pdf", width=5*length(sims), height=5)
#~ 	layout(t(1:length(sims)))
#~ 	for (si in seq_along(sim.sets)) {
#~ 		plot(NULL, xlim=range(coff), ylim=c(-0.3,0.3), xlab="Cutoff", ylab="Slope", main=names(sim.sets)[si])
#~ 		ds <- do.call(rbind, lapply(seq_along(trees[[si]]), function(i) sapply(coff, function(co) { 
#~ 					ns <- branchsuccess(trees[[si]][[i]], cutoff=co)
#~ 					if (nrow(ns) > 1) {
#~ 						ll <- lm(ns[,"branchlength"] ~ ns[,"distance"])
#~ 						coef(ll)[2]
#~ 					} else NA
#~ 				})))
#~ 		if (nrow(ds) > 0) {
#~ 			for (i in 1:nrow(ds)) lines(coff, ds[i,], col="gray")
#~ 			lines(coff, colMeans(ds, na.rm=TRUE), col="black")
#~ 		}
#~ 	}
#~ dev.off()

pdf("fig-trees-branchlengths.pdf", width=5*length(sims), height=10)
	layout(rbind(1:length(sims), (length(sims)+1):(2*length(sims))))
	for (si in seq_along(sim.sets)) {
		branchlengths <- unlist(lapply(trees[[si]], branchlengthdist))
#~ 		plot(density(unlist(branchlengths)), main=names(sim.sets)[si])
#~ 		branchlengths <- branchlengths[branchlengths > 2]
		hist(branchlengths, breaks=100, main=names(sim.sets)[si], freq=FALSE)
		mm <- mean(unlist(branchlengths))
		curve(exp(-x/mm)/mm, add=TRUE, col="red", lwd=2)
		
	}
	for (si in seq_along(sim.sets)) {
		branchlengths.mut <- unlist(lapply(trees.mut[[si]], branchlengthdist))
#~ 		plot(density(unlist(branchlengths)), main=names(sim.sets)[si])
#~ 		branchlengths.mut <- branchlengths.mut[branchlengths.mut > 0]
		hist(branchlengths.mut, breaks=100, main=names(sim.sets)[si], freq=FALSE)
		mm <- mean(unlist(branchlengths.mut))
		curve(exp(-x/mm)/mm, add=TRUE, col="red", lwd=2)
	}
dev.off()
