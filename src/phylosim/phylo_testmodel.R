#!/usr/bin/env Rscript

source("./phylosim.R")

library(parallel)
mc.cores <- detectCores() - 1

todo <- c(
    transpo = TRUE,
    regul   = FALSE,
    degrad  = FALSE
)

# The purpose is to test the properties of the phylogenetic model

replicates <- 100
maxG <- 20
maxTEpop <- 1000

################ Transposition / Deletion

sim.td <- list(
    c(u=0.1, v=0),
    c(u=0.2, v=0),
    c(u=0.3, v=0.1)
)

if (todo["transpo"]) {    
    pdf("fig-test-td.pdf", width=length(sim.td)*5, height=5)
    layout(t(1:length(sim.td)))
    for (td in sim.td) {
        ans <- do.call(rbind, mclapply(1:replicates, function(r) phylosim(n0=1, G=maxG, u=td["u"], v=td["v"], maxTEpop=maxTEpop)$n, mc.cores=mc.cores))
        plot(NULL, ylim=c(0, max(ans, na.rm=TRUE)), xlim=c(0, ncol(ans)), xlab="Time steps", ylab="Copy number", main=paste0("u=", td["u"], ", v=", td["v"]))
        invisible(apply(ans, 1, function(nn) points(nn, pch=1, col="gray")))
        points(colMeans(ans, na.rm=TRUE), pch=1, col="black")
        lines(do.call(theor.dyn, as.list(td)), lwd=3, col="red")
    }
    dev.off()
}

################# Regulation

maxG <- 100
maxTEpop <- 200
replicates <- 100

sim.reg <- list(
    c(u=0.2, v=0.1, k.regul=0, k.pi=0),
    c(u=0.2, v=0.1, k.regul=0.03465736, k.pi=0), # expected: 20 copies at equilibrium
    c(u=0.2, v=0.1, k.regul=0.01386294, k.pi=0)  # expected: 50 copies at equilibrium
)

if (todo["regul"]) {
    pdf("fig-test-reg.pdf", width=length(sim.td)*5, height=5)
    layout(t(1:length(sim.reg)))
    for (td in sim.reg) {
        ans <- do.call(rbind, mclapply(1:replicates, function(r) phylosim(n0=1, G=maxG, u=td["u"], v=td["v"], k.regul=td["k.regul"], k.pi=td["k.pi"], maxTEpop=maxTEpop)$n, mc.cores=mc.cores))
        plot(NULL, ylim=c(0, max(ans, na.rm=TRUE)), xlim=c(0, ncol(ans)), xlab="Time steps", ylab="Copy number", main=paste0("k.regul=", round(td["k.regul"], digits=3)))
        invisible(apply(ans, 1, function(nn) points(nn, pch=1, col="gray")))
        points(colMeans(ans, na.rm=TRUE), pch=1, col="black")
        lines(do.call(theor.dyn, as.list(td)), lwd=3, col="red")
    }
    dev.off()
}


sim.pi <- list(
    c(u=0.2, v=0.1, k.regul=0, k.pi=0.05),
    c(u=0.2, v=0.1, k.regul=0, k.pi=0.04), # expected: 20 copies at equilibrium
    c(u=0.2, v=0.1, k.regul=0, k.pi=0.03)  # expected: 50 copies at equilibrium
)

if (todo["regul"]) {
    pdf("fig-test-pi.pdf", width=length(sim.td)*5, height=5)
    layout(t(1:length(sim.pi)))
    for (td in sim.pi) {
        ans <- do.call(rbind, mclapply(1:replicates, function(r) phylosim(n0=1, G=maxG, u=td["u"], v=td["v"], k.regul=td["k.regul"], k.pi=td["k.pi"], maxTEpop=maxTEpop)$n, mc.cores=mc.cores))
        plot(NULL, ylim=c(0, max(ans, na.rm=TRUE)), xlim=c(0, ncol(ans)), xlab="Time steps", ylab="Copy number", main=paste0("k.pi=", round(td["k.pi"], digits=3)))
        invisible(apply(ans, 1, function(nn) points(nn, pch=1, col="gray")))
        points(colMeans(ans, na.rm=TRUE), pch=1, col="black")
#~         lines(do.call(theor.dyn, as.list(td)), lwd=3, col="red")
    }
    dev.off()
}


############### Degradation



maxG <- 200
maxTEpop <- 200
replicates <- 50

sim.deg <- list(
    c(u=0.2, v=0.1, deg.rate=0, deg.effect=0),
    c(u=0.2, v=0.1, deg.rate=1, deg.effect=0.01), 
    c(u=0.2, v=0.1, deg.rate=0.1, deg.effect=0.1)
)

if (todo["degrad"]) {
    pdf("fig-test-deg.pdf", width=length(sim.deg)*5, height=5)
    layout(t(1:length(sim.deg)))
    for (td in sim.deg) {
        ans <- do.call(rbind, mclapply(1:replicates, function(r) phylosim(n0=1, G=maxG, u=td["u"], v=td["v"], deg.rate=td["deg.rate"], deg.effect=td["deg.effect"], maxTEpop=maxTEpop)$n, mc.cores=mc.cores))
        plot(NULL, ylim=c(0, max(ans, na.rm=TRUE)), xlim=c(0, ncol(ans)), xlab="Time steps", ylab="Copy number", main=paste0("deg.rate=", td["deg.rate"], ", deg.effect=", td["deg.effect"]))
        invisible(apply(ans, 1, function(nn) points(nn, pch=1, col="gray")))
        points(colMeans(ans, na.rm=TRUE), pch=1, col="black")
#~         lines(do.call(theor.dyn, as.list(td)), lwd=3, col="red")
    }
    dev.off()
}
