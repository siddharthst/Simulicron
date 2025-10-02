library(parallel)
library(reticulate) 
reticulate::use_python("/usr/bin/python3", required=TRUE)

source("../figures/common.R")

mc.cores <- min(12, detectCores() - 1) # No need for more, parallelizing file I/O is limited

sapply <- function (X, FUN, ..., simplify = TRUE, USE.NAMES = TRUE) 
{
    FUN <- match.fun(FUN)
    answer <- mclapply(X = X, FUN = FUN, ..., mc.cores=mc.cores)
    if (USE.NAMES && is.character(X) && is.null(names(answer))) 
        names(answer) <- X
    if (!isFALSE(simplify) && length(answer)) 
        simplify2array(answer, higher = (simplify == "array"))
    else answer
}

res.dir <- "../results/Results-fig1"

get.HTgen <- function(filenames) {
    as.numeric(sapply(regmatches(filenames, regexec("HT(\\d+)", filenames)), "[", 2))
}

get.eta <- function(filenames) {
    as.numeric(sapply(regmatches(filenames, regexec("eta(\\d\\.\\d\\d)", filenames)), "[", 2))
}

get.sel <- function(filenames) {
    as.numeric(sapply(regmatches(filenames, regexec("sel(\\d\\.\\d+)", filenames)), "[", 2))
}


get.copy.number <- function(filenames) {
    lapply(filenames, function(f) {
        r <- reticulate::py_load_object(file.path(res.dir, f))
        list(alpha=unlist(r$TEfamilyCN[[1]]), beta=unlist(r$TEfamilyCN[[2]]))
    })
}

max.copy.number <- function(filenames) {
    t(sapply(filenames, function(f) {
        r <- reticulate::py_load_object(file.path(res.dir, f))
        c(max(unlist(r$TEfamilyCN[[1]])), max(unlist(r$TEfamilyCN[[2]])))
    }))
}

select.dyn <- function(filenames, HTgen, eta) {
    H <- get.HTgen(filenames)
    e <- get.eta  (filenames)
    ff <- filenames[!is.na(e) & !is.na(H) & HTgen == H & eta == e]
    dyn <- get.copy.number(ff)
    list(
        alpha = t(sapply(dyn, "[[", "alpha")), 
        beta  = t(sapply(dyn, "[[", "beta" ))
    ) 
}

plot.dyn <- function(dyn, ylim=c(0, max(dyn$alpha)), err.bars  = TRUE, xlab="Generations", ylab="Copy number", ...) {
    plot(colMeans(dyn$alpha), type="l", col=col.alpha, ylim=ylim, xlab=xlab, ylab=ylab, ...)
    lines(colMeans(dyn$beta), col=col.beta)
    if (err.bars) {
        xxa <- round(seq(1,ncol(dyn$alpha), length.out=16))
        xxb <- xxa[-1] - diff(xxa[1:2])/2
        arrows(x0=xxa, y0=(colMeans(dyn$alpha)-apply(dyn$alpha, 2, sd))[xxa], y1=(colMeans(dyn$alpha)+apply(dyn$alpha, 2, sd))[xxa], col=adjustcolor(col.alpha, 0.3), length=0)
        arrows(x0=xxb, y0=(colMeans(dyn$beta)-apply(dyn$beta, 2, sd))[xxb], y1=(colMeans(dyn$beta)+apply(dyn$beta, 2, sd))[xxb], col=adjustcolor(col.beta, 0.3), length=0)
    }
}

pch   <- c(1, 0, 19)

ff <- list.files(path=res.dir, pattern="*.pickle")

max.cp <- max.copy.number(ff)

dd <- data.frame(
    file = ff,
    HTgen = get.HTgen(ff),
    eta   = get.eta(ff), 
    sel   = get.sel(ff),
    max.alpha = max.cp[,1],
    max.beta  = max.cp[,2]
)

dd$max.beta[dd$HTgen == max(dd$HTgen)] <- NA # No need to plot the zeros
names(pch) <-  as.character(sort(unique(dd$HTgen)))

ylim.panels <- c(0,40)
xlim.panels <- c(0, 1500)

etas <- c("0", "0.25", "0.5", "0.75", "1")

pdf("Figure1.pdf", width=page.width, height=fig.height, pointsize=fontsize-1)
	layout(matrix(1:(2*length(etas)), ncol=length(etas)))
	par(mar=c(1,1,1,1), oma=c(4,4,1,0), cex=1)
	for (eta in etas) {
		plot.dyn(select.dyn(ff, HTgen=200, eta=eta), xlim=xlim.panels, ylim=ylim.panels, main=substitute("H=200, "*eta*"="*ee, list("ee"=eta)), xaxt="n", yaxt="n", xlab="", ylab="")
		if (eta == etas[1]) {
			axis(2, xpd=NA)
			mtext("Copy number", 2, line=2.5, cex=1.2)
		}
		plot.dyn(select.dyn(ff, HTgen=0, eta=eta), xlim=xlim.panels, ylim=ylim.panels, main=substitute("H=0, "*eta*"="*ee, list("ee"=eta)), xaxt="n", yaxt="n", xlab="", ylab="")
		if (eta == etas[1]) {
			axis(2, xpd=NA)
			mtext("Copy number", 2, line=2.5, cex=1.2)
		}
		axis(1, xpd=NA)
		mtext("Generations", 1, line=2.5, cex=1.2)
	}
	
dev.off()
