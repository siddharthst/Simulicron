library(parallel)
library(reticulate) 

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




res.dir <- "./Results"

get.HTgen <- function(filenames) {
    as.numeric(sapply(regmatches(filenames, regexec("HT(\\d+)", filenames)), "[", 2))
}

get.eta <- function(filenames) {
    as.numeric(sapply(regmatches(filenames, regexec("eta(\\d\\.\\d\\d)", filenames)), "[", 2))
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
    ff <- filenames[HTgen == H & eta == e]
    dyn <- get.copy.number(ff)
    list(
        alpha = t(sapply(dyn, "[[", "alpha")), 
        beta  = t(sapply(dyn, "[[", "beta" ))
    ) 
}

plot.dyn <- function(dyn, ylim=c(0, max(dyn$alpha)), err.bars  = TRUE, ...) {
    plot(colMeans(dyn$alpha), type="l", col=col["alpha"], ylim=ylim, xlab="Generations", ylab="Copy number", ...)
    lines(colMeans(dyn$beta), col=col["beta"])
    if {err.bars) {
        xx <- round(seq(1,ncol(dyn$alpha), length.out=16))
        arrows(x0=xx, y0=(colMeans(dyn$alpha)-apply(dyn$alpha, 2, sd))[xx], y1=(colMeans(dyn$alpha)+apply(dyn$alpha, 2, sd))[xx], col=adjustcolor(col["alpha"], 0.3), length=0)
         arrows(x0=xx, y0=(colMeans(dyn$beta)-apply(dyn$beta, 2, sd))[xx], y1=(colMeans(dyn$beta)+apply(dyn$beta, 2, sd))[xx], col=adjustcolor(col["beta"], 0.3), length=0)
	}
}

col   <- c(alpha="blue", beta="red")
pch   <- c(1, 0, 19)
xshft <- c(-0.02, 0, 0.02)

ff <- list.files(path=res.dir, pattern="*.pickle")

max.cp <- max.copy.number(ff)

dd <- data.frame(
    file = ff,
    HTgen = get.HTgen(ff),
    eta   = get.eta(ff), 
    max.alpha = max.cp[,1],
    max.beta  = max.cp[,2]
)

dd$max.beta[dd$HTgen == max(dd$HTgen)] <- NA # No need to plot the zeros
names(pch) <- names(xshft) <- as.character(sort(unique(dd$HTgen)))

pdf("Figure2.pdf", width=5, height=4)
    par(mar=c(5,4,1,1))
    
    plot(NULL, xlim=range(dd$eta), ylim=c(0, max(c(dd$max.alpha, dd$max.beta), na.rm=TRUE)), xlab=expression("Cross regulation "*(eta)), ylab="Max copy number")
    
    # All data points
    # points(dd$eta, dd$max.alpha, col=adjustcolor(col["alpha"], 0.5), pch=pch[as.character(dd$HTgen)])
    # points(dd$eta, dd$max.beta,  col=adjustcolor(col["beta"], 0.5),  pch=pch[as.character(dd$HTgen)])

    # Means
    mm.alpha <- tapply(dd$max.alpha, INDEX=list(dd$eta, dd$HTgen), FUN=mean)
    mm.beta  <- tapply(dd$max.beta, INDEX=list(dd$eta, dd$HTgen), FUN=mean)
    sd.alpha <- tapply(dd$max.alpha, INDEX=list(dd$eta, dd$HTgen), FUN=sd)
    sd.beta  <- tapply(dd$max.beta, INDEX=list(dd$eta, dd$HTgen), FUN=sd)
    
    for (g in colnames(mm.alpha)) {
        points(as.numeric(rownames(mm.alpha))+xshft[g], mm.alpha[,g], type="p", col=col["alpha"], pch=pch[g])
        points(as.numeric(rownames(mm.beta))+xshft[g],  mm.beta[,g],  type="p", col=col["beta"], pch=pch[g])
        arrows(x0=as.numeric(rownames(mm.alpha))+xshft[g], y0=mm.alpha[,g]-sd.alpha[,g], y1=mm.alpha[,g]+sd.alpha[,g], code=3, angle=90, length=0.0, col=adjustcolor(col["alpha"], 0.3))
        arrows(x0=as.numeric(rownames(mm.beta))+xshft[g],  y0=mm.beta[,g] -sd.beta[,g] , y1=mm.beta[,g] +sd.beta[,g],  code=3, angle=90, length=0.0, col=adjustcolor(col["beta"], 0.3))
    }
    
    legend("bottomleft", pch=c(rev(pch), NA, NA), lty=c(rep(0, 3), rep(1,2)), col=c(rep("darkgray", 3), col), legend=c("No HT", "HT gen 300", "HT gen 0", expression(alpha*" (resident)"), expression(beta*" (invading)")))
dev.off()

pdf("Figure1bis.pdf", width=12, height=8)
    layout(rbind(1:3,4:6))
    
    ylim <- c(0,40)
    
    plot.dyn(select.dyn(ff, HTgen=300, eta=0),   ylim=ylim, main=paste0("H=300, eta=0"))
    plot.dyn(select.dyn(ff, HTgen=300, eta=0.5), ylim=ylim, main=paste0("H=300, eta=0.5"))
    plot.dyn(select.dyn(ff, HTgen=300, eta=1),   ylim=ylim, main=paste0("H=300, eta=1"))
    
    plot.dyn(select.dyn(ff, HTgen=0, eta=0),     ylim=ylim, main=paste0("H=0, eta=0"))
    plot.dyn(select.dyn(ff, HTgen=0, eta=0.5),   ylim=ylim, main=paste0("H=0, eta=0.5"))
    plot.dyn(select.dyn(ff, HTgen=0, eta=1),     ylim=ylim, main=paste0("H=0, eta=1"))

dev.off()
