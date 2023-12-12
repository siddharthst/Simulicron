library(parallel)
library(reticulate) 

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

res.dir <- "../results/Results-fig3"

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

browser()

pdf("Figure3A.pdf",  width=page.width / 2, height=fig.height, pointsize=fontsize)
    xshft.2A <- setNames(c(-0.02, 0, 0.02), nm=names(pch))
    
    par(mar=c(4.5,4,1,1), cex=1)
    
    dd.eta <- dd[!is.na(dd$eta),]
    
    plot(NULL, xlim=range(dd.eta$eta), ylim=c(0, max(c(dd.eta$max.alpha, dd.eta$max.beta), na.rm=TRUE)), xlab=expression("Cross regulation "*(eta)), ylab="Max copy number")
    rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = col.bg)

    # Means
    mm.alpha <- tapply(dd.eta$max.alpha, INDEX=list(dd.eta$eta, dd.eta$HTgen), FUN=mean)
    mm.beta  <- tapply(dd.eta$max.beta, INDEX=list(dd.eta$eta, dd.eta$HTgen), FUN=mean)
    sd.alpha <- tapply(dd.eta$max.alpha, INDEX=list(dd.eta$eta, dd.eta$HTgen), FUN=sd)
    sd.beta  <- tapply(dd.eta$max.beta, INDEX=list(dd.eta$eta, dd.eta$HTgen), FUN=sd)
    
    for (g in colnames(mm.alpha)) {
        points(as.numeric(rownames(mm.alpha))+xshft.2A[g], mm.alpha[,g], type="p", col=col.alpha, pch=pch[g])
        points(as.numeric(rownames(mm.beta))+xshft.2A[g],  mm.beta[,g],  type="p", col=col.beta, pch=pch[g])
        arrows(x0=as.numeric(rownames(mm.alpha))+xshft.2A[g], y0=mm.alpha[,g]-sd.alpha[,g], y1=mm.alpha[,g]+sd.alpha[,g], code=3, angle=90, length=0.0, col=adjustcolor(col.alpha, 0.3))
        arrows(x0=as.numeric(rownames(mm.beta))+xshft.2A[g],  y0=mm.beta[,g] -sd.beta[,g] , y1=mm.beta[,g] +sd.beta[,g],  code=3, angle=90, length=0.0, col=adjustcolor(col.beta, 0.3))
    }
    
    legend("bottomleft", pch=c(rev(pch)), lty=c(rep(0, 3)), col=c(rep("black", 3)), legend=c("No HT", "HT gen 300", "HT gen 0"), bty="n")
    legend("topright", pch=c(NA, NA), lty=c(0,0), col=c(col.alpha, col.beta), text.col=c(col.alpha, col.beta), legend=c(expression(alpha*" (resident)"), expression(beta*" (invading)")), bty="n")
dev.off()

pdf("Figure3B.pdf",  width=page.width / 2, height=fig.height, pointsize=fontsize)
    xshft.2B <- setNames(c(-0.0001, 0, 0.0001), nm=names(pch))
    
    par(mar=c(4.5,4,1,1), cex=1)
    
    dd.sel <- dd[!is.na(dd$sel),]
    
    xlim <- c(0,0.01) # range(dd.sel$sel)
    ylim <- c(0, 175) # c(0, max(c(dd.sel$max.alpha, dd.sel$max.beta), na.rm=TRUE))
    
    plot(NULL, xlim=xlim, ylim=ylim, xlab=expression("Selection coeffcient (s)"), ylab="Max copy number")
    rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = col.bg)

    # Means
    mm.alpha <- tapply(dd.sel$max.alpha, INDEX=list(dd.sel$sel, dd.sel$HTgen), FUN=mean)
    mm.beta  <- tapply(dd.sel$max.beta, INDEX=list(dd.sel$sel, dd.sel$HTgen), FUN=mean)
    sd.alpha <- tapply(dd.sel$max.alpha, INDEX=list(dd.sel$sel, dd.sel$HTgen), FUN=sd)
    sd.beta  <- tapply(dd.sel$max.beta, INDEX=list(dd.sel$sel, dd.sel$HTgen), FUN=sd)
    
    for (g in colnames(mm.alpha)) {
        points(as.numeric(rownames(mm.alpha))+xshft.2B[g], mm.alpha[,g], type="p", col=col.alpha, pch=pch[g])
        points(as.numeric(rownames(mm.beta))+xshft.2B[g],  mm.beta[,g],  type="p", col=col.beta, pch=pch[g])
        arrows(x0=as.numeric(rownames(mm.alpha))+xshft.2B[g], y0=mm.alpha[,g]-sd.alpha[,g], y1=mm.alpha[,g]+sd.alpha[,g], code=3, angle=90, length=0.0, col=adjustcolor(col.alpha, 0.3))
        arrows(x0=as.numeric(rownames(mm.beta))+xshft.2B[g],  y0=mm.beta[,g] -sd.beta[,g] , y1=mm.beta[,g] +sd.beta[,g],  code=3, angle=90, length=0.0, col=adjustcolor(col.beta, 0.3))
    }
    
    legend("topright", pch=c(rev(pch), NA, NA), lty=c(rep(0, 5)), col=c(rep("black", 3), col.alpha, col.beta), text.col=c(rep("black", 3), col.alpha, col.beta), legend=c("No HT", "HT gen 300", "HT gen 0", expression(alpha*" (resident)"), expression(beta*" (invading)")), bty="n")
dev.off()




