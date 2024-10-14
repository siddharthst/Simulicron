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

res.dir <- "../results/Results-fig2"

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

constrain <- function(x, range) {
    ifelse( x > range[2], range[2], ifelse(x < range[1], range[1], x))
}

plot.dots <- function(eta, HTgen, max.copy, grad, zlim, ...) {
    layout(t(1:2), width=c(0.8,0.2))
    par(cex=1, mar=c(4.5,4,2,0.5))
    
    colplot <- grad[1+round(ncol*((constrain(max.copy, zlim)-zlim[1])/diff(zlim)))]
    
    plot(eta, HTgen, xlab=expression("Cross-regulation coefficient ("*eta*")"), ylab=expression("Introduction of TE "*beta*"(H)"), ...)
    rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = col.bg)
    points(eta, HTgen, pch=19, col=colplot, cex=2)
    
    par(mar=c(4.5,2,2,0.5))
    image(x=1, y=seq(zlim[1], zlim[2], length=ncol), z=t(seq(zlim[1], zlim[2], length=ncol)), col=grad, xlab="", ylab="", xaxt="n")
    mtext(side=2, line=-2, col="white", text=expression("Max copy number"))
}

convolve.mat <- function(m, n=7) {
    k <- outer(exp(-seq(-2,2,length.out=n)^2), exp(-seq(-2,2,length.out=n)^2))
    
    stopifnot(nrow(k)%%2 == 1, ncol(k)%%2 == 1)
    
    k <- k/sum(k)
        
    mki <- (nrow(k)-1)/2
    mkj <- (ncol(k)-1)/2
    
    ma <- matrix(NA, ncol=ncol(m), nrow=nrow(m))
    m2 <- matrix(0, ncol=ncol(m)+2*mki, nrow=nrow(m)+2*mkj)
    m2[(mki+1):(nrow(m2)-mki),(mkj+1):(ncol(m2)-mkj)] <- m
    m2[1:mki,1:mkj] <- m[1,1]
    m2[(nrow(m)+mki):nrow(m2), 1:mkj] <- m[nrow(m),1]
    m2[1:mki, (nrow(m)+mki):nrow(m2)] <- m[1,ncol(m)]
    m2[(nrow(m)+mki):nrow(m2),(nrow(m)+mki):nrow(m2)] <- m[nrow(m),ncol(m)]
    for (i in 1:mki) {
        m2[i,(mki+1):(ncol(m2)-mki)] <- m[1,]
        m2[nrow(m2)+1-i,(mki+1):(ncol(m2)-mki)] <- m[nrow(m),]
    }
    for (j in 1:mkj) {
        m2[(mkj+1):(ncol(m2)-mkj),j] <- m[,1]
        m2[(mkj+1):(ncol(m2)-mkj),ncol(m2)+1-j] <- m[,ncol(m)]
    }

    for (i in seq_len(nrow(m))) {
        for (j in seq_len(ncol(m))) {
            ma[i,j] <- sum(k*m2[(i):(i+2*mki),(j):(j+2*mkj)])
        }
    }
    
    ma
}

plot.image <- function(eta, HTgen, max.copy, grad, zlim, ...) {
    par(mar=c(4.5, 4, 2, 0.5))
    
    ee <- sort(unique(eta))
    gg <- sort(unique(HTgen))
    mm <- matrix(NA, nrow=length(gg), ncol=length(ee), dimnames=list(as.character(gg), as.character(ee)))
    
    
    for (i in seq_along(max.copy))
        mm[as.character(HTgen[i]), as.character(eta[i])] <- max.copy[i]
        
    mm.rot <- t(mm)
    
    image(x=ee, y=gg, z=mm.rot, col=grad, zlim=zlim, xlab=expression("Cross-regulation coefficient ("*eta*")"), ylab=expression("Introduction of TE "*beta*"(H)"), ...)
    contour(x=ee, y=gg, z=convolve.mat(mm.rot), zlim=zlim, levels=if (min(mm.rot) < 10) c(1, 10, 20, 30) else c(20, 25, 30), add=TRUE)
}


zlim <- list(
    logratio = c(-0.5, 3), 
    alpha    = c(0, 50),
    beta     = c(0, 50))
    
ncol <- 1024

grad <- list(
    logratio = colorRampPalette(c("yellow", "darkred"))(ncol+1),
    alpha    = colorRampPalette(c("white", col.alpha))(ncol+1),
    beta     = colorRampPalette(c("white", col.beta))(ncol+1))


ff <- list.files(path=res.dir, pattern="*.pickle")

max.cp <- max.copy.number(ff)

dd <- data.frame(
    file = ff,
    HTgen = get.HTgen(ff),
    eta   = get.eta(ff), 
    max.alpha = max.cp[,1],
    max.beta  = max.cp[,2]
)
dd$logratio <- dd$max.alpha / dd$max.beta

# in case of replicates: take the mean
ddm <- aggregate(. ~ HTgen + eta, dd[,-1], mean)

mode <- "image" # Alternative: "image"

pdf("Figure2A.pdf", width=page.width / 2, height=fig.height, pointsize=fontsize)

    if (mode == "dots") {
        plot.dots(ddm$eta, ddm$HTgen, ddm$max.alpha, grad[["alpha"]], zlim[["alpha"]], main=expression("Resident TE ("*alpha*")"))
    } else { # mode = image
        plot.image(ddm$eta, ddm$HTgen, ddm$max.alpha, grad[["alpha"]], zlim[["alpha"]], main=expression("Resident TE ("*alpha*")"))
    }

dev.off()

pdf("Figure2B.pdf", width=page.width / 2, height=fig.height, pointsize=fontsize)
    if (mode == "dots") {
        plot.dots(ddm$eta, ddm$HTgen, ddm$max.beta, grad[["beta"]], zlim[["beta"]], main=expression("Invading TE ("*beta*")"))
    } else { # mode = image
        plot.image(ddm$eta, ddm$HTgen, ddm$max.beta, grad[["beta"]], zlim[["beta"]], main=expression("Invading TE ("*beta*")"))
    }
dev.off()
