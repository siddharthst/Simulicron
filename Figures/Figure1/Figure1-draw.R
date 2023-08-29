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

constrain <- function(x, range) {
    ifelse( x > range[2], range[2], ifelse(x < range[1], range[1], x))
}

zlim <- list(
    logratio = c(-0.5, 3), 
    alpha    = c(0, 50),
    beta     = c(0, 50))
    
ncol <- 1024

grad <- colorRampPalette(c("yellow", "darkgreen"))(ncol+1)

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

# Main panel
pdf("Figure1.pdf", width=8, height=5)
    layout(t(1:2), width=c(0.8,0.2))
    
    colplot.logratio <- grad[1+round(ncol*((constrain(ddm$logratio, zlim[["logratio"]])-zlim[["logratio"]][1])/diff(zlim[["logratio"]])))]
    
    plot(ddm$eta, ddm$HTgen, xlab=expression("Cross-regulation coefficient ("*eta*")"), ylab=expression("Introduction of TE "*beta*"(H)"), pch=19, col=colplot.logratio, cex=2, main=expression("Max("*alpha*")/Max("*beta*")"))
    
    image(x=1, y=seq(zlim[["logratio"]][1], zlim[["logratio"]][2], length=length(grad)-1), z=t(seq(zlim[["logratio"]][1], zlim[["logratio"]][2], length=length(grad)-1)), col=grad, xlab="", ylab=expression("Max("*alpha*")/Max("*beta*")"), xaxt="n", yaxt="n")
    axis(2, at=log(c(1, 2, 5, 10)), labels=c(1,2,5,10))
dev.off()

pdf("Figure1-alpha.pdf", width=8, height=5)
    layout(t(1:2), width=c(0.8,0.2))
    
    colplot.alpha <- grad[1+round(ncol*((constrain(ddm$max.alpha, zlim[["alpha"]])-zlim[["alpha"]][1])/diff(zlim[["alpha"]])))]
    
    plot(ddm$eta, ddm$HTgen, xlab=expression("Cross-regulation coefficient ("*eta*")"), ylab=expression("Introduction of TE "*beta*"(H)"), pch=19, col=colplot.alpha, cex=2, main=expression("Resident TE"))
    
    image(x=1, y=seq(zlim[["alpha"]][1], zlim[["alpha"]][2], length=length(grad)-1), z=t(seq(zlim[["alpha"]][1], zlim[["alpha"]][2], length=length(grad)-1)), col=grad, xlab="", ylab=expression("Max("*alpha*")"), xaxt="n")
dev.off()

pdf("Figure1-beta.pdf", width=8, height=5)
    layout(t(1:2), width=c(0.8,0.2))
    
    colplot.beta <- grad[1+round(ncol*((constrain(ddm$max.beta, zlim[["beta"]])-zlim[["beta"]][1])/diff(zlim[["beta"]])))]

    plot(ddm$eta, ddm$HTgen, xlab=expression("Cross-regulation coefficient ("*eta*")"), ylab=expression("Introduction of TE "*beta*"(H)"), pch=19, col=colplot.beta, cex=2, main=expression("Invading TE"))
    
    image(x=1, y=seq(zlim[["beta"]][1], zlim[["beta"]][2], length=length(grad)-1), z=t(seq(zlim[["beta"]][1], zlim[["beta"]][2], length=length(grad)-1)), col=grad, xlab="", ylab=expression("Max("*beta*")"), xaxt="n")
dev.off()
