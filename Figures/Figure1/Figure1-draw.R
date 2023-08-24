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

zlim <- c(-3, 3)
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
dd$logratio <- dd$max.beta / dd$max.alpha

# in case of replicates: take the mean
ddm <- aggregate(. ~ HTgen + eta, dd[,-1], mean)

toplot <- dd$logratio
toplot[toplot < zlim[1]] <- zlim[1]
toplot[toplot > zlim[2]] <- zlim[2]

colplot <- grad[1+round(ncol*((toplot-zlim[1])/diff(zlim)))]

# Main panel

plot(dd$eta, dd$HTgen, xlab=expression(eta), ylab="H", pch=19, col=colplot, cex=2)
