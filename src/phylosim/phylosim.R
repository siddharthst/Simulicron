
phylo.dist <- function(TE.table, TEs, cur.gen, keep.root=FALSE) {
    # Quite unefficient, for many reasons:
    # - distance matrix calculated twice (because of symmetry)
    # - lineages computed every generation (should memorize, since the lineage of a given TE is constant)
    # - common ancestor finding suboptimal (which(X %in% Y) computes much more things that we actually need)
    # - computes lineages up to the initial TE, while keeping track of the MRCA would be possible
    .lineage <- function(i) {
        parent <- TE.table[i,"parent"]
        if (is.na(parent)) 
            return(i)
        c(i, .lineage(parent))
    }
    if (keep.root && !"1" %in% TEs)
        TEs <- c("1", TEs)
    lineages <- lapply(setNames(nm=TEs), .lineage)
    common.parent <- outer(TEs, TEs, function(te1, te2) mapply(te1, te2, FUN=function(tte1, tte2) if (tte1==tte2) NA else lineages[[tte1]][which(lineages[[tte1]] %in% lineages[[tte2]])[1]]))
    
    p.dist <- matrix(vapply(common.parent, FUN=function(pp) {if (is.na(pp)) return(NA); 2*(cur.gen-TE.table[pp, "birth"])}, FUN.VALUE=numeric(1)), ncol=length(TEs), dimnames=list(TEs,TEs))
}

transpo.rates <- function(genome.content, TE.table, u, k.regul, k.pi, cur.gen) {
    a <- TE.table[genome.content, "activity"]
    n <- length(genome.content)
    f.regul <- exp(-k.regul*n)
    f.pi    <- if (k.pi > 0) {
            dd <- phylo.dist(TE.table, genome.content, cur.gen)
            closest <- apply(dd, 2, min, na.rm=TRUE)
            1-exp(-closest*k.pi) 
        } else {
            1
        }
    a*u*f.regul*f.pi
}

phylosim <- function(n0=1, G=1000, u=0.1, v=0, deg.rate=0, deg.effect=1, k.regul=0, k.pi=0, allow.loss=FALSE, maxTEpop=100, maxTEdb=10000) {
    stopifnot(
        n0         >= 1,
        G          >= 1,
        u          >= 0,
        v          >= 0, v          <= 1,
        deg.rate   >= 0, deg.rate   <= 1,
        deg.effect >= 0, deg.effect <= 1,
        k.regul    >= 0,
        k.pi       >= 0
        )
    TE.table <- data.frame(parent=NA, birth=0, activity=1, row.names=as.character(1:n0))
    genome.content <- as.character(1:n0)
    dynam <- list(n=c(n0, rep(NA, G)), a=c(1, rep(NA, G)), u=rep(NA, G+1))
    for (gg in 1:G) {
	if (length(genome.content) == 0)
		break
	if (length(genome.content) > maxTEpop || nrow(TE.table) > maxTEdb)
		break
        # Transposition
        tr.rates <-transpo.rates(genome.content, TE.table, u, k.regul, k.pi, gg)
        transpo  <- rep(genome.content, rpois(length(genome.content), lambda=tr.rates))
        if (length(transpo) > 0) {
        maxID <- as.numeric(rownames(TE.table)[nrow(TE.table)])
        newIDs <- as.character((maxID+1):(maxID+length(transpo)))
	names(transpo) <- newIDs
        TE.table <- rbind(TE.table, data.frame(parent=transpo, birth=gg, activity=TE.table[transpo,"activity"], row.names=newIDs))
        } 
        # new TEs will be added later to the genome (so that they are not deleted nor degraded immediately)

        # Deletion
        dels <- runif(length(genome.content)) < v
        if (all(dels) && !allow.loss) 
            dels[sample(seq_along(dels),1)] <- FALSE
        genome.content <- genome.content[!dels]
        # Degradation
        activities <- TE.table[genome.content,"activity"]
        deg.events <- runif(length(activities)) < deg.rate
        activities[deg.events] <- activities[deg.events] - deg.effect
        activities[activities < 0] <- 0
        TE.table[genome.content,"activity"] <- activities

        # don't forget to add new TEs
        genome.content <- c(genome.content, names(transpo))

        dynam$n[gg+1] <- length(genome.content)
        dynam$a[gg+1] <- mean(activities)
        dynam$u[gg] <- mean(tr.rates)
    }
    list(
        dist = phylo.dist(TE.table, genome.content, G, keep.root=TRUE),
        n    = dynam$n,
        a    = dynam$a,
        u    = dynam$u
    )
}

theor.dyn <- function(u=0.1, v=0, k.regul=0, G=1000, n0=1, ...) {
    ans <- c(n0, rep(NA, G))
    for (gg in 1:G)
        ans[gg+1] <- ans[gg] * (1 + u*exp(-k.regul*ans[gg]) - v)
    ans
}
