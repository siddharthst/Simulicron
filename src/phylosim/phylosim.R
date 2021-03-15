

make.symmatrix <- function(M) {
	M[lower.tri(M)] <- t(M)[lower.tri(M)]
	M
}

phylo.dist <- function(TE.table, TEs, cur.gen, keep.root=FALSE, cache.dist=list()) {
    # Still room for performance improvement: 
    # - common ancestor finding suboptimal (match computes more things that we actually need)
    # - computes lineages up to the initial TE, while keeping track of the MRCA would be possible
    .myintersect <- function(x, y) {
		y[match(x, y, 0L)][1]
	}
    .lineage.compute <- function(i) {
        parent <- TE.table[i,"parent"]
        if (is.na(parent)) 
            return(i)
        c(i, .lineage.compute(parent))
    }
    .lineages <- function(vi) {
		cached <- vi %in% names(cache.dist)
		vi.cache <- vi[cached]
		lineages.uncached <- lapply(setNames(nm=vi[!cached]), .lineage.compute)
		cache.dist[vi[!cached]] <<- lineages.uncached
		c(cache.dist[vi.cache], lineages.uncached)
	}
    if (keep.root && !"1" %in% TEs)
        TEs <- c("1", TEs)
    lineages <- .lineages(TEs)
    common.parent <- outer(TEs, TEs, function(te1, te2) mapply(te1, te2, FUN=function(tte1, tte2) 
							if (as.numeric(tte1) >= as.numeric(tte2)) 
								NA 
							else 
								.myintersect(lineages[[tte1]],lineages[[tte2]])[1]
						))
	p.dist <- make.symmatrix(matrix(2*(cur.gen-TE.table[common.parent,"birth"]), ncol=ncol(common.parent)))
    attr(p.dist, "cache.dist") <- cache.dist
	return(p.dist)
}

transpo.rates <- function(genome.content, TE.table, u, k.regul, k.pi, cur.gen, cache.dist=list()) {
    a <- TE.table[genome.content, "activity"]
    n <- length(genome.content)
    f.regul <- exp(-k.regul*n)
    f.pi    <- if (n > 1 && k.pi > 0) {
            dd <- phylo.dist(TE.table, genome.content, cur.gen, cache.dist=cache.dist)
            closest <- apply(dd, 2, min, na.rm=TRUE)
            1-exp(-closest*k.pi) 
        } else {
            1
        }
    ans <- a*u*f.regul*f.pi
    if (exists("dd") && !is.null(dd))
		attr(ans, "cache.dist") <- attr(dd, "cache.dist")
    ans
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
    cache.dist <- list()
    genome.content <- as.character(1:n0)
    dynam <- list(n=c(n0, rep(NA, G)), a=c(1, rep(NA, G)), u=rep(NA, G+1))
    for (gg in 1:G) {
	if (length(genome.content) == 0)
		break
	if (length(genome.content) > maxTEpop || nrow(TE.table) > maxTEdb)
		break
        # Transposition
        tr.rates <-transpo.rates(genome.content, TE.table, u, k.regul, k.pi, gg, cache.dist=cache.dist)
        cache.dist <- attr(tr.rates, "cache.dist")
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
    finaldist <- phylo.dist(TE.table, genome.content, G, keep.root=TRUE, cache.dist=cache.dist)
    diag(finaldist) <- 0
    attr(finaldist, "cache.dist") <- NULL
    colnames(finaldist) <- rownames(finaldist) <- c('1', genome.content)
    list(
        dist = finaldist,
        n    = dynam$n,
        a    = dynam$a,
        u    = dynam$u
    )
}

tree.keepright <- function(tree, cutoff) {
	maxdepth <- max(ape::node.depth.edgelength(tree))
	# trivial = TRUE to get the subtrees with only one node
	# dirty trick to avoid slicing at exactly a node position (which does not behave well) 
	phytools::treeSlice(tree, cutoff*maxdepth-0.0000001, trivial=TRUE)
}

tree.keepleft <- function(tree, cutoff) {
	if (is.null(tree$node.label)) 
		tree$node.label <- paste0("N", 1:tree$Nnode)
	maxdepth <- max(ape::node.depth.edgelength(tree))
	sliced <- tree.keepright(tree, cutoff)
	newtree <- ape::drop.tip(phy=tree, tip=unlist(lapply(sliced, function(tt) tt$tip[-1])))
	
	# Keeping the connection between nodes of tree and newtree is not trivial.
	tiplab <- newtree$tip.label
	newtiplab <- sapply(sliced, function(x) if (ape::Ntip(x) > 1) x$node.label[1] else x$tip.label)
	names(newtiplab) <- sapply(sliced, function(x) tiplab[match(x$tip.label, tiplab, 0L)])
	newtree$tip.label <- setNames(nm=newtiplab[tiplab])
	
	tips <- apply(newtree$edge, 1, function(x) any(1:ape::Ntip(newtree) %in% x))
	newtree$edge.length[tips] <- newtree$edge.length[tips] - (1-cutoff)*maxdepth
	newtree
}

branchsuccess <- function(tree, cutoff) {
	tree$node.label <- paste0("N", 1:tree$Nnode)
	mytree <- tree.keepleft(tree, cutoff)

		# Get the x (distance to the closest neighbor)
	dd <- ape::cophenetic.phylo(mytree) # Back to distances
	x <- apply(as.matrix(dd), 1, function(ll) min(ll[ll>0]))
	
		# Get the y (remaining branch lengths)
	nm.tree <- c(tree$tip.label, tree$node.label)
	raw.branchlengths <- sapply(mytree$tip.label, function(lb) { id <- which(nm.tree == lb); tree$edge.length[tree$edge[,2]==id] })
	y <- raw.branchlengths - sapply(mytree$tip.label, function(lb) { mytree$edge.length[mytree$edge[,2]==which(mytree$tip.label==lb)] } )
	return(cbind(distance=x, branchlength=y))
}

branchlengthdist <- function(tree) {
	ntips <- ape::Ntip(tree)
	tree$edge.length[! tree$edge[,2] %in% 1:ntips]
}

dist2tree <- function(dist.matrix, mut.rate=1e3) {
	# mut rate is per time step
	
	muts      <- make.symmatrix(matrix(rpois(length(dist.matrix), lambda=mut.rate*dist.matrix), ncol=ncol(dist.matrix))/mut.rate)
	colnames(muts) <- rownames(muts) <- colnames(dist.matrix) # Important to keep track of the node labels as "1" is the outgroup.
	njtree    <- ape::nj(muts)
	rootedtree<- ape::root(njtree, outgroup="1", resolve.root=TRUE)
	rootedtree<- ape::drop.tip(rootedtree, "1")
	ultratree <- phytools::force.ultrametric(rootedtree)
}


theor.dyn <- function(u=0.1, v=0, k.regul=0, G=1000, n0=1, ...) {
    ans <- c(n0, rep(NA, G))
    for (gg in 1:G)
        ans[gg+1] <- ans[gg] * (1 + u*exp(-k.regul*ans[gg]) - v)
    ans
}


