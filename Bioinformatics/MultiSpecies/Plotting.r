# Load the required libs
suppressPackageStartupMessages({
require(tidyverse)
require(ape)
require(ggtree)
})

options(repr.plot.width = 10, repr.plot.height = 15)

# Define minimum number of piRNA required for a TE
minPiRNA <- 200

# Define TEs to keep - which should be shared between all species
NAfraction <- 0.4

# Read TE stats
TEstats <- read.table("./Results/Stats.tsv", header = TRUE, row.names = 1)

# Remove TEs which do not have enough piRNA
TEstats <- TEstats[sapply(TEstats, function(x) max(x, na.rm = T) > minPiRNA)]

# Also remove R12 as we cannot find any copy in DM 
TEstats <- subset(TEstats, select = -c(R12_DM_NonLTR_Retrotransposon))                          

# Read tree
tree <- ape::read.tree("./BigGene_.tre")

# Keep only those tips which are present in the result dataframe
tree <- keep.tip(tree, c(row.names(TEstats)))

# Define outgroup
outgroup <- c("C.costata", "L.varia")

# Root the tree
tree <- root(tree, outgroup)

# Rotate the tree
tree <- rotateConstr(tree, c("D.melanogaster", tree$tip[!tree$tip %in% c("D.melanogaster",outgroup)], outgroup))

# Get the plotting order
tips.plotordered <- tree$tip[tree$edge[tree$edge[,2] <= length(tree$tip),2]]

# Rearrange dataframe
TEstats <- TEstats[match(rev(tips.plotordered), row.names(TEstats)), ] 

# Remove NA columns
TEstats <- TEstats[, which(colMeans(!is.na(TEstats)) > NAfraction)]

# Convert to matrix
TEstats <- as.matrix(TEstats)

# Normalise by DM reads
TEstats <- t( t(TEstats) / TEstats["D.melanogaster",])

# Make >1 to 1
TEstats[TEstats > 1] <- 1.0

# Make 0 to -1
TEstats[TEstats == 0] <- - 1.0                          

# Fix column names
colnames(TEstats) <- gsub("_.*", "", colnames(TEstats))


# Prepare plotting matrix
# PlottingMatrix <- TEstats %>% replace(is.na(.), -1)
                          
breaks <- c(-1.0, seq(from = 0.0, to = 1, by = 0.3), 1)

p1 <- ggtree(tree, ladderize = FALSE) + geom_tiplab(size = 4, align=TRUE, linesize=.5)
gheatmap(p1, TEstats, offset=0.2, width=5, colnames_angle=90, color="black", colnames_offset_y=-3) + 
scale_fill_viridis_b(option = "E", direction = -1, na.value = 'white', name="Normalized piRNA hits") + ggtree::vexpand(.1, -1)

p2 <- ggtree(tree, ladderize = FALSE) + geom_tiplab(size = 4, align=TRUE, linesize=.5)
gheatmap(p1, TEstats, offset=0.2, width=5, colnames_angle=90, color="white", colnames_offset_y=-3) + 
scale_fill_viridis_b(option = "A", direction = -1, na.value = 'gray90', name="Normalized piRNA hits", breaks = breaks) + ggtree::vexpand(.1, -1)

ggsave("./Results/Plot.pdf", width = 50, height = 50, units = "cm", limitsize = FALSE)
ggsave("./Results/Plot_binned.pdf", width = 50, height = 50, units = "cm", limitsize = FALSE)
# scale_fill_gradient2(high="blue", mid="navy", low="white", midpoint=0, limits=c(0,1))
# scale_fill_continuous(type = "gradient", name="Normalized piRNA hits", na.value = 'white')
# ggsave("./Results/Plot.pdf", width = 50, height = 50, units = "cm", limitsize = FALSE)


