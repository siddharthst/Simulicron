library(phytools)
library(ape)

#source("Scripts/AllSpecies.Functions.R")
source("Scripts/AllSpecies.Data.R")

# Get colors for clades based on taxo file

unique_groups <- unique(taxo$group)
group_colors <- rainbow(length(unique_groups))
names(group_colors) <- unique_groups
taxo$col<- group_colors [taxo$group]
vectcol<- taxo$col
names(vectcol)<-taxo$Species

# Get the Newick tree
TEfiles <- list.files(path="AllSpecies/Tree/", pattern=".SRR11846566.3.fasta.tre", full.names=TRUE)

# Arrange the tree and save a pdf if more than 5 tips
for (file in TEfiles){
	TEtree <- ape::read.tree(file)
	TE<-strsplit(basename(file),"\\.")[[1]][1]
	TEtree <-midpoint.root(TEtree)	
	tip.list<-TEtree$tip
	tip_colors <- vectcol[tip.list]
	if (length(tip.list)>5){
		pdf(paste("AllSpecies/SupFig/Tree.",TE,".tre.pdf",sep=''), height=10, width=6.5)
		plot(TEtree,tip.color= tip_colors,edge.width=2,label.offset=0.01,cex = 1,main=TE)
		add.scale.bar()
		dev.off()
		}
	}
