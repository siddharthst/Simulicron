library(phytools)
library(ape)
library(reshape2)

setwd("/Users/huavan/Dropbox/Siddharth/Simulicron-master/CrossRegulation")
source("Scripts/AllSpecies.Data.R")
source("Scripts/AllSpecies.Functions.R")


# Get the phylogenetic distance to melanogaster
SpeciesDistance <- cophenetic(SpeciesTree)
SpeciesDistance2Dmel <-SpeciesDistance[which(rownames(SpeciesDistance) =="D.melanogaster"),]
SpeciesDistance2Dmel[order(names(SpeciesDistance2Dmel))]

# Get the number of TE copies present in each species
FiltCopyNb <- read.table("AllSpecies/StatsSeq/All.CpNb.filtered.txt")
FiltCopyNb2 <- na.omit(melt(FiltCopyNb)) #useless??

# Get the number of TE families (with more than 1 filtered copies) present in each species
PresentTE<-table(FiltCopyNb2 $V1)
PresentTE<-PresentTE[order(names(PresentTE))]

# Combine
Table4plot<-as.data.frame(cbind(SpeciesDistance2Dmel, PresentTE =PresentTE[names(SpeciesDistance2Dmel)]))


# Correlation
correlation<-cor.test(Table4plot[,1],Table4plot[,2])
pval<-round(correlation$p.value,10)
r<-round(correlation$estimate,4)

# plot
pdf("AllSpecies/SupFig/SupFig.Corr_Phylo_TENb.pdf", width=7,height=7)
colors<-GetColorsPerGroup(taxo)
plot(Table4plot,col=colors[rownames(Table4plot)],pch=19,xlab="Distance to D. melanogaster", ylab="Number of TE families")
mtext(text=paste("r = ",r), side=3, line=-2)
mtext(text=paste("p = ",pval) ,side=3, line=-3)
dev.off()