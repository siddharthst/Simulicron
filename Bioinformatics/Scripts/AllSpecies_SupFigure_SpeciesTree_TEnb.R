require(tidyverse)
library(viridis)
require(ape)
require(ggtree)
library(ggplot2)
library(phytools)

setwd('~/Dropbox/Siddharth/Simulicron-master/CrossRegulation')

source("Scripts/AllSpecies.Functions.R")
source("Scripts/AllSpecies.Data.R")

### Script for making figures (OK)
### AllSpecies/SupFig/Nb.TE.Species.pdf
### AllSpecies/SupFig/Tree.Nb.TE.Species.pdf
### AllSpecies/SupFig/Nb.TE.families.pdf

# CpNb1: All Species-TE pairs, with initial and filtered copy number

# Class TE by class and family, and compute nb of species
# CpNb1: copy number with NA when no copies
# TabT1: for each TE, the number of species carrying related copies + TE type
tabT <-as.data.frame(table(CpNb1$TE))
FreqT<-setNames(tabT$Freq,tabT$Var1)
tabT1<-merge(tabT,TypeTE,by.x='Var1',by.y='V1')

tabS <-as.data.frame(table(CpNb1$Species))
tabS <- tabS[match(rev(tips.plotordered), tabS$Var1), ] 
FreqS<-setNames(tabS$Freq,tabS$Var1)

# Number of reads for each TE pair = TEstats3
# Calculate the percent reads relative to Dmel, using MakeLong(, AllTE,Normalized by Dmel)
TabCR <-MakeLong(TEstats3,TRUE,TRUE)
# Select TE-Species pairs for which the number of mapping reads is >= 70 % of what is observed for melanogaster
TabCR <-na.omit(TabCR[TabCR $percent>=0.7,])
dim(TabCR)

# How much TE fam with more than 70% of reads in the different species
FreqCR<-table(TabCR$Species)
FreqCR<-FreqCR[rev(tips.plotordered)]

df<-TEstats3
All=FALSE
Norm <-TRUE
df<-FilterLowDmel(df,1000)
df <- FilterLowSpecies(df,15)
test3Select<-PrepareBoxplot(df,All,Norm)

##### Barplot Number of each Type of TE shared by species
pdf("AllSpecies/SupFig/Nb.TE.Species.pdf", width=6, height=9)
par(mfrow=c(1,1))
hist(FreqT, breaks=seq(0,95+5,5),xlab="Number of sharing species",ylab="Number of TE families", main="Distribution of TE families among species")
dev.off()

##### Barplot Number of each Type of TE shared by species
tabT1$V3 <- factor(tabT1$V3, levels = unique(tabT1$V3))
colnames(tabT1)<-c("TE","Freq","V1","Subclass")
custom_breaks <- seq(0, 100, 5)
plot<-ggplot(tabT1, aes(x = Freq, fill = Subclass)) +
  geom_histogram(binwidth = 1, breaks = custom_breaks, position = "stack") +
  geom_histogram(binwidth = 1, breaks = custom_breaks, position = "stack", color = "white", linewidth = 0.2) +
  scale_fill_viridis_d(option = "E") +  # Use viridis color scale
  labs(x = "Number of shared species", y = "Number of TEs", title = "Distribution of TEs among species") +
  theme_minimal()+
  theme(panel.grid = element_blank(),  # Remove grid lines
        axis.line = element_line(color = "black"),  # Add axis lines
        axis.ticks = element_line(color = "black"),  # Add tick marks on the axes
        panel.background = element_blank(),
        legend.position = "none") +  # Adjust legend position
  scale_x_continuous(breaks = seq(0, 100, by = 5)) +  # Specify x-axis breaks
  scale_y_continuous(breaks = seq(0, 40, by = 5))  + # 
  coord_flip()  # Flip coordinates to make it horizontal

pdf("AllSpecies/SupFig/Nb.TE.Species.2.pdf", width = 6, height = 9)
print(plot + theme(legend.position = c(0.8, 0.8),  # Position legend inside the plot area
            legend.justification = "center",  # Justify legend to the top
            legend.box = "horizontal",  # Arrange legend items horizontally
            legend.margin = margin(0)))  # Remove any margin around the legend box

dev.off()

zap<-unique(colors[grepl('Z.' , names(colors))])
#"#00CCFF"

######## Plot the tree + the number of TE per species All
Freq<-setNames(tabS$Freq,tabS$Var1)
cw<-reorder(SpeciesTree,"cladewise")
plotTree(cw,tips=NULL,plot=FALSE)
colors<-GetColorsPerGroup(taxo)
zap<-colors[9]

obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
xlim<-obj$x.lim
ylim<-obj$y.lim
offset<-0.2*diff(xlim)
offset<-0.02
xlim[2]<-xlim[2]+offset

pdf("AllSpecies/SupFig/SupFig1B.Nb.TE.families.Tree.pdf", width=6, height=9)
layout(matrix(c(1,2),1,2),widths=c(0.5,0.5))
par(fg="black",mar=c(2,0,0,0))
plot.new()
plot.window(xlim=xlim,ylim=ylim)
text(rep(max(obj$xx[1:Ntip(cw)])+offset,Ntip(cw)),obj$yy[1:Ntip(cw)],
    labels=gsub("_"," ",cw$tip.label),font=3, col=colors[cw$tip.label] ,pos=4,cex=0.6,offset=-0.05)
for(i in 1:Ntip(cw)) lines(c(obj$xx[i],max(obj$xx[1:Ntip(cw)])+offset),
    rep(obj$yy[i],2),lty="dotted")

col=viridis(2, option='E')
par(fg="transparent",mar=c(2,0,0,0))
plotTree(cw,fsize=0.6,lwd=2,add=TRUE,tip.labels=FALSE,xlim=xlim,ylim=ylim, mar=c(2,0,0,0))
par(fg="black",mar=c(2,0,0,2))
barplot(rev(FreqS),horiz=TRUE,names.arg=FALSE,col=col[1],xaxt='n')
barplot(rev(FreqCR),horiz=TRUE,add=TRUE,names.arg=FALSE,col=col[2],xaxt='n')
axis(1,at=seq(0,160,50),labels=seq(0,160,50),line=-1)
dev.off()


cw<-reorder(SpeciesTree,"cladewise")
plotTree(cw,tips=NULL,plot=FALSE)

obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
xlim<-obj$x.lim
ylim<-obj$y.lim
offset<-0.2*diff(xlim)
offset<-0.02
xlim[2]<-xlim[2]+offset

vectcol <- GetColorsPerGroup(taxo)
vectgroup<- unique(taxo$group)
groupcol<- unique(vectcol)

pdf("AllSpecies/SupFig/SupFig1D.Legend.rainbow.Tree.pdf", width=5, height=1.2)
par(mar=c(0,0,1,0))
plot.new()
title(main="Clade color legend")
legend("topleft",legend=vectgroup[1:5],text.col= groupcol[1:5] , bty="n",cex=0.8)
legend("top",legend=vectgroup[6:10],text.col= groupcol[6:10] , bty="n",cex=0.8)
legend("topright",legend=vectgroup[11:15],text.col= groupcol[11:15] , bty="n",cex=0.8)
dev.off()



pdf("AllSpecies/SupFig/SupFig1C.Nb.TE.families.rainbow.Tree.pdf", width=6, height=9)
layout(matrix(c(1,2),1,2),widths=c(0.5,0.5))
par(fg="black",mar=c(2,0,0,0))
plot.new()
plot.window(xlim=xlim,ylim=ylim)
text(rep(max(obj$xx[1:Ntip(cw)])+offset,Ntip(cw)),obj$yy[1:Ntip(cw)],
    labels=gsub("_"," ",cw$tip.label),font=3,pos=4,cex=0.6,offset=-0.05, col= vectcol[obj$yy[1:Ntip(cw)]])
for(i in 1:Ntip(cw)) lines(c(obj$xx[i],max(obj$xx[1:Ntip(cw)])+offset),
    rep(obj$yy[i],2),lty="dotted")

col=viridis(2, option='E')
par(fg="transparent",mar=c(2,0,0,0))
plotTree(cw,fsize=0.6,lwd=2,add=TRUE,tip.labels=FALSE,xlim=xlim,ylim=ylim, mar=c(2,0,0,0))
par(fg="black",mar=c(2,0,0,2))
plot(test3Select$percent,1:nrow(test3Select),las=1,cex=0.5,xlim=c(0,1),cex.axis=0.6,yaxt='n',bty='n', ylab='Fraction of mapping reads',col= vectcol[test3Select$Species])
legend("topright",legend=c("mismatch = 3","Selected TEs"),bty="n",cex=0.5)
dev.off()

################# CORRELATION PLOT
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
pdf("AllSpecies/SupFig/SupFig1A.Corr_Phylo_TENb.pdf", width=4,height=7)
par(mar=c(5,2,1,2))
colors<-GetColorsPerGroup(taxo)
plot(Table4plot,col=colors[rownames(Table4plot)],pch=19,xlab="Distance to D. melanogaster", ylab="Number of TE families")
mtext(text=paste("r = ",r), side=3, line=-2)
mtext(text=paste("p = ",pval) ,side=3, line=-3)
dev.off()

