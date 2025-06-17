#!/usr/bin/env Rscript
require(ape)
require(phytools)
library(viridis)
require(ggtree)
library(adephylo)

wd=getwd()
setwd(wd)

##### Arguments passed to R scripts
args <- commandArgs(trailingOnly = TRUE)
sra <- args[1]
mismatch <- args[2]
piAnnot <- args[3]
strain <- args[4]

distance=250
print(paste(sra,distance,mismatch,piAnnot,strain,sep='  '))
#Functions
## Transparent colors for TE divergence
## Mark Gardener 2015
## www.dataanalytics.org.uk
## (modified from)
t_col <- function(color,  percent) {
  #      color = color name
  #    percent = % transparency
  #       name = an optional name for the color

## Get RGB values for named color
rgb.val <- col2rgb(color)
## Make new color using input color as base and alpha set by transparency
t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
             max = 255,
             alpha =  max(0,(percent-65)) * 255 / 35
             )

## Save the color
invisible(t.col)
}
## END


# colors for legend of heatmap
ncols <- 101 # number of colors for the image display
# colors for the pi and non pi copies
col.pi <- "darkorange"
col.te <- "darkorchid"           
#col.pi <- magma(4)[3]
#col.te <- magma(4)[2]
          
# Read the file with summary for all TEs
MaxPi.file <- read.table(paste0(wd,'/ResultsDmel.',strain,'/AllTEs/Dmel.', strain, '.', piAnnot,'.',sra,'.',distance,'.',mismatch,'.TE.summary.txt'), header=TRUE)
Max.Pi<-max(na.omit(MaxPi.file$NbReadsKept))
print(Max.Pi)

### list all tree file
listOffiles <- list.files(paste(wd,'/ResultsDmel.',strain,'/Tree',sep=''), pattern=paste(piAnnot,sra,distance,mismatch,'tre',sep='.'), full.names=TRUE)
#print(listOffiles)

#listOffiles<-list(paste(wd,'/ResultsDmel.',strain,'/Tree/',"PROTOP_A.",piAnnot,'.',sra,'.',distance,'.',mismatch,'.tre',sep=''))

for(te in listOffiles){
	
	cat('\n')
	### Get info about TE
	# file fulltpath
	print(te) 
	getTE=unlist(strsplit(basename(te), split=sra))
	TEw=unlist(strsplit(getTE[1], split='[.]'))

	# TE name
	TE=TEw[1]
	
	# Get path
	halfpath=paste(TE,piAnnot,sra,distance,mismatch,sep='.')
	halfpathExt=paste(TE,strain,sra,piAnnot,distance,mismatch,sep='.')
	print(halfpathExt)
	

	# Abort if the tree file is empty
	if (file.info(te)$size == 0){
		print('The file is empty')
		next}
	
	
	##########################################################
	#### MAKE ULTRAMETRIC TREE (FROM .tree file)
	##########################################################
	# open and midpoint root the tree
	treeN<-ape::read.tree(te)
	treeN <- midpoint.root(treeN)
	treeN.size <- mean(adephylo::distRoot(treeN))
	
	# Make tree ultrametric
	mycal <- ape::makeChronosCalib(treeN, age.min=treeN.size)
	treeU <- ape::chronos(treeN, lambda=1, calibration=mycal, quiet=TRUE)	
	
	# Write the ultrametric tree (ape)
	write.tree(treeU,paste0(wd,'/ResultsDmel.',strain,'/Tree/',halfpath,'.ultra.tre'))
	print("Tree done")

	##########################################################
	#### Get Max Nb of Reads for the TE (Kept)
	#### Get Max shared (if 0, next)
	##########################################################
	# From the summary file, Get the number of piMax for the TE
	# Get the total nb of pi that are shared, for the TE	
	MaxPi.TE<-MaxPi.file[MaxPi.file$TE==TE,]
	pi<-MaxPi.TE$NbReadsKept
	print("pi, and ratio")
	print(pi)
	ratio<-round(MaxPi.TE$ReadsInNonPi/MaxPi.TE$ReadsInPi,2)
	print(ratio)
	shared<-MaxPi.TE$Shared
	if (shared==0){	next}
	print(shared)
	
	
	

	####################################
	############## Read BedMat (matrix) which contain a matrix with piRNA counts for each copy of the TE
	####################################
	ReadsDataFrame <- read.table(paste0(wd,'/ResultsDmel.',strain,'/Bowtie/Dmel.', strain, '.',halfpath,'.bed.mat'), header=TRUE)
	
	### Do not analyse if discarded (No euchromatic copies)
	if (file.exists(paste0(wd,'/ResultsDmel.',strain,'/Bowtie/Dmel.', strain, '.',halfpath,'.bed.discarded.txt'))){next}
	
	ReadsDataFrame <-as.data.frame(ReadsDataFrame)
	print(dim(ReadsDataFrame))

	########################
	#### Read the TE tree
	####################################
	tree <- read.tree(paste0(wd,'/ResultsDmel.',strain,'/Tree/',halfpath,'.ultra.tre'))

	##############################################
	#### Read the TE specific summary file on copy
	##############################################
	# Read the TE summary that contains information for each TE copy
	summary <- read.table(paste0(wd,'/ResultsDmel.',strain,'/Bed/',TE,'.',distance,'.', piAnnot,'.summary.txt'))
	### add copy names as rownames (V15)
	rownames(summary) <-summary$V15 # crappy but no better solution since there is no header in the file

	##############################################
	#### Read the coord.txt file containing the coordinates of the fragments, for plotting deletion and truncation
	##############################################
	CoordSeg <- read.table(paste0(wd,'/ResultsDmel.',strain,'/Tree/',halfpath,'.coord.txt'), col.names=c('CopyName','xleft','xright','TE'), header=FALSE,stringsAsFactors = FALSE)


	####################################
	# Get from the tree tips the copyNames that are annotated as in piCluster
	piRNAs <- grepl("piRNA", tree$tip.label)

	# Remove from the tree the copies that are not in the matrix
	tip2remove<-tree$tip[! tree$tip %in% rownames(ReadsDataFrame)]
	tree <- drop.tip(tree, tip2remove)

	# rotate the tree so that a euchromatic copies are at bottom and piCluster on top
	tree <- rotateConstr(tree, c(tree$tip.label[!piRNAs], tree$tip.label[piRNAs]))

	# We need the same order in both figures. Not a better way to get the order from the tree? 
	# Get the copies order as in the tree
	tip.ordered <- tree$tip.label[tree$edge[tree$edge[,2] <= length(tree$tip),2]]
	# Order the piCluster copies as in the tree
	tip.pirev<-rev(tip.ordered[grep("piRNA", tip.ordered)])


	# Continue only if euchromatic copies exists and if more than one pi
	if (length(tip.ordered)==length(tip.pirev)){
		print("==> No euc copies, exiting")
		next
		}
	if (length(tip.pirev)==1){
		print("==> Not enough pi copies, exiting")
		next
		}

	# Count the pi in the matrix                                      
	piRNA.counts <- as.data.frame(ReadsDataFrame[,, drop=FALSE])
	print(length(tip.ordered))
	print(length(tip.pirev))
	print(dim(piRNA.counts))
	
	# Reorder the matrix according to the tree order
	piRNA.counts <- piRNA.counts[tip.ordered, tip.pirev, drop=FALSE]
	colnames(piRNA.counts)<- tip.pirev

	
	### Sub-matrices for getting the max

	# Matrix with only sharing between pi and Nonpi
	piNonpi.counts <-piRNA.counts[! grepl("piRNA", rownames(piRNA.counts)), , drop=FALSE]
	#print(head(piNonpi.counts))
	# Get the max between pi and non pi
	max.piNonpi<-max(piNonpi.counts)
	
	# Get the max between pi 
	max.pipi<-max(piRNA.counts[grepl("piRNA", rownames(piRNA.counts)),])

	#####################################################
	######## Prepare panel with length and divergence to consensus
	###### manage the summary file to get % id and length
	#####################################################
	
	# reduce to keep %len, div, pi ()
	#10=div #7=%start, 13=%len ,15=Cpname

	keptCopies<-summary[summary$V15 %in% tip.ordered,c(10,7,8,13,15) ]

	# Calculate the length of the copies
	keptCopies$len<-keptCopies$V8 - keptCopies$V7 +1
	colnames (keptCopies)<- c("div","start","end","ConsLength","CpName", "len")

	# create a column with colors for the length/%id plot
	# All black...
	keptCopies$col1 <- ifelse(grepl("piRNA" , keptCopies$CpName), "black", "black")

	# Get transparency level based on divergence (it is actually similarity, not div!)
	keptCopies$col <-mapply(t_col,keptCopies$col1, keptCopies$div) 

	# Create a column with % length with a maximum of 1 (not nymore in % of consensus)
	keptCopies$len<- ifelse(keptCopies$len <=keptCopies$ConsLength, keptCopies$len, keptCopies$ConsLength) 

	# For the coordinates
	# Keep the name of the copies to map
	CoordSeg$ybot <-match(CoordSeg$CopyName,tip.ordered)
	# Define the color
	CoordSeg$col <-keptCopies$col[match(CoordSeg$CopyName,keptCopies$CpName)]
	#print(CoordSeg)
	#print(dim(CoordSeg))

	#####################################################
	####### normalisation of the matrix relative to the max piRNA cell
	#####################################################

	#norm.counts <- piRNA.counts / max(piRNA.counts)
	norm.counts <- as.data.frame(piRNA.counts)
	
	# when intersection is between pi and pi, put it as negative, and avoid 0 by adding -1e-6
	norm.counts[grep("piRNA", rownames(norm.counts)),] <- -norm.counts[grep("piRNA", rownames(norm.counts)),]-1e-6

	#####################################################
	#####################################################
	#####################################################
	####### 			Draw the plot
	#####################################################
	#####################################################
	#####################################################

	# Draw the plot
	#pdf(paste(wd,"/Figures/Fig4.",TE,".plot.pdf",sep=''), width=0.5*ncol(ReadsDataFrame), height=0.5*nrow(ReadsDataFrame))
	pdf(paste0(wd,"/Figures4/Fig4.", halfpathExt,".pdf"), width=15, height=10)

	layout(cbind(c(1,1),c(2,2),c(3,3),c(4,5)), widths=c(1, 1 ,0.5,0.3))
	par(mar=c(3, 1, 10, 1), oma=c(0,0,1,0), cex=1,bg="transparent")

	#####################################################
	### Draw the tree (panel 1)
	plot(tree, 
		x.lim           = 1.25*max(node.depth.edgelength(tree)),
		y.lim           = c(0.5, length(tree$tip.label)+0.5), 
		yaxs            = "i",
		tip.color       = ifelse(grepl("piRNA", tree$tip.label), col.pi, col.te), 
		align.tip.label = 0,
		cex			   = 0.8
	)
 
	# TODO : Find a way to put the scale below the tree 
	ape::add.scale.bar()
	
	### Add information
	title(TE, 3, line=4.5,adj=0)

	mtext(paste(pi,"mapped piRNAs",sep=' '), 3, line=3,adj=0)
	mtext(paste("Ratio NonPi-Pi :",ratio,sep=' '), 3, line=2,adj=0)
	#mtext(paste("Max ",max.piNonpi,"pi-Nonpi","Max PiPi",max.pipi,sep=' '), 3, line=1,adj=0)



	#####################################################
	### Draw the heatmap (panel 2)
	# We can use Max.pi but correspond to the max.pi of all TE, so it is better to use pi: Total nb of kept read for the TE
	# The best is to use the max(piRNA.counts) which are the max pi shared between pi and pi or pi and non pi

	#M <-max(piRNA.counts)
	M <-max.piNonpi
	print(M)
	zlim <- c(-M, M) + c(-0.01, 0.01)

	# Denormalize
	#norm.counts <- norm.counts * M

	norm.counts[norm.counts > M] <- M
	norm.counts[norm.counts < -M] <- -M

	image(
    		x      = seq_along(colnames(norm.counts)), 
    		y      = seq_along(rownames(norm.counts)), 
    		z      = t(norm.counts), 
    		axes   = FALSE, 
    		xlab   = "", 
    		ylab   = "", 
    		zlim   = zlim, 
    		#col    = c(colorRampPalette(c(col.pi,gray(0.8)))(floor(ncols/2)), colorRampPalette(c("white", col.te))(ncols-floor(ncols/2))), 
    		col    = c(colorRampPalette(c(col.pi,"gray80"))(floor(ncols/2)), colorRampPalette(c("white", col.te))(ncols-floor(ncols/2))), 
    		breaks = c(seq(min(zlim), -1e-6, length=floor(ncols/2+1)), seq(0, max(zlim), length=ncols-floor(ncols/2)))
    		)

	axis(3, at=seq_along(colnames(norm.counts)), colnames(norm.counts), col.axis=col.pi,las=2, tick=FALSE, cex.axis= min(1,10.0/length(colnames(norm.counts))) )
	mtext("piClusters", col = col.pi, side = 3, line=5)

	#####################################################
	# Draw the barplot for id (div) and length of each copy (panel 3)
	par(mar=c(3,1,10,1))
	lenCons<-keptCopies$ConsLength[1]
	# to replace % per length of consensus, modifiy below

	# 1- white bars bordered by black
	barplot(height= keptCopies$ConsLength, width=.5, yaxs = "i",horiz=TRUE, space = 1, ylim=c(0.25,length(keptCopies$len)+0.25),col='white',border='gray',cex.axis=0.5)

	#barplot(height= rep(1, length(keptCopies$len)), width=.5, yaxs = "i",horiz=TRUE, space = 1, ylim=c(0.25,length(keptCopies$len)+0.25),col='white',border='gray')

	# 2- colored bars from 0 to end of copy
	#barplot(height= keptCopies$len+keptCopies$start, width=.5, yaxs = "i",border=NA, col    = keptCopies$col,horiz=TRUE, space = 1, ylim=c(0.25,length(keptCopies$len)+0.25), add=TRUE, xaxt='n')

	# 3- white bars boredered by black from 0 to start of copy
	#barplot(height= keptCopies$start, width=.5, yaxs = "i",border="NA", col    = 'white',horiz=TRUE, space = 1, ylim=c(0.25,length(keptCopies$len)+0.25), add=TRUE, xaxt='n')

	# 4- solution with gap and rectangle
	#(xyxy)
	nbcopie<-length(keptCopies$len)
	rect(CoordSeg$xleft,CoordSeg$ybot-0.5,CoordSeg$xright,CoordSeg$ybot,col    = CoordSeg$col,border=NA)

	text("Alignment and similarity to consensus", y=length(keptCopies$len)+4,x=lenCons/2,xpd=NA)
	text(c("60%","100%"), x=c(0, lenCons),y=length(keptCopies$len)+2.6,cex=0.7,xpd=NA)

	# Similarity bar: TO IMPROVE
	rect(seq(0, lenCons, lenCons/100),length(keptCopies$len)+1.5,seq(lenCons/100, lenCons+lenCons/100, lenCons/100),length(keptCopies$len)+2, 	xpd=NA,col=colorRampPalette(c("gray80","black"))(101),border=NA)

	#####################################################
	##Draw the legend for the piRNA counts (panel 4 and 5)
	yy.pi <- seq(0, -zlim[1]*max(piRNA.counts), length.out=floor(ncols/2)+1)
	yy.te <- seq(0, zlim[2]*max(piRNA.counts), length.out=floor(ncols/2)+1)
	yy.pi <- yy.te <- seq(0, M, length.out=floor(ncols/2)+1)

	par(mar=c(3,1,10,6))
	##Draw the legend top
	image(
    		y= yy.pi,
    		z=t(as.matrix(yy.pi)), 
    	xlab="", ylab="", axes=FALSE,
    	col=colorRampPalette(c(gray(0.8),col.pi))(floor(ncols/2))
	)

	axis(4)
	mtext("pi-Cluster hits", 4, line=2.5)

	# Draw the legend bottom   
	image(
    		y= yy.te,
    		z=t(as.matrix(yy.te)), 
    		xlab="", ylab="", axes=FALSE,
    		col=colorRampPalette(c("white",col.te))(floor(ncols/2))
	)
	axis(4)
	mtext("TE hits", 4, line=2.5)
   
   
	dev.off()
}
