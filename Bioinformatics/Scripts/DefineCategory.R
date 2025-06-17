#!/usr/bin/env Rscript
require(ape)
require(phytools)
require(ggtree)
require(gridExtra)
library(ggplot2)
library("ggpubr")
library(adephylo)
library(castor)

### Siddharth S. Tomar
### Arnaud Le Rouzic
### Aurélie Hua-Van
# Server version

#####################################################################
### This script analyses TE trees (if any, from kept or discarded TE families), various summary files, and the piRNA matrix
### It uses (for each TE family) a tree file (if any), a matrix file with piRNA count
### It needs the Classification file (in Data)
### It needs a TE summary file (data for copy number and for piRNA, for each TE in one line)
### It needs a summary file (data for each copy of each TE in one line)
### It determines automatically if there is old pi sequences and young euchromatic copies
### It determines automatically the regulatory status, the presence of bursts
### It generates a table with some statistics for each TE families (AllTEs/...summary.txt.2.txt)
#####################################################################

# Comment this for debug/verbose output
# print <- function(...) {invisible()}

wd=getwd()
setwd(wd)

#source("Dmelfunctions")
#####################################################################
############## FUNCTIONS #############################
####### To put into a file Dmelfunctions.R ##############
#####################################################################

CheckParent<-function(x, nod,trees){
  y<- tree$edge[which(x==tree$edge[,2])]
  z<- getDescendants(trees,x)
  w<-z[z>length(tree$tip)] # other internal nodes
  v<-z[z %in% tip.nodes]
  u<-tip.nodes[tip.nodes %in% v ]
  f<- length(u[names(u %in% NbEucYoung)])	    
  if (y==mrca & length(u)>2)	 return(length(u))
}

# function to determine if several waves when NbEucYoung >=4 (method2) 
# function to generate two list of names to compare (find mrca and their depth) 
f1 <- function(x) setNames(as.data.frame(t(combn(x, 2))), c("x", "y"))

# function to determine if several waves when NbEucYoung >=4 (method1) 
case <-function(x){
	y=x-2 
	return((y^2-y)/2)
	}

# Function to get the Names of old pi copies (external to mrca of young copies and > 0.01)
GetNamesOldPi<- function(x,y){
		### Get all copies descendant of mrca of young euk copies
 		AllDescendantMRCA<-getDescendants(tree,mrcaYoung)		
		### Get the names
		NamesAllDescendantMRCA <-names(tip.nodes[tip.nodes %in% AllDescendantMRCA])
		### Get Name of all Pi copies descendant of mrca of young euk copies
		NamesPiDescendantMRCA <-NamesAllDescendantMRCA[grepl('piRNA',NamesAllDescendantMRCA)]
		### Get Name of all Pi copies descendant of mrca of young euk copies
		NamesEucDescendantMRCA <-NamesAllDescendantMRCA[! grepl('piRNA',NamesAllDescendantMRCA)]
		### Get Name of Externalcopies
		NamesExternalCopies <-tree$tip[! tree$tip %in% NamesAllDescendantMRCA]
		### Get Depth of Externalcopies
		ExternalCopiesDepth <- tip.depths[names(tip.depths) %in% NamesExternalCopies]		
		### Get Depth of PiExternalcopies
		ExternalPiCopiesDepth <- ExternalCopiesDepth[grepl('piRNA', names(ExternalCopiesDepth))]		
		### Get Names of old PiExternalcopies
		NamesOldPi <- names(ExternalPiCopiesDepth[ExternalPiCopiesDepth > Youth2])	
		NamesOtherPi <- names(ExternalPiCopiesDepth[! names(ExternalPiCopiesDepth) %in% NamesOldPi]  )	
		return(list('NamesOldPi'= NamesOldPi, 'NamesPiDescendantMRCA' =NamesPiDescendantMRCA, 'NamesEucDescendantMRCA' = NamesEucDescendantMRCA, 'NamesOtherPi'= NamesOtherPi))		
	}

#####################################################################
#####################################################################

##### Arguments passed to R scripts
args <- commandArgs(trailingOnly = TRUE)
sra <- args[1]
mismatch <- args[2]
piAnnot <- args[3]
strain <- args[4]
distance=250

##############################
#### Parameters
Youth <- 0.004 # Max Branch length for considering young insertion
Treg <- 0.33 # regulation factor (fraction YE relative to old Pi reads)
Texp <- 0.5 # regulation factor (fraction YO relative to old Pi reads)
TF <- 0.3 # for calculating the proportion of young euchromatic copies (Not used)
min4burst <- 6
min4Highburst <- 15
Youth2 <- Youth * 1 # Max Branch length for considering young insertion in piRNA cluster


### number of piRNA in the different dataset 
NbReadSRA<-c('SRR11846566'=15679203, 'SRR5687217'=3646851, 'SRR25922470'=11326065)
# Total number of decollapsed reads. 
minReads<-NbReadSRA[sra]/15000

### Read the Classification table 
Classif<-read.table(paste(wd,"/Data/RepBaseDmel.Classif.txt",sep=''))
colnames(Classif)<-c("TE","Superfamily","Class")

### Read the Copy summary table (all the copies, and their status)
Copyfilename<-paste(wd,'/ResultsDmel.',strain,'/AllTEs/Dmel.',strain,'.',piAnnot,'.',sra,'.',distance,'.',mismatch,'.Cp.summary.txt',sep='')
CopyTab <- read.table(Copyfilename, header=TRUE, fill=TRUE)
print(dim(CopyTab))

### Read the summary table
Summaryfilename<-paste(wd,'/ResultsDmel.',strain,'/AllTEs/Dmel.',strain,'.',piAnnot,'.',sra,'.',distance,'.',mismatch,'.TE.summary.txt',sep='')
tab <- read.table(Summaryfilename, header=TRUE, fill=TRUE)

# List of TE to analyse
listTE<- tab$TE

# Create empty vector to gather all new calculations
vec<-c()
excluded<-c()

####################################################################################
##################### TE by TE #####################################################
####################################################################################

for(TE in sort(listTE)){
	cat("\n")
	print(TE)

	# combined parameters for saving file
	halfpath=paste(TE,piAnnot,sra,distance,mismatch,sep='.')
	halfpathExt=paste(TE,strain,sra,piAnnot,distance,mismatch,sep='.')
	
	#################################################
	#### Get the type of element (Class, superfamily)
	#### From the classification file
	#################################################	
	TEclass <- Classif$Class[Classif$TE==TE]
	TEfam <- Classif$Superfamily[Classif$TE==TE]
	
	#######################################
	#### Processing of the table for Status
	#### From the summary file
	#######################################	
	tabTE<-tab[which(tab$TE==TE),]
	# Select columns (2:TotCpNb,3:FiltCpNb,4:KeptCpNb,5:PiCpNb,6:EucCpNb,9:TotalNbReads,10:Freqminus,12:ReadsInNonPi,13:ReadsInPi,18:Status)
	Col2Keep<-c(2,3,4,5,6,9,10,12,13,18)
	
	tabTEv <- c(TE, TEclass, TEfam, as.vector(tabTE[Col2Keep]))
	names(tabTEv)<-c('TE', 'Class', 'Superfamily', colnames(tabTE[Col2Keep]))
	# Change the name of the TE status not to nr confused with the regulation status
	names(tabTEv)[length(tabTEv)]<-"StatusTE"
	# 13 columns
	
	#### Detect if TE has been kept of discarded (for the tree)
	if (tabTE$Status=='kept') pattern<-paste(TE,piAnnot,sra,distance,mismatch,'tre',sep='.')
	if (tabTE$Status=='discarded') pattern<-paste(TE,piAnnot,sra,distance,mismatch,'Discarded.tre',sep='.')
	
	################################################
	#### Processing of the Copy table for Divergence
	################################################
	TECopyTab<-CopyTab[which(CopyTab$TE==TE),]
	print(dim(TECopyTab))
	
	
	########################################
	# Defaults values for variables
	Status <- 'NA'	
	burst <- 'NA'	
	NamesOldPi<-c()
	naseries<-rep('NA',10)
	

	RatioNonPiPi<-tabTEv$ReadsInNonPi/tabTEv$ReadsInPi
	
	MaxTPTP<-'NA'
	MaxTETP<-'NA'
	MaxOPOP <-'NA'
	MaxYEOP<-'NA'
	MaxYEYP<-'NA'
	MaxYPOP<-'NA'
	MaxVec<-c(RatioNonPiPi, MaxTPTP, MaxTETP, MaxOPOP , MaxYEOP, MaxYEYP,'NA','NA')
	names(MaxVec)<-c("ratioNonPiPi","MaxTPTP", "MaxOPOP" , "MaxYEOP", "MaxYEYP","ratioOPTP","ratioOPYP")
	#### Detect the cause of the exclusion (if kept, will be 'NA', to be filled later)
	
	### Excluded <- EucCpNb==0 -> NoBurst  # The tree (if any) will have only pi Copies --> abort
	### No pi copies <- PiCpNb + ExcPi ==0 # The tree (if any) will have only euc Copies --> check if burst or not
	### Not enough piRNAs <- TotalNbReads/PiCpNb<1000 --> check if burst or not
	### Not enough piRNAs <- TotalNbReads<10000 --> check if burst or not
	
	# No euc copy
	if(tabTE$EucCpNb==0) {
		Status<-"Excluded"
		burst<-"NoBurst"
		excluded<-c(excluded, tabTEv, MaxVec,naseries)
		next
	# There is Euc copies
	# But no pi copy
	} else if(tabTE$PiCpNb==0 & tabTE$ExcPi==0) {
		Status<-"No pi copy"

	# There is euc and pi copies but not expressed
	} 
	

	############################
	##### PROCESSING OF THE TREE
	############################

	# We check if the tree file exists
	if (file.exists(paste0(wd,'/ResultsDmel.',strain,'/Tree/',pattern)) ==FALSE){
	print('The file does not exist')
	excluded <-c(excluded, tabTEv, MaxVec ,naseries)
	next}
		
	# We do not want empty files
	if (file.info(paste0(wd,'/ResultsDmel.',strain,'/Tree/',pattern))$size ==0){
	print('The file is empty')
	excluded <-c(excluded, tabTEv, MaxVec ,naseries)
	next}
	
	# Otherwise, Read the tree
	tree <- read.tree(paste0(wd,'/ResultsDmel.',strain,'/Tree/',pattern))

	# We do not want tree with less than 4 branches
	if (length(tree$tip)<4){
	print('The tree has less than 4 branches')
	#NbEukTot <- length(tree$tip.label[! grepl('piRNA',tree$tip.label)])	
	#if(NbEukTot==0) Status<-'Excluded'
	excluded <-c(excluded, tabTEv, MaxVec ,naseries)
	next}
	
	# Reroot the tree (AHV: no idea why it is done twice )
	tree <- midpoint.root(tree)
	#tree <- phytools::midpoint.root(unroot(tree))

	#######################################################################################
	#######################################################################################
	#######################################################################################
	#######################################################################################
	###############
	##### Stat calculation from the tree
	###############

	#### 0. tips info 
	tip.names<- tree$tip
	Ntips<-length(tip.names)
	tip.nodes<-sapply(tip.names,function(x,y) which(y==x),y=tree$tip.label)
	tip.depths<-setNames(tree$edge.length[sapply(tip.nodes, function(x,y) which(y==x),y=tree$edge[,2])],names(tip.nodes))
    #print(tip.nodes)
    
	#### 1. nodes info
	node.depths<-get_all_node_depths(tree) ### == length tips - 1
	Nnodes<-length(node.depths)
	a<-as.data.frame(get_all_node_depths(tree, as_edge_count=FALSE))
	nodes.nodes<-as.numeric(rownames(a))+ Ntips
	a$nodes.nodes<-nodes.nodes
	colnames(a)<- c('a','nodes.nodes')
	names(node.depths)<-nodes.nodes

	
	# Relative distance between the root and the next node
	# If large, the next node is 
	Others<-a$a[which(a$a !=max(a$a))]
	rati<-1-(max(Others)/max(a$a))	
		
	#### Get the names and depths of all euchromatic tips 
	EucTot<-tip.depths[! grepl("piRNA",names(tip.depths))]
	NEucTot<- length(names(EucTot))
	
	#### Get the names of all piCopies 	
	PiTot <-tip.depths[grepl("piRNA",names(tip.depths))]
	PiTot.names <- names(PiTot)
	NPiTot<- length(PiTot.names)
	print(paste('NPiTot', NPiTot))
	
	############################
	##### LENGTH EDGES OF YOUNG EUCHROMATIC COPIES
	##### FILTER YOUNG EUCHROMATIC COPIES BASED ON DIVERGENCE TO CONSENSUS
	############################

	#### 1. Select the young copies based on terminal edge length (depth) (named vector)
	Young.depths <- tip.depths[tip.depths<= Youth]

	#### 2. Select euchromatic copies (names) with low divergence from the copy table and not piRNA
	TECopyYoung <- TECopyTab$CopyName[TECopyTab$div >90 & ! grepl("piRNA",TECopyTab$CopyName)]
	
	#### 3. Select the euchromatic young copies (named vector)
	EucYoung1.depths <- Young.depths[! grepl("piRNA",names(Young.depths))]
	EucYoung.depths <- Young.depths[names(Young.depths) %in% TECopyYoung]
	
	# Save the name
	YoungEuc.names <-names(EucYoung.depths)
    #print("EucYoung.names")
    #print(YoungEuc.names)
    
    # Calculate CpNb for each category
	# NbYoung <- length(Young.depths)
	NbEucYoung <- length(EucYoung.depths) # <= Youth and >90%
	NbEucOld   <- length(names(EucTot))-length(names(EucYoung1.depths)) # Tot - Young (<= Youth)
	NbEucOther <-length(names(EucYoung1.depths)) - NbEucYoung # Young but <= 90 % (may indicate recent wave of divergent element)
	
	#### 4. Proportion of young euchromatic copies among all euchromatic copies
	PEucYoung <- length(EucYoung.depths)/length(names(EucTot))

		
	###########################################
	########## MRCA OF YOUNG EUCHROMATIC COPIES
	########## NODE AND DEPTH
	###########################################
	
	### 1. Get the node (number) of mrca of the young euchromatic tips
	# need to have at least one copy. If NbEucYoung ==0 and NEucTot==1, take this one, and 
	# consider all piCopies as external
	if (NbEucYoung>=2){
		mrcaYoung<-get_mrca_of_set(tree, names(EucYoung.depths))
		print(mrcaYoung)
		# Get the depth of the node
		mrcaYoungDepth<-node.depths[mrcaYoung-Ntips] #2 if it does not exist
		#print(mrcaYoungDepth)
	# if NbEucYoung ==1, then get_mrca_of_set returns the tip, we take its length
	} else if (NbEucYoung==1){
		mrcaYoung<-get_mrca_of_set(tree, names(EucYoung.depths))
		# Get the length of the terminal branch
		mrcaYoungDepth<-tip.depths[mrcaYoung] #2 if it does not exist		
	# if there is only one euc copie, with a long branch
	} else if (NbEucYoung==0 & NEucTot==1){
		# We consider it as a young copy
		mrcaYoung<-tip.nodes[names(EucTot)]
		# Get the length of the terminal branch
		mrcaYoungDepth<-tip.depths[names(EucTot)] #2 if it does not exist
	} else { # NbEucYoung==0 & NEucTot > 1
		mrcaYoung=0
		mrcaYoungDepth<- 0
	}
	
	# Correct the mrcaYoungDepth (when >1)
	mrcaYoungDepth <-ifelse(length(mrcaYoungDepth)>1,0, mrcaYoungDepth)
	### 2. Relative distance to root: The larger the most recent
	propmrca<- mrcaYoungDepth/max(a$a)
	print(paste("mrcaYoung: ", mrcaYoung))
	
	
	###########################################
	########## MRCA OF ALL PI COPIES
	###########################################
	if (NPiTot>=2){
		mrcaPi<-get_mrca_of_set(tree, PiTot.names)
		# Get the depth of the node
		mrcaPiDepth<-node.depths[mrcaPi-Ntips] #2	
	} else if (NPiTot ==1){
		mrcaPi <-get_mrca_of_set(tree, PiTot.names)
		# Get the length of the terminal branch
		mrcaPiDepth <-tip.depths[mrcaPi] #2 if it does not exist
	} else { # No pi Copy. Status must be "No pi copy" or "Not enought piRNAs" or 'NA' (if lots of piRNAs)
		mrcaPiDepth <- 0		
	}
	mrcaPiDepth <-ifelse(length(mrcaPiDepth)>1,0, mrcaPiDepth)
			
	

	##############################################
	### DETERMINE IF THERE IS NEW AND OLD PI COPIES 
	##############################################

    # Defaults values (empty vectors for names)
	NamesOldPi<-c()
	NamesPiDescendantMRCA<-c()
	NamesEucDescendantMRCA<-c()
	NamesOtherPi<-c()	
	
	NamesExternalCopies<-c()
	ExternalCopiesDepth<-c()
	ExternalPiCopiesDepth<-c()

	
	# 1- mrcaYoung is the root: No external copy
	# In this case, there is no external copy
	if (mrcaYoung - Ntips == 1 & Status =='NA'){
		Status<-'No old pi copy'
		res <- GetNamesOldPi(tree, mrcaYoung)
	    NamesOldPi<-res[[1]]
	    NamesPiDescendantMRCA <-res[[2]]
	    NamesEucDescendantMRCA <-res[[3]]
	    NamesOtherPi <-res[[4]]


	# 2- if mrcaYoung is among the tip (meaning there is only one young copy)
	# mrcaYoung==0 
	# In this case all piCopies are external
	} else if (mrcaYoung <= Ntips) {
		ExternalCopiesDepth <- tip.depths		
		ExternalPiCopiesDepth <- ExternalCopiesDepth[grepl('piRNA', names(ExternalCopiesDepth))]	
		NamesOldPi <- names(ExternalPiCopiesDepth[ExternalPiCopiesDepth > Youth2])	### 2- NEW COPIES (WITHIN THE BURST)	
		NamesOtherPi <- names(ExternalPiCopiesDepth[ExternalPiCopiesDepth <= Youth2])	

	# 3- if mrcaYoung is not in the tip and not the root
	### The good stuff
	} else if (mrcaYoung - Ntips != 1){	
	
		res <- GetNamesOldPi(tree, mrcaYoung)
		NamesOldPi<-res[[1]]
		NamesPiDescendantMRCA <-res[[2]]
		NamesEucDescendantMRCA <-res[[3]]
		NamesOtherPi <-res[[4]]

	}

	##############################################
	### Determine if there is recent burst or not
	##############################################

	if (NbEucYoung >= min4Highburst) {
			burst<-"RecentHighBurst"
	} else if (NbEucYoung >= min4burst) {
			burst<-"RecentBurst"
	} else if (NbEucYoung < min4burst & PEucYoung >= TF) {
			burst<-"NoBurst"
	} else if (NbEucYoung < min4burst & PEucYoung < TF) {
			burst<-"OldBurst"
	}
		
	##################################################
	######## Fill regulation status	(1)
	##################################################
	if (length(YoungEuc.names) ==0){   		
     	Status="No young Euc" # Cat E
 	} else if (length(PiTot.names)==0){
		Status='No pi copy' # Cat D
	} else if (length(NamesOldPi) <= 0){
		Status='No old pi copy' # Cat C
	} else if (Status=='NA' & length(NamesOldPi)>0){
		Status<-"Old pi copy" # Cat A, B
	} else if (length(NamesOldPi)>0){
		Status<-"Old pi copy"			
	} 
		
	##################################################
	##### PROCESSING OF THE MATRIX WHEN "OLD PI COPIES" 
	##################################################
	# For old copies determined by the phylogeny, how much piRNA?
	
	ReadsDataFrame <- as.data.frame(read.table(paste0(wd,'/ResultsDmel.',strain,'/Bowtie/Dmel.', strain, '.',halfpath,'.bed.mat'), header=TRUE))
	
	#### Selection of copies that are present in the matrix
	#### Some copies not mapped at all by piRNA are not in the matrix. This concerns TE families with very few piRNA of course!
	#### The tree may have been drawed while the TE family was initially discarded
	
	# Select old pi copies that are kept in the matrix (the tree may have been drawed while the TE family was initially discarded)
	PiTot.names <-PiTot.names[PiTot.names %in% colnames(ReadsDataFrame)]
	# Select Euc copies that are kept in the matrix (the tree may have been drawed while the TE family was initially discarded)
	EucTot.names <- names(EucTot[names(EucTot) %in% rownames(ReadsDataFrame)])

	FinalNamesOldPi<-NamesOldPi
	# Same as NamesOldPi
	OldPi.names <-FinalNamesOldPi[FinalNamesOldPi %in% colnames(ReadsDataFrame)]
	YoungPi.names <-NamesPiDescendantMRCA[NamesPiDescendantMRCA %in% colnames(ReadsDataFrame)]
	YoungEuc.names <-YoungEuc.names[YoungEuc.names %in% rownames(ReadsDataFrame)]
	print(paste(OldPi.names, PiTot.names, EucTot.names ,sep=' '))
	
	
	MaxYPOP   <- max(as.data.frame(ReadsDataFrame[YoungPi.names, OldPi.names, drop=FALSE]))	
	MaxOPOP   <- max(as.data.frame(ReadsDataFrame[OldPi.names, OldPi.names, drop=FALSE]))
	MaxYEOP   <- max(as.data.frame(ReadsDataFrame[YoungEuc.names, OldPi.names, drop=FALSE]))
	MaxYEYP   <- max(as.data.frame(ReadsDataFrame[YoungEuc.names, YoungPi.names, drop=FALSE]))

	### Calculate Max Total shared
	MaxTPTP   <- max(as.data.frame(ReadsDataFrame[PiTot.names, PiTot.names, drop=FALSE]))
	MaxTETP  <- max(as.data.frame(ReadsDataFrame[EucTot.names, PiTot.names, drop=FALSE]))

	MaxVec<-c(RatioNonPiPi,MaxTPTP, MaxOPOP , MaxYEOP, MaxYEYP, round(MaxOPOP/MaxTPTP,2),round(MaxYEOP/MaxYEYP,2))
	#Avoid -Inf in the final table
	MaxVec[MaxVec==-Inf]<-0
	MaxVec[MaxVec==Inf]<-0
	

	
	#Fill regulation status	(2)
	##### Determination of Status when there is "Old pi copies"
	if (length(OldPi.names) ==0 & length(FinalNamesOldPi)>0 & Status=='Old pi copy') {
		# There is old pi from the tree but not in the matrix
		Status<-"Old pi not expressed"
    # At least one old pi copy
	} else if (length(FinalNamesOldPi) ==0 & Status=='Old pi copy') {
		print("FinalNamesOldPi == O")
		Status<-"No old pi copy"
	} else if (nrow(ReadsDataFrame)>0 & ncol(ReadsDataFrame)>0 & Status=='Old pi copy' ){
		# enough pi
   		if (MaxTPTP < minReads){Status <-"Old pi not expressed"# Old pi sequences cross-mapped by more than 0.3 of the young pi cross-mapped piRNA
     		} else if (MaxYEOP >= MaxOPOP *Treg & MaxOPOP  >= Texp * MaxTPTP){Status <-"Old pi regulating"
			} else if (MaxOPOP  < Texp * MaxTPTP ) {Status <- "Old pi not expressed"
			} else if (MaxOPOP  >= Texp * MaxTPTP & MaxYEOP < MaxOPOP *Treg ) {Status <- "Old pi expressed but not regulating"
			} else {Status <- "Old pi but unsure"} # Should not happen
    }
    
	    
	vec<-c(vec,tabTEv, MaxVec,max(a$a),NbEucYoung , NbEucOther, NbEucOld, length(NamesOldPi), length(NamesOtherPi) ,length(NamesPiDescendantMRCA), mrcaYoungDepth,Status, burst)	
	

} # end of the for loop

names(MaxVec)<- c("RatioNonPiPi","MaxTPTP", "MaxOPOP" , "MaxYEOP","MaxYEYP","ratioOPTP","ratioOPYP")
Names3<-c('RootDepth','NbEucYoung','NbEucOther', 'NbEucOld','NbPiOld', 'NbPiOther','NbPiYoung', 'mrcaDepth','Status','Burst')
FinalNames<-c(names(tabTEv),names(MaxVec),Names3)

# Write the final summary file with the status and burst for each family
Summary2<-as.data.frame(matrix(vec, ncol=length(FinalNames), byrow=TRUE))
Summary2 <-apply(Summary2,2,as.character)
colnames(Summary2)<- FinalNames
write.table(Summary2, paste(Summaryfilename,'.2.txt',sep=''),row.names=FALSE,col.names=TRUE,quote=FALSE, sep='\t')