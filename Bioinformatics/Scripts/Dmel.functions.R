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
		NamesOldPi <- names(ExternalPiCopiesDepth[ExternalPiCopiesDepth > Youth*2])	
		NamesOtherPi <- names(ExternalPiCopiesDepth[! names(ExternalPiCopiesDepth) %in% NamesOldPi]  )	
		return(list('NamesOldPi'= NamesOldPi, 'NamesPiDescendantMRCA' =NamesPiDescendantMRCA, 'NamesEucDescendantMRCA' = NamesEucDescendantMRCA, 'NamesOtherPi'= NamesOtherPi))		
	}
