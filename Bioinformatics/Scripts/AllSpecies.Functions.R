library(ggtree)
library(viridis)
require(tidyverse)
library(ggplot2)
library(phytools)
require(ape)

#################################################################################
#### prepare df for boxplot with filtering or not
#################################################################################
PrepareBoxplot <- function(df,All,Norm){
	test1<- MakeLong(df,All,Norm)
	test2 <-merge(test1,taxo, by='Species', all.x=TRUE)
	test3 <- merge(test2,CpNb1)
	test3 <- na.omit(test3[order(test3$id),])
	test3$group <-factor(test3$group,levels=unique(test3$group))	
	return(test3)
}

#################################################################################
#### Prepare df with both level of mismatch for boxplot
#################################################################################
PrepareFig<- function(New0,New3){
	New0$mismatch <-c(0)
	New3$mismatch <-c(3)
	New <- rbind(New3,New0)
	New<-na.omit(New[order(c(New$id,New$mismatch)),])
	#New<-New[order(c(New$id,New$mismatch)),]
	New$id<- 1:nrow(New)
	return(New)	
}

#################################################################################
##### Generate the dataframe for agregated percent per species
#################################################################################
DoAll<-function(df,All,Norm){
    row_col_pairs0<-MakeLong(df,All,Norm)
    df0<-Mergetab1(row_col_pairs0, CpNb1)
    a0<-aggregate(percent~Species, data= df0,mean)
    b0<-aggregate(percent~Species, data= df0,length)	
    New0 <- Mergetab2(taxo,a0,b0)
    return(New0)
}

#################################################################################
##### Get sra and mismatch from data name
#################################################################################
GetSRAmismatch <- function(x){
	mismatch<-substring(x,nchar(x),nchar(x))
	strain<-substring(x,1,nchar(x)-1)
	return(c(strain, mismatch))
}

GetSRA <- function(x){
	if(x=='TEstats') sra<-"SRR11846566"
	if(x=='w1118') sra<-"SRR14569563"
	if(x=='CanS') sra<-"SRR5687217"
	if(x=='OreR') sra<-"SRR25922470"
	return(sra)
}
#################################################################################
##### Select Cross Regulated data
#################################################################################

GetTableCrossReg <- function(df,All, minpercent,minreads){
	
	# Get the "Dmel.Ref" row
	colRef<-MakeDfRef(df)
	# Transform the df and normalize
	dfN<-MakeLong(df, TRUE,TRUE, minreads)
	# Transform the df and do not normalize
	df<-MakeLong(df, TRUE,FALSE, minreads)
	# Merge Normalize and not Normalized
	df <-merge(dfN, df, all.x=TRUE, by=c("Species","TE"))
	# Merge with colRef (piRNA number in Dmel)
	df <- merge(df,colRef, by='TE', all.x=TRUE)
	# Merge with CpNb
	df <- merge(df,CpNb, all.x=TRUE, by=c("Species","TE"))
	# Merge with taxo
	df <- merge(df,taxo, by='Species', all.x=TRUE)
	# Save a df with only value >=0.7
	dfcr<-df[df$percent>= minpercent,]
	# Remove melanogaster subgroup
	dfcr <- dfcr[dfcr$group != 'melanogaster subgroup',]

	return(dfcr)
}


#################################################################################
##### Make dataframe with Dmel.Ref row
#################################################################################

MakeDfRef<-function(df){
	Ref<-df["Dmel.Ref",]
	colRef<-as.data.frame(t(Ref))
	colnames(colRef)<-c('Dmel.piRNA')
	rownames(colRef) <- gsub("\\.", "-", rownames(colRef))
	colRef$TE <- rownames(colRef)
	rownames(colRef) <- NULL  # Remove the row names from the dataframe
	return(colRef)	
}


#################################################################################
##### Get rid of row Dmel.Ref
#################################################################################
GetRidOfRef<- function(df){
	df <-df[-which(rownames(df)=="Dmel.Ref"),]
	Species<-rownames(df)
	return(df)
	}

#################################################################################
######### Exclude TEs that are detected in less than 15 species
######### Output == Filtered input
#################################################################################
FilterLowSpecies<-function(df,th) df[,which(colSums(!is.na(df))>th)]

#################################################################################
######### Exclude TEs that less than 1/700 in Dmel.Ref
######### Output == Filtered input
#################################################################################
FilterLowDmel<-function(df,th) {
	thresholdDmel=rowSums(df['Dmel.Ref',])
	df<-df[,df['Dmel.Ref',]> thresholdDmel/th]
	return(df)
	}

#######################################################################################
#######    Transform the dataframe (All species + Dmel.Ref) into 2 columns entries
#######	   Output the long dataframe
#######    Option All: all TEs or filtering for 1/700
#######    Option Norm : Normalize by Dmel.Ref or not
#######################################################################################
MakeLong <-function(df,All,Norm,minreads){
	
# Remove TE for which max piRNA in Dmel.Ref is less than 1/700 of the total mapping reads
# All species map poorly with this TEs
if (All==FALSE) {
	df <- df[sapply(df, function(x) max(x, na.rm = T) > sum(df["Dmel.Ref",])/minreads)]
	}

# Remove Dmel.Ref row and save it
RefAll<-df["Dmel.Ref",]
df <-df[-which(rownames(df)=="Dmel.Ref"),]
Species<-rownames(df)

# Normalise by Dmel.Ref nb of reads
if (Norm==TRUE) {
	vecnormAll<-as.vector(RefAll)
	df <- mapply('/', df, vecnormAll,SIMPLIFY = TRUE)
	df <-as.data.frame(df)
	rownames(df)<-Species
	}

#### ALL TEs All Species
non_na_pairsALL <- which(!is.na(df), arr.ind = TRUE)
row_col_pairsALL <- data.frame(
  Species = rownames(df)[non_na_pairsALL[, 1]],
  TE = colnames(df)[non_na_pairsALL[, 2]],
  percent = df[non_na_pairsALL]
)
# Replace '.' by '-' for TEs
row_col_pairsALL $TE <- gsub("\\.", "-", row_col_pairsALL $TE)
if (Norm==FALSE) {
	colnames(row_col_pairsALL)<-c("Species","TE","PiNb")
}
dim(row_col_pairsALL)
#4606
length(unique(row_col_pairsALL$TE))
#175
return(row_col_pairsALL)
}

#################################################################################
############ Merge the long dataframe with TE copy number
#################################################################################
Mergetab1<- function(df,CpNb1){
	dfAll <-merge(df, CpNb1, all=TRUE, by=c("Species","TE"))
	print(head(dfAll))
	return(dfAll)
}

#################################################################################
############ Merge the long dataframe by adding the taxonomy data
#################################################################################
Mergetab2<- function(taxo,a,b){
	taAll <-merge(a,b,by='Species')
	#taAll <- table(df $Species)[match(rev(tips.plotordered), names(table(df $Species))) ] 
	taAll<-  as.data.frame(taAll)
	colnames(taAll)<-c("Species","percent","NbCase")
	New <- merge(taxo,taAll,by=c("Species"))
	New<-New[order(New$id),]
	New$group <-factor(New$group,levels=unique(New$group))
	return(New)
}

#################################################################################
######### Get threshold in Nb of reads that correspond do 1/1000 of Dmel.Ref reads
######### For printing threshold in plots
#################################################################################
GetThreshold <-function(df,th) sum(df["Dmel.Ref",])/th

#################################################################################
######### Get colors for clades based on taxo file
#################################################################################

GetColorsPerGroup<- function(taxo){
	unique_groups <- unique(taxo$group)
	group_colors <- rainbow(length(unique_groups))
	names(group_colors) <- unique_groups
	taxo$col<- group_colors [taxo$group]
	vectcol<- taxo$col
	names(vectcol)<-taxo$Species
	return(vectcol)
}
