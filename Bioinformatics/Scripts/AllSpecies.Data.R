library(ape)
#setwd("/Users/huavan/Dropbox/Siddharth/Simulicron-master/CrossRegulation/")

iso1<-"SRR11846566"
w1118<-"SRR14569563"
CanS<-"SRR5687217"
OreR<-"SRR25922470"
	
taxo <- read.table("Data/TaxoDroso.txt",fill=TRUE,sep='\t')
colnames(taxo)<-c('Species','group')
taxo$id<-1:nrow(taxo)

CpNb0<- read.table("AllSpecies/StatsSeq/AllSpecies.CpNb.txt",fill=TRUE)
colnames(CpNb0)<-c("Species",'TE',"InitialCpNb","NA","FilteredCpNb")
CpNb0=CpNb0[,c(1,2,3,5)]

# Put 0 for NA
CpNb <-CpNb0
CpNb[is.na(CpNb)]<-0

# Put NA for 0
CpNb1 <-CpNb0
CpNb1[CpNb1==0]<-NA
#CpNb1==CpNb0

TypeTE<-read.table("Data/RepBaseDmel.Classif.txt")



TEstats0 <- read.table(paste0("AllSpecies/Table/",iso1,'.',0,".final_result_readcount.txt"), header = TRUE,  sep="\t")
dim(TEstats0)#105 142

TEstats3 <- read.table(paste0("AllSpecies/Table/",iso1,'.',3,".final_result_readcount.txt"), header = TRUE,  sep="\t")
dim(TEstats0)#105 142

TEw11180 <- read.table(paste0("AllSpecies/Table/",w1118,'.',0,".final_result_readcount.txt"), header = TRUE,  sep="\t")
dim(TEstats0)#105 142

TEw11183 <- read.table(paste0("AllSpecies/Table/",w1118,'.',3,".final_result_readcount.txt"), header = TRUE,  sep="\t")
dim(TEstats0)#105 142

TECanS0 <- read.table(paste0("AllSpecies/Table/",CanS,'.',0,".final_result_readcount.txt"), header = TRUE,  sep="\t")
dim(TEstats0)#105 142

TECanS3 <- read.table(paste0("AllSpecies/Table/", CanS,'.',3,".final_result_readcount.txt"), header = TRUE,  sep="\t")
dim(TEstats0)#105 142

TEOreR0 <- read.table(paste0("AllSpecies/Table/", OreR,'.',0,".final_result_readcount.txt"), header = TRUE,  sep="\t")
dim(TEstats0)#105 142

TEOreR3 <- read.table(paste0("AllSpecies/Table/", OreR,'.',3,".final_result_readcount.txt"), header = TRUE,  sep="\t")
dim(TEstats0)#105 142

######## Tree
tree <- ape::read.tree("Data/BigGene_.tre")
# Keep only those tips which are present in the result dataframe

SpeciesList<-c(row.names(TEstats3))
SpeciesList <-SpeciesList[SpeciesList!="Dmel.Ref"]

tree <- keep.tip(tree, SpeciesList)
# Define outgroup
outgroup <- c("C.costata", "L.varia")
# Root the tree
tree <- root(tree, outgroup)
# Rotate the tree
SpeciesTree <- rotateConstr(tree, c("D.melanogaster", tree$tip[!tree$tip %in% c("D.melanogaster",outgroup)], outgroup))

# Get the plotting order
tips.plotordered <- SpeciesTree $tip[SpeciesTree $edge[SpeciesTree $edge[,2] <= length(SpeciesTree $tip),2]]
