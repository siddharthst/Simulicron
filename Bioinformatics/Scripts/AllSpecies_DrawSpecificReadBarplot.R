## Siddharth S. Tomar
## Arnaud Le Rouzic
## Aurelie Hua-Van

# This script is to be used after the python spript GetSpecificReads.py
# It generates a barplot with shared or specific reads for TE-species pairs considered as cross-regulated


# Load the required libs
library(viridis)

#### if no argument provided, use these default
Project="/home/huavan/CrossRegulation"
#sra='SRR14569563'
#sra='SRR5687217'
#sra='SRR25922470'
sra='SRR11846566'
mismatch=3

#### Arguments passed to R scripts
args <- commandArgs(trailingOnly = TRUE)
if (exists(args[1])) Project <- args[1]
if (exists(args[2])) method <- args[2]
if (exists(args[3])) sra <- args[3]
if (exists(args[4])) mismatch <- args[4]

setwd(Project) 

source("Scripts/AllSpecies.Functions.R")
source("Scripts/AllSpecies.Data.R")


col=viridis(4, option='E')

crTE <- read.table(paste0(Project,"/AllSpecies/Table/Species.TE.",sra,".",mismatch,".CrossRegulated.txt"),header=T, sep="\t")
head(crTE)

spe <- read.table(paste0(Project,"/AllSpecies/Stats/AllSpecies.",sra,".",mismatch,".SpeReads.txt"),header=T, sep="\t")
spe<-merge(spe, crTE, by.x=c("Species","TE"), all.x=TRUE)
head(spe)


spe <-spe[order(spe$TE),]
spe$pos<- seq(1,length(spe$TE),1)

#spe$col<-ifelse(spe$FilteredCpNb>1,"black",col[3])
spe$col<-ifelse(spe$FilteredCpNb>1,"black","firebrick")
spe<-spe[spe$NbDmelD>=4000,]

write.table(spe,paste0(Project,"/AllSpecies/Table/Species.TE.",sra,".",mismatch,".CrossRegulated.ReadsNb.Final.txt"),quote=FALSE,col.names=TRUE)

TEspe<-as.data.frame(table(spe$TE))
#spe$NbSpecies<-TEspe$Freq[which(TEspe$Var1==spe$TE)]
dim(TEspe)


pdf('AllSpecies/SupFig/SupFig2_SpecificReads.pdf',width=10, height=14)
par(mar=c(5,10,1,1))
mp<-barplot(spe$NbSpeD+spe$NbDmelD, horiz=T,las=1,cex.names=0.3,col=col[1],cex.axis=0.5,yaxt='n',xlab="Number of mapping small RNAs")
spe$mp=mp
head(spe)
TEpos<-aggregate(mp~TE,data=spe,mean)
dim(TEpos)
TEpos2<- merge(TEpos, TEspe, by.x="TE",by.y="Var1")
TEpos2$col<-ifelse(TEpos2$Freq>1,"black","firebrick")
head(TEpos2)
TEmax<-aggregate(mp~TE,data=spe,max)
print(TEpos$V1)
print(TEmax$V1)
barplot(spe$NbTotD,add=TRUE,col=col[3],horiz=T,xaxt='n',yaxt='n')
barplot(spe$NbSpeD,add=TRUE,col=col[4],horiz=T,xaxt='n',yaxt='n')
axis(2,tick=TRUE,at=TEmax$V1+1/2,labels=FALSE,pos=0)

mtext(TEpos$TE,side=2,line=0.5,at=TEpos$V1,las=2,cex=0.45,col=TEpos2$col)
mtext(spe$Species,side=2,line=3,at=mp,las=2,cex=0.45,col=spe$col)
legend("topright",legend=c("not mapping D. melanogaster copies", "mapping TEs in both species", "mapping D. melanogaster copies only"),pch=15, col=c(col[4], col[3], col[1]))
dev.off()