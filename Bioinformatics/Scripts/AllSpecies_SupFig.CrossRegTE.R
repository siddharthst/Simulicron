## Siddharth S. Tomar
## Arnaud Le Rouzic
## Aurelie Hua-Van

# Load the required libs
suppressPackageStartupMessages({
require(tidyverse)
require(ape)
require(ggtree)
})
#setwd("~/CrossRegulation")
#setwd("~/Dropbox/Siddharth/Simulicron-master/CrossRegulation")
#setwd("/mnt/Poles/Genomes/huavan/CrossRegulation")
setwd("/home/huavan/CrossRegulation")

source("Scripts/AllSpecies.Functions.R")
source("Scripts/AllSpecies.Data.R")

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

options(repr.plot.width = 10, repr.plot.height = 15)

# minimum species to be kept
minSpe<-20
minRead<-700
minpercent<-0.7

################################
# For selecting CrossRegulated TEs
################################
# Table with piRNA counts
TEstats <- read.table(paste0("AllSpecies/Table/",sra,'.',mismatch,".final_result_readcount.txt"), header = TRUE,  sep="\t")
TEstatsF<- GetTableCrossReg(TEstats, TRUE, minpercent, minRead)

# Save CrossRegulated table
write.table(TEstatsF,paste("AllSpecies/Table/Species.TE.",sra,'.',mismatch,'.',"CrossRegulated.txt",sep=""),quote=FALSE,col.names=TRUE, row.names=FALSE,sep='\t')

SpeciesPerTE<-MakeLong(TEstats,TRUE,TRUE,0)
pdf("AllSpecies/SupFig/SupFig3.pdf", width=4,height=7)
par(mfrow=c(3,1))
head(SpeciesPerTE)
dim(SpeciesPerTE)
barplot(sort(table(SpeciesPerTE$TE)))
barplot(sort(table(SpeciesPerTE$Species)))
dev.off()
#### END
