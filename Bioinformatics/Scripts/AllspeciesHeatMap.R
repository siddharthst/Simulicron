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


#sra='SRR14569563'
#sra='SRR5687217'
#sra='SRR25922470'
sra='SRR11846566'
mismatch=3

#### Arguments passed to R scripts
args <- commandArgs(trailingOnly = TRUE)
Project <- args[1]
method <- args[2]
sra <- args[3]
mismatch <- args[4]
setwd(Project) # give an error

options(repr.plot.width = 10, repr.plot.height = 15)


# parameters
minSpe<-20
minRead<-700
minpercent<-0.7
print(sra)

#################################################################################################
### piRNA mapping on the best copy: Create the dataframe
####################################################################################################
# List all your files
files <- list.files(path=paste0("AllSpecies/Bowtie.",sra,".",mismatch), pattern=".*\\.readcount.txt", full.names=TRUE)

# Read and merge files using Reduce and merge
result <- Reduce(function(x, y) merge(x, y, by = "rowname", all = TRUE, sort = FALSE),
                 lapply(files, function(file) {
                   data <- read.table(file, header=FALSE, col.names=c("rowname", "contig", "value"), stringsAsFactors=FALSE)
                   setNames(data[, c("rowname", "value")], c("rowname", basename(file)))
                 }))

# Transpose the result
colnames(result) <- sub("\\.readcount\\.txt", "", colnames(result))
rownames(result) <- result$rowname
result <- result[,-1 ]


# Write the transposed result to a file with row names
write.table(t(result), file=paste0("AllSpecies/Table/",sra,'.',mismatch,".final_result_readcount.txt"), quote=FALSE, row.names=TRUE, col.names=TRUE, sep="\t")

####################################################################################################
## Read TE stats
####################################################################################################
# no need to reopen: = t(result)
TEstats0 <- read.table(paste0("AllSpecies/Table/",sra,'.',mismatch,".final_result_readcount.txt"), header = TRUE,  sep="\t")
TEstats <- TEstats0
dim(TEstats)#95 179

##### Filter for those cases that are cross-regulated (automatically exclude the melanogaster subgroup and TE with less reads than DmelTot/700)
TEstatsF<- GetTableCrossReg(TEstats, TRUE, minpercent, minRead)
dim(TEstatsF)
TEstatsF<- TEstatsF[TEstatsF$PiNb>=minRead,]
dim(TEstatsF)

# Save CrossRegulated table
write.table(TEstatsF,paste("AllSpecies/Table/Species.TE.",sra,'.',mismatch,'.',"CrossRegulated.txt",sep=""),quote=FALSE,col.names=TRUE, row.names=FALSE,sep='\t')


###### FILTERING ########
# Remove TE for which max piRNA in Dmel.Ref is less than 1/700 of the total mapping reads
TEstats <- TEstats[sapply(TEstats, function(x) max(x, na.rm = T) > sum(TEstats["Dmel.Ref",])/minRead)]
dim(TEstats) #95 120

# And we want at least minSpe species sharing TE to keep it
TEstats <-TEstats[,which(colSums(!is.na(TEstats))> minSpe)]
dim(TEstats) #95 68

###### PREPARING NORMALIZATION ########
# Reorder columns according to the decreasing number of species bearing copies
TEstats <-TEstats[,order(-colSums(!is.na(TEstats)))]
dim(TEstats) #94 68

# Exclude the ref and put it in ref
Ref<-TEstats["Dmel.Ref",]
TEstats <-TEstats[-which(rownames(TEstats)=="Dmel.Ref"),]
dim(TEstats) #94 68 now 118

# Save Species names for later
Species<-rownames(TEstats)

###################################################################################################
#########.   Combine with the tree
###################################################################################################
# Read tree
tree <- ape::read.tree("Data/BigGene_.tre")

# Keep only those tips which are present in the result dataframe
tree <- keep.tip(tree, c(row.names(TEstats)))

# Define outgroup
outgroup <- c("C.costata", "L.varia")

# Root the tree
tree <- root(tree, outgroup)

# Rotate the tree
tree <- rotateConstr(tree, c("D.melanogaster", tree$tip[!tree$tip %in% c("D.melanogaster",outgroup)], outgroup))

# Get the plotting order
tips.plotordered <- tree$tip[tree$edge[tree$edge[,2] <= length(tree$tip),2]]

# Rearrange dataframe rows according to the tree order
TEstats <- TEstats[match(rev(tips.plotordered), row.names(TEstats)), ] 
Species<-rownames(TEstats)

# Normalize by total number of pi in Dmel.Ref (must be done on df)
vecnorm<-as.vector(Ref)
TEstatsN <- mapply('/', TEstats, vecnorm,SIMPLIFY = TRUE)
TEstatsN<-as.data.frame(TEstatsN)

########################################
########### For constructing the figure 4
########################################
TEstats2 <- TEstatsN

### Color specific groups)
melsubgroup=c("D.melanogaster","D.sechellia","D.mauritiana","D.simulans","D.erecta","D.yakuba","D.teissieri.2733")
colors<-rep("black",94)
names(colors)<-tree$tip
colors[grepl("Z.", names(colors))]<-"#00CCFF"
colors[names(colors) %in% melsubgroup]<-"red"


# Make >1 to 1 to keep range between 0 and 1
TEstats2[TEstats2 > 1] <- 1.0
rownames(TEstats2)<- Species
                          
breaks <- c(0, 0.01, seq(from = 0.1, to = 1, by = 0.3), 1)

p1 <- ggtree(tree, ladderize = FALSE) + geom_tiplab(size = 4.7, align=TRUE, linesize=.5,color=colors)
gheatmap(p1, TEstats2, offset=0.2, width=5, colnames_angle=90, color="white", colnames_offset_y=-3) +
scale_fill_viridis_c(option = "E", direction = -1, na.value = 'gray', name="Normalized piRNA hits") + ggtree::vexpand(.1, -1)+theme(legend.position='top',legend.key.size = unit(1, 'cm'), legend.text = element_text(size = 12),)
ggsave(paste0("AllSpecies/Figures/",sra,".",mismatch,".Fig4.Heatmap.pdf"), width = 50, height = 55, units = "cm", limitsize = FALSE)

print("Figure Done")
