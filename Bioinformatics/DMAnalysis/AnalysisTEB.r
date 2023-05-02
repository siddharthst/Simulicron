# Load the required libs
suppressPackageStartupMessages({
require(tidyverse)
require(ape)
require(phytools)
require(Biostrings)
library(seqinr)}
)

# Make sure blast, seqtk, bedtools, bedops and mafft are installed and accessible via path
# set additional parameters below
genome         <- "~/Documents/Projects/Genomes/Drosophila_melanogaster.BDGP6.32.dna.toplevel.fa"
TEdbOrig       <- "./TEDB.fa"
SmallRNAReads  <- "~/Documents/Projects/SequencingLib/SRR14569563/reads_trimmed.fq"
pathTopiDB     <- "./piCluster.dmel6.bed"

# Distance threshold for merging TE fragments
TEfragmentDistance <- 250

# Minimum consensus length
minimumConLength <- 1500

# Minimum piCluster count
minPiClusterCount = 2

# Minimum number of TE copies
minimumCopyNumber = 6

# Minumum number of piRNA reads required to keep the piRNA insertion
minPiRNAReads = 200

# Create a filtered consensus
TEdbFile <- "TEDBFiltered.fa"
cmd = paste("seqtk seq -L ", minimumConLength, " ", TEdbOrig, " > ", TEdbFile, sep="")
system(cmd)

# Create dirs
dir.create(file.path(getwd(), "temp"), showWarnings = FALSE)
dir.create(file.path(getwd(), "genomeDB"), showWarnings = FALSE)
dir.create(file.path(getwd(), "NCBIBLResults"), showWarnings = FALSE)
dir.create(file.path(getwd(), "BEDFiles"), showWarnings = FALSE)
dir.create(file.path(getwd(), "Results"), showWarnings = FALSE)


# Threshold for minimum length
minTElen <- 0.7

# If only specific TEs are required to be analysed
# Otherwise analyse all
# Uncomment as needed!
# TEsToAnalyze <- c("IDM", "GYPSY", "COPIA", "BATUMI")
TEsToAnalyze <- names(read.fasta(file = TEdbFile, seqtype = "DNA"))
# TEsToAnalyze <- c("I_DM_I", "COPIA_DM_Copia")

###################-----------------------------------------------------------------###################

# First build a blast database
cmd = paste("makeblastdb -in ", genome, " -dbtype nucl -out ./genomeDB/genomeDB", sep="")
system(cmd)

# Now run BLAST on the sequences
cmd = paste("blastn -db ./genomeDB/genomeDB -query ", TEdbFile, " -outfmt 6 -evalue 20 -out ./NCBIBLResults/results.tab", sep="")
system(cmd)

# Convert the output to a BED file
# https://github.com/nterhoeven/blast2bed
cmd = "bash blast2bed.sh ./NCBIBLResults/results.tab"
system(cmd)

# Copy the resulting file to BED folder
fs::file_copy("./NCBIBLResults/results.tab.bed", "./BEDFiles/results.tab.bed", overwrite = TRUE)

# Sort the bed file in proper order
cmd = "sort-bed ./BEDFiles/results.tab.bed > ./BEDFiles/results.tab.bed.sorted"
system(cmd)
unlink("./BEDFiles/results.tab.bed")
file.rename("./BEDFiles/results.tab.bed.sorted", "./BEDFiles/results.tab.bed")

# Create lists for storing statistics
RecordTEName              <- c()
RecordTELength            <- c()
RecordTECopies            <- c()
RecordTECopiesAfterFilter <- c()
RecordTEFiltered          <- c()
RecordPiClusters          <- c()
RecordPiFiltered          <- c()
RecordFullLengthPi        <- c()

# Working on each TE individually 
ConsensusSequence <- readDNAStringSet(TEdbFile)
for(TE in TEsToAnalyze){
    FullLengthCopyNumber <- 0
    FullLengthInPi       <- 0
    RecordTEName <- c(RecordTEName, TE)
    dir.create(paste(file.path(getwd(), "Results/"), TE, "/",sep=""), showWarnings = FALSE)
    TELocalDIR = paste(file.path(getwd(), "Results/"), TE, "/",sep="")
    cmd = paste("grep \t", TE, "\t ./BEDFiles/results.tab.bed | bedtools sort > ", TELocalDIR, TE, ".bed", sep="")
    system(cmd)
    
    # Merge TEs
    cmd = paste("bedtools merge -s -d ", TEfragmentDistance, " -c 4,5,6 -o distinct,mean,distinct -i ", TELocalDIR, TE, ".bed > ", TELocalDIR, TE, ".merged.bed", sep="")
    system(cmd)
    
    # Fix the name column
    bedFrame <- read.table(file = paste(TELocalDIR, TE, ".merged.bed", sep=""), sep = '\t', header = FALSE)
    bedFrame$V4 <- seq.int(nrow(bedFrame))
    write.table(bedFrame, file = paste(TELocalDIR, TE, ".merged.bed", sep=""), row.names=FALSE, col.names=FALSE, sep="\t", quote = FALSE)
    
    # Check which TEs overlap with piRNA
    cmd = paste("bedtools intersect -wo -a ", TELocalDIR, TE, ".merged.bed -b ./piCluster.dmel6.bed > ", TELocalDIR, TE, ".merged.piRNA.bed", sep="")
    system(cmd)
    
    # Identify piRNA insertions and label them as such
    fileSize = file.info(paste(TELocalDIR, TE, ".merged.piRNA.bed", sep=""))$size
    if (fileSize == 0)
        {
        RecordPiClusters <- c(RecordPiClusters, 0)
        RecordPiFiltered <- c(RecordPiFiltered, 0)
        RecordTELength <- c(RecordTELength, 0)
        RecordTECopies <- c(RecordTECopies, 0)
        RecordFullLengthPi <- c(RecordFullLengthPi,0)
        RecordTECopiesAfterFilter <- c(RecordTECopiesAfterFilter, 0)
        RecordTEFiltered <- c(RecordTEFiltered, "Yes/NPi")
       next 
    }
    piFrame <- read.table(file = paste(TELocalDIR, TE, ".merged.piRNA.bed", sep=""), sep = '\t', header = FALSE)
    bedFrame <- read.table(file = paste(TELocalDIR, TE, ".merged.bed", sep=""), sep = '\t', header = FALSE)
    listOfPiInsertions <- c(unique(piFrame$V4))
    for (piInsertions in listOfPiInsertions){
        bedFrame$V4[piInsertions] <- paste(piInsertions, "_piRNA", sep="")
    }
    write.table(bedFrame, file = paste(TELocalDIR, TE, ".merged.bed", sep=""), row.names=FALSE, col.names=FALSE, sep="\t", quote = FALSE)
    
    # Extract the sequences
    cmd = paste("bedtools getfasta -nameOnly -s -fi ", genome, " -bed ", TELocalDIR, TE, ".merged.bed > ", TELocalDIR, TE, ".fa", sep="")
    system(cmd)
    
    # Perform MSA
    ConsensusLength <- nchar(toString(ConsensusSequence[TE]))
    RecordTELength <- c(RecordTELength, ConsensusLength)
    minReqLength <- minTElen * ConsensusLength
    fastaLocation = paste(TELocalDIR, TE, ".fa", sep="")
    FastaSubset <- readDNAStringSet(fastaLocation)
    RecordTECopies <- c(RecordTECopies, length(FastaSubset))
    MSAset <- c(ConsensusSequence[TE])
    Alignmentset <- c(ConsensusSequence[TE])
    outputNames <- paste(TELocalDIR, TE, ".msa",sep="")
    tempFasta = tempfile(pattern = TE, tmpdir = "./temp")
    tempMSA = tempfile(pattern = TE, tmpdir = "./temp")
    for (i in 1:length(FastaSubset)){
        # Write the object into fasta file
        writeXStringSet(c(ConsensusSequence[TE],FastaSubset[i]), tempFasta, append=FALSE, compress=FALSE, format="fasta")
        # Perform MSA using mafft
        cmd = paste("mafft --auto --preservecase ",  tempFasta," > ", tempMSA, sep="")
        system(cmd)
        # Read the results into another DNAStringSet
        MSA <- readDNAStringSet(tempMSA)
        # For insertions in TE relative to consensus
        InsertionSites <- gregexpr(pattern ='-',MSA[1])
        if (unlist(InsertionSites)[1] == -1)
            {
            # Nothing to do
        }
        else{
            # Delete the insetions in genomic TE copy
            MSA[[2]] <- MSA[[2]][-c(unlist(InsertionSites))]
            # Now delete the gaps in consenus
            MSA[[1]] <- MSA[[1]][-c(unlist(InsertionSites))]
        }
        # Append the results to MSA if it follows the min length criteria
        if (nchar(gsub("-", "", MSA[[2]])[1]) > minReqLength)
            {
                FullLengthCopyNumber <- FullLengthCopyNumber + 1
                MSAset <- append(MSAset, MSA[2])
                if (grepl("piRNA", names(MSA[2])[1]) == TRUE)
                {
                    FullLengthInPi <- FullLengthInPi + 1
                }
                # Append the results to Alignment file by replacing - with N
                # MSA[[2]] <- gsub("-", "N", MSA[[2]])
                # Alignmentset <- append(Alignmentset, MSA[2])
        }
        else if (grepl("piRNA", names(MSA[2])[1]) == TRUE)
        {
            MSAset <- append(MSAset, MSA[2])  
        }
        # Remove the temp files
        unlink(tempFasta)
        unlink(tempMSA)
    }
    if (FullLengthCopyNumber > minimumCopyNumber)
        {
            MSAset <- MSAset[2:length(MSAset)]
            RecordTECopiesAfterFilter <- c(RecordTECopiesAfterFilter, FullLengthCopyNumber)
            writeXStringSet(MSAset, outputNames, append=FALSE, compress=FALSE, format="fasta")
            # Before generating trees, we need to rename the fasta headers as fasttree doesn't like :
            # cmd = paste("seqtk rename ", TELocalDIR, TE, ".msa > ", TELocalDIR, TE, ".msa.fixed", sep="")
            # system(cmd)
            # Generate the tree
            cmd = paste("fasttree -nt -gtr < ", TELocalDIR, TE, ".msa > ", TELocalDIR, TE, ".tree", sep="")
            system(cmd)
            # Load the tree and plot it
            tree <- read.tree(paste(TELocalDIR, TE, ".tree", sep=""))
            tree <- midpoint.root(tree)
            pdf(paste(TELocalDIR, TE, ".pdf", sep=""))
            plot(tree)
            add.scale.bar()
            dev.off()
            
            # Perform smallRNA alignement 
            # Creating a bowtie index
            cmd=paste("bowtie-build -f ", TELocalDIR, TE, ".msa ./temp/TE", sep="")
            system(cmd)
            # Perform alignment
            # Align
            pathToAlignment    = paste(TELocalDIR, TE, ".sam", sep="")
            pathToBAM          = paste(TELocalDIR, TE, ".bam", sep="")
            pathToAlignmentBED = paste(TELocalDIR, TE, ".Alignment.bed", sep="")
            cmd = paste("bowtie -a -v 0 -p 8 --sam --no-unal -q -x ./temp/TE ", SmallRNAReads, " > ", pathToAlignment, sep="")
            system(cmd)
            # Converting SAM to BAM
            cmd = paste("samtools sort ", pathToAlignment, " > ", pathToBAM, sep="")
            system(cmd)
            # Generating index
            cmd = paste("samtools index ", pathToBAM, sep="")
            system(cmd)
            # Converting BAM to BED12
            cmd = paste("bedtools bamtobed -i ", pathToBAM, " > ",  pathToAlignmentBED, sep="")
            system(cmd)
            # Read the file into dataframe
            AlignmentDataFrame         <- as.data.frame(read.table(pathToAlignmentBED, header = FALSE, sep="\t",stringsAsFactors=FALSE, quote=""))
            AlignmentDataFrame$V1      <- gsub("\\(-\\)","",as.character(AlignmentDataFrame$V1))
            AlignmentDataFrame$V1      <- gsub("\\(\\+\\)","",as.character(AlignmentDataFrame$V1))
            ReadsDataFrame             <- data.frame(matrix(nrow = length(tree$tip.label), ncol = length(c(grep("piRNA", tree$tip.label, value = TRUE)))))
            colnames(ReadsDataFrame)   <- grep("piRNA", tree$tip.label, value = TRUE)
            rownames(ReadsDataFrame)   <- tree$tip.label

            # Find shared reads between two loci
            for (piRNAName in colnames(ReadsDataFrame)){
                for (TEName in rownames(ReadsDataFrame)) {
                    subsetDataPiRNA          <- AlignmentDataFrame[AlignmentDataFrame$V1 == piRNAName,]
                    subsetDataPiRNAReads     <- subsetDataPiRNA$V4
                    subsetDataTarget         <- AlignmentDataFrame[AlignmentDataFrame$V1 == TEName,]
                    subsetDataTargetpiRNA     <- subsetDataTarget[subsetDataTarget$V4 %in% subsetDataPiRNAReads,]
                    ReadsDataFrame[TEName, piRNAName] = NROW(subsetDataTargetpiRNA)
                }
            } 
            
            # Remove piRNA insertions with low number of reads
            ReadsDataFrame <- ReadsDataFrame[sapply(ReadsDataFrame, function(x) max(x, na.rm = T) > minPiRNAReads)]
            
            # Continue if no piRNA insertion left
            checkPIRNA <- length(grep("piRNA", colnames(ReadsDataFrame), value = TRUE))
            originalPiClusters  <- grep("piRNA", tree$tip.label, value = TRUE)
            remainingPiClusters <- grep("piRNA", colnames(ReadsDataFrame), value = TRUE)
            RecordPiClusters <- c(RecordPiClusters, length(originalPiClusters))
            RecordPiFiltered <- c(RecordPiFiltered, length(remainingPiClusters))
            if (checkPIRNA == 0){
                RecordFullLengthPi <- c(RecordFullLengthPi, FullLengthInPi)
                RecordTEFiltered <- c(RecordTEFiltered, "Yes/NoClusterAfterFilter")
                file.create(paste(TELocalDIR, "Failed.TE",sep=""))
                next
            }
            else if (checkPIRNA < minPiClusterCount){
                RecordFullLengthPi <- c(RecordFullLengthPi, FullLengthInPi)
                RecordTEFiltered <- c(RecordTEFiltered, "Yes/MinPiCLuster")
                file.create(paste(TELocalDIR, "Failed.TE",sep=""))
                next
            }  
                                                    
            # Also remove those respective rows and leafs! 
            if (length(originalPiClusters) != length(remainingPiClusters)){
                removedClusters <- c(setdiff(remainingPiClusters, originalPiClusters), setdiff(originalPiClusters, remainingPiClusters))
                ReadsDataFrame <- ReadsDataFrame[!(row.names(ReadsDataFrame) %in% removedClusters),]
                tree <- ape::drop.tip(tree, removedClusters)
            }
            RecordFullLengthPi <- c(RecordFullLengthPi, FullLengthInPi)
            RecordTEFiltered <- c(RecordTEFiltered, "No")
            
                                                    
            # Plot using Arnaud's method
            # Assuming no normalisation is applied 
            pdf(paste(TELocalDIR, TE, ".plot.pdf"), width=10, height=10)

            tree <- phytools::midpoint.root(unroot(tree))
            ncols <- 101 # number of colors for the image display
            col.pi <- "darkred"
            col.te <- "darkblue"           
                                                    
            piRNA.counts <- ReadsDataFrame
            piRNAs <- grepl("piRNA", tree$tip.label)
            tree <- rotateConstr(tree, c(tree$tip.label[!piRNAs], tree$tip.label[piRNAs]))

            # We need the same order in both figures. Not a better way to get the order from the tree? 
            tip.ordered <- tree$tip.label[tree$edge[tree$edge[,2] <= length(tree$tip),2]]
            piRNA.counts <- piRNA.counts[tip.ordered, rev(tip.ordered[grep("piRNA", tip.ordered)]), drop=FALSE]

            norm.counts <- piRNA.counts / max(piRNA.counts)
            norm.counts[grep("piRNA", rownames(norm.counts)),] <- -norm.counts[grep("piRNA", rownames(norm.counts)),]-1e-6
            
            # Legend
            layout(cbind(c(1,1),c(2,2),c(3,4)), widths=c(1,1,0.3))
            par(mar=c(1, 1, 6, 1), oma=c(0,0,2,0), cex=1)
            plot(tree,
            x.lim           = 1.25*max(node.depth.edgelength(tree)),
            y.lim           = c(0.5, length(tree$tip.label)+0.5), 
            yaxs            = "i",
            tip.color       = ifelse(grepl("piRNA", tree$tip.label), col.pi, col.te), 
            align.tip.label = 0)
            ape::add.scale.bar()
            zlim <- range(norm.counts) + c(-0.01, 0.01)
            
            image(
                x      = seq_along(colnames(norm.counts)), 
                y      = seq_along(rownames(norm.counts)), 
                z      = t(norm.counts), 
                axes   = FALSE, 
                xlab   = "", 
                ylab   = "", 
                zlim   = zlim, 
                col    = c(colorRampPalette(c(col.pi,gray(0.8)))(floor(ncols/2)), colorRampPalette(c("white", col.te))(ncols-floor(ncols/2))), 
                breaks = c(seq(min(zlim), -1e-6, length=floor(ncols/2+1)), seq(0, max(zlim), length=ncols-floor(ncols/2))))

            axis(3, at=seq_along(colnames(norm.counts)), colnames(norm.counts), las=2, tick=FALSE)
            
            yy.pi <- seq(0, -zlim[1]*max(piRNA.counts), length.out=floor(ncols/2)+1)
            yy.te <- seq(0, zlim[2]*max(piRNA.counts), length.out=floor(ncols/2)+1)

            par(mar=c(1,1,6,4))

            image(
                y= yy.pi,
                z=t(as.matrix(yy.pi)), 
                xlab="", ylab="", axes=FALSE,
                col=colorRampPalette(c(gray(0.8),col.pi))(floor(ncols/2))
            )
            axis(4)
                
            mtext("pi-Cluster hits", 4, line=2.5)
                
            image(
                y= yy.te,
                z=t(as.matrix(yy.te)), 
                xlab="", ylab="", axes=FALSE,
                col=colorRampPalette(c("white",col.te))(floor(ncols/2))
            )
            axis(4)
                
            mtext("TE hits", 4, line=2.5)
                
            title(TE, outer=TRUE)
                
            dev.off()
                                                    
            # Write the data
            write.table(ReadsDataFrame, file=(paste(TELocalDIR, TE, ".tsv", sep="")), quote=FALSE, sep='\t')
            scaledReadsDataFrame <- apply(ReadsDataFrame, 2, function(x) x / max(x))
            write.table(scaledReadsDataFrame, file=(paste(TELocalDIR, TE, ".max.norm.tsv", sep="")), quote=FALSE, sep='\t')
    }
    else
        {
        RecordTECopiesAfterFilter <- c(RecordTECopiesAfterFilter, FullLengthCopyNumber)
        RecordTEFiltered <- c(RecordTEFiltered, "Yes/MinimumTEFailed")
        RecordPiClusters <- c(RecordPiClusters, 0)
        RecordPiFiltered <- c(RecordPiFiltered, 0)
        RecordFullLengthPi <- c(RecordFullLengthPi, 0)
        file.create(paste(TELocalDIR, "Failed.TE",sep=""))
    }
    # Delete all files in temp folder
    f <- list.files("./temp", include.dirs = F, full.names = T, recursive = T)
    # remove the files
    file.remove(f)
}

# Create a dataframe with statistics
statsFrame <- data.frame(RecordTEName, RecordTELength, RecordTECopies, RecordTECopiesAfterFilter, RecordFullLengthPi, RecordTEFiltered, RecordPiClusters, RecordPiFiltered)

# Save the dataframe on file
write.table(statsFrame, file="./Results/stats.tsv", quote=FALSE, sep='\t', row.names=FALSE)

####
