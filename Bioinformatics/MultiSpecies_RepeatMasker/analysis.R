# Load required packages
suppressPackageStartupMessages({
library(R.utils)
library(tidyverse)
library(Biostrings)
library(GenomicRanges)
library(GenomicAlignments)
library(phylotools)
library(stringr)
library(ggtree)
library(gplots)
library(superheat)
library(ggtreeExtra)
library(aplot)
library(ape)
library(adephylo)
library(phylobase)
library(phylosignal)
library(phytools)
library(foreach)
library(RColorBrewer)
library(doParallel)
library(bio3d)
library(DECIPHER)
library(purrr)
})

# Prepare core files 
# Parameters
urlToSmallRNA = "https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos4/sra-pub-run-25/SRR14569563/SRR14569563.1"
smallRNASRAid = "SRR14569563"
# Create folders for storing reads
mainDir <- getwd()
dir.create(file.path(getwd(), "Reads"), showWarnings = TRUE)
dir.create(file.path(getwd(), "temp"), showWarnings = TRUE)
# Getting the SRA archives for small RNA data
print("Getting small-RNA fastq")
cmd=paste("fasterq-dump", "--threads 8", smallRNASRAid, "-o", "./Reads/reads.fq", sep=' ')
system(cmd)

# Performing adapter trimming
print("Performing adapter trimming")
cmd=paste("trim_galore ./Reads/reads.fq -o ./Reads")
system(cmd)

# Performing de-duplication of sequences
print("Performing de-duplication and dusting")
# cmd=paste("seqkit rmdup -s < ./Reads/reads_trimmed.fq > ./Reads/reads_deDup.fq")
cmd=paste("fqtrim -D -C ./Reads/reads_trimmed.fq -o reads_collapsed.fq --outdir ./Reads")
system(cmd)

# Performing tRNA, rRNA and miRNA filtering
print("Performing tRNA, rRNA and miRNA filtering")
#cmd=paste("bowtie-build -f filters.fa ./temp/filters")
#system(cmd)
#cmd=paste("bowtie -a -p 8 --sam -q --un ./Reads/reads_filtered.fq -x ./temp/filters ./Reads/reads_trimmed.fq > ./temp/filterCatch.sam")
#system(cmd)
cmd="bbduk.sh k=18 minlength=18 maxlength=30 in=./Reads/reads_trimmed.fq out=./Reads/reads_filtered.fq ref=filters.fa"
#cmd="bbduk.sh k=18 minlength=18 maxlength=30 in=./Reads/reads_trimmed.reads_collapsed.fq out=./Reads/reads_filtered.fq ref=filters.fa"
system(cmd)

# Run repeat masker on the TE file
print ("Running RepeatMasker with included TE.fa as library")
dir.create(file.path(getwd(), "RepeatMasker"), showWarnings = TRUE)
cmd = "RepeatMasker -pa 8 -s -gff -no_is -nolow -norna -div 40 -dir ./RepeatMasker -lib ./TE.fa ./Genomes/*.fasta"
system(cmd)

# Convert the generated GFF file to BED file using bedops
print ("Converting RepeatMasker GFF to BED")
cmd = "bash ./GFFtoBED.sh"
system(cmd)

# Read TE fasta file into dataframe and calculate the length of each TE
TEDataBase <- readDNAStringSet("TE.fa")
seq_name = names(TEDataBase)
sequence = paste(TEDataBase)
TEDataBaseFrame <- data.frame(seq_name, sequence)file.exists("leaflet.R")
TEDataBaseFrame$LengthDB <- str_count(TEDataBaseFrame$sequence)

# Process the output
dir.create(file.path(getwd(), "Results"), showWarnings = FALSE)
dir.create(file.path(getwd(), "Shortlist"), showWarnings = FALSE)
RepeatMaskerFiles <- list.files("./RepeatMasker/")
bedFiles <- grep(".*bed", RepeatMaskerFiles, value=T)
for (i in bedFiles){
    # Load the file into R
    TEtable <- read.table(file = paste("./RepeatMasker/",i, sep=""))
    # List of unique TE records
    TEs <- as.list(unique(TEtable[c("V11")]))$V11
    # Create an empty list to store results 
    resultList <- list()
    # Create empty table to store results
    ResultTable <- data.frame()
    # Strings to remove from RepeatMasker results
    deletionString <- 'Motif:'
    # Subset bed file based on TE and then process each record individually
    for (TE in TEs){
        tempTable <- subset(TEtable, V11==TE)
        write.table(tempTable, file="./temp/merge.tab", quote=FALSE, sep='\t', col.names = FALSE, row.names = FALSE)
        # Merging TE fragments which are closs enough (enough = 250bp)
        # print (paste("Merging fragmented genomic copies for :", TE))
        cmd = ("bedtools merge -s -d 250 -c 4,5,6 -o distinct,mean,distinct -i ./temp/merge.tab > ./temp/merged.tab")
        system(cmd)
        # Reading the merged file
        tempTable <- read.table(file = './temp/merged.tab')
        tempTable["V4"] <- gsub(deletionString, "", TE)
        ResultTable <- rbind(ResultTable, tempTable)
        write.table(ResultTable, file=paste("./Results/All_TEs_",i, sep=""), quote=FALSE, sep='\t', col.names = FALSE, row.names = FALSE)
    }
    # Read file into datframe
    TEoverlapTable <- read.table(file = paste("./Results/All_TEs_",i, sep=""))
    TEoverlapTable["Length_TE"] <- TEoverlapTable["V3"] - TEoverlapTable["V2"]
    TEoverlapTable["seq_name"] <- as.character(TEoverlapTable[,4])
    # Check for length 
    LengthCutOff <- 0.7
    # Add information about length
    TEoverlapTable <- TEoverlapTable %>% left_join(TEDataBaseFrame)
    TEoverlapTable <- subset(TEoverlapTable, select = -c(sequence))
    TEoverlapTable <- TEoverlapTable[!(TEoverlapTable$Length_TE < (TEoverlapTable$LengthDB * LengthCutOff)),]
    write.table(TEoverlapTable, file=paste("./Shortlist/Shortlisted_TEs_",i, sep=""), quote=FALSE, sep='\t', col.names = FALSE, row.names = FALSE)
    # Create another frame
    TEhistogramFrame <- as.data.frame(table(TEoverlapTable["V4"]))
    # Split the transposons 
    spl <- strsplit(as.character(TEhistogramFrame$Var1), "_")
    TEhistogramFrame["Superfamily"] <- sapply(lapply(spl, tail, 1), paste, collapse="_")
    TEhistogramFrame["Family"] <- sapply(lapply(spl, head, -1), paste, collapse="_")
    # Remove DM from Family
    TEhistogramFrame %>% mutate_all(~gsub("_DM", "", .))
    write.table(TEhistogramFrame, file=paste("./Results/All_Tables_",i, sep=""), quote=FALSE, sep='\t', col.names = TRUE, row.names = FALSE)
}

# Run the python script to find the top 10 shared TEs
cmd = "python PythonAnalysis.py"
system(cmd)

# For the shortlisted TEs, run grep on all files
cmd = "bash GREP.sh"
system(cmd)

# Extract fasta files 
dir.create(file.path(getwd(), "Fasta"), showWarnings = FALSE)
# Get the base file names for the genome files
bed2Records <- list.files("./Shortlist/")
bed2Files <- grep(".*bed2", bed2Records, value=T)
bed2Files <- str_remove(bed2Files, '^Shortlisted_TEs_')
hostGenomes <- sub('\\.out.bed.bed2$', '', bed2Files) 
for (i in hostGenomes){
    #cmd =  paste("perl", "bed2faidxsta.pl", "--fastaIsUncompressed", paste("--fastaDir=", "./Genomes/",i, " < ./Shortlist/Shortlisted_TEs_", i ,".out.bed.bed2 > ./Fasta/", i,sep=""), 
    #             sep=" ")
    cmd = paste("bedtools getfasta -name -s -fi", paste("./Genomes/", i, " -fo ./Fasta/", i, " -bed ./Shortlist/Shortlisted_TEs_", i, ".out.bed.bed2",sep="") ,sep=" ")
    system(cmd)
    # print (cmd)
    }

# With the fasta files ready
# align the short RNA sequence
# to all the files seperately 
# But before that, correct files 
# Delete any file which is empty
cmd = "find ./Fasta/ -size 0 -delete"
system(cmd)
# Convert the files to uppercase
fastaFiles <- list.files("./Fasta/")
for (i in fastaFiles){
    cmd = paste("reformat.sh in=./Fasta/", i, " out=./Fasta/", i, ".fa", " tuc -Xmx1g", sep="")
    system(cmd)
}

# Setup a parallel backend
cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)
options(warn=-1)

# Setup variables
fastaFiles <- list.files("./Fasta/")
fastaFaFiles <- grep(".fa$", fastaFiles, value=TRUE)
# fastaFaFiles <- c("D.mauritiana.fasta.fa")
shortListedTEs <- scan("shortlisted_TEs.txt", character(), quote = "")
ConsensusSequence <- readDNAStringSet("./TE.fa") 
ConsensusSequence <- ConsensusSequence[shortListedTEs]

fastaFiles <- list.files("./Fasta/")
fastaFaFiles <- grep(".fa$", fastaFiles, value=TRUE)
# fastaFaFiles <- c("D.mauritiana.fasta.fa")
shortListedTEs <- scan("shortlisted_TEs.txt", character(), quote = "")
ConsensusSequence <- readDNAStringSet("./TE.fa") 
ConsensusSequence <- ConsensusSequence[shortListedTEs]
foreach(fastaI=fastaFaFiles, .errorhandling="pass") %dopar%{
    library(Biostrings)
# for(fastaI in fastaFaFiles){
    print(fastaI)
    fastaLocation = paste("./Fasta/", fastaI, sep="")
    FastaSequence <- readDNAStringSet(fastaLocation) 
    MSAset <- c(ConsensusSequence[1])
    Alignmentset <- c(ConsensusSequence[1])
    outputNames <- paste("./Fasta/", fastaI, ".msa",sep="")
    tempFasta = tempfile(pattern = fastaI, tmpdir = "./temp")
    tempMSA = tempfile(pattern = fastaI, tmpdir = "./temp")
    if (file.exists(outputNames) == TRUE){
        next
    }
    for(ConsensusI in 1:length(ConsensusSequence)){
        Consensus <- (names(ConsensusSequence[ConsensusI]))[1]
        isMatching <- grep(Consensus, names(FastaSequence), value=TRUE)
        if(length(isMatching) == 0){
            next
        }
        FastaSubset <- FastaSequence[sapply(Consensus, function(x) grep(x, names(FastaSequence)))]
        for (i in 1:length(FastaSubset)){
            # Write the object into fasta file
            writeXStringSet(c(ConsensusSequence[ConsensusI],FastaSubset[i]), tempFasta, append=FALSE, compress=FALSE, format="fasta")
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
            # Append the results to MSA
            MSAset <- append(MSAset, MSA[2])
            # Append the results to Alignment file by replacing - with N
            MSA[[2]] <- gsub("-", "N", MSA[[2]])
            Alignmentset <- append(Alignmentset, MSA[2])
            # Remove the temp files
            unlink(tempFasta)
            unlink(tempMSA)
            }        
    }
    MSAset <- MSAset[2:length(MSAset)]                                        
    writeXStringSet(MSAset, outputNames, append=FALSE, compress=FALSE, format="fasta")
}

options(warn=-1)
fastaFiles <- list.files("./Fasta/")
fastaFaFiles <- grep(".fa.msa$", fastaFiles, value=TRUE)
shortListedTEs <- scan("shortlisted_TEs.txt", character(), quote = "")
ConsensusSequence <- readDNAStringSet("./TE.fa") 
ConsensusSequence <- ConsensusSequence[shortListedTEs]
# fastaFaFiles <- c("D.mauritiana.fasta.fa")
foreach(fastaI=fastaFaFiles) %dopar%{
    library(tidyverse)
    library(Biostrings)
    library(ape)
    library(bio3d)
    library(DECIPHER)
    # https://gist.github.com/joelnitta/6f30a7c0f1c83d78c76a5469e935d56f
    DNAbin_to_DNAstringset <- function (seqs, remove_gaps = FALSE) {
      if(isTRUE(remove_gaps)) {
      seqs %>% as.list() %>% as.character %>% 
          lapply(.,paste0,collapse="") %>% 
          lapply( function (x) gsub("-", "", x)) %>% 
          unlist %>% Biostrings::DNAStringSet()
      } else {
        seqs %>% as.list() %>% as.character %>% 
          lapply(.,paste0,collapse="") %>% 
          unlist %>% Biostrings::DNAStringSet()
      }
    }
    #----------------------------------------------
    fastaLocation = paste("./Fasta/", fastaI, sep="")
    FastaSequence <- readDNAStringSet(fastaLocation) 
    MSAset <- c(ConsensusSequence[1])
    outputNames <- paste("./Fasta/", fastaI, ".msa.cor",sep="")
    if (file.exists(outputNames) == TRUE){
        next
    }
    for(ConsensusI in 1:length(ConsensusSequence)){
        Consensus <- (names(ConsensusSequence[ConsensusI]))[1]
        isMatching <- grep(Consensus, names(FastaSequence), value=TRUE)
        if(length(isMatching) == 0){
            next
        }
        FastaSubset <- FastaSequence[sapply(Consensus, function(x) grep(x, names(FastaSequence)))]
        # Define base identity
        baseIdentity = 0.5
        tempFasta = tempfile(pattern = fastaI, tmpdir = "./temp")
        # Find ungapped length
        sequenceLengths <- width(DNAbin_to_DNAstringset(del.gaps(c(ConsensusSequence[ConsensusI], FastaSubset))))
        consensusLength <- sequenceLengths[1]
        # Convert to %
        sequencePerctns <- sequenceLengths/consensusLength
        # Keep only 80% length
        subset <- c(ConsensusSequence[ConsensusI], FastaSubset)[sequencePerctns>0.80]
        if(length(subset)==1){
            next
        }
        # Write the object into fasta file
        # writeXStringSet(subset, tempFasta, append=FALSE, compress=FALSE, format="fasta")
        # Perform identity calculation
        # mat <- bio3d::seqidentity(bio3d::read.fasta(tempFasta))
        mat <- DistanceMatrix(subset,type = "matrix",
                       includeTerminalGaps = TRUE,
                       penalizeGapLetterMatches = TRUE,
                       penalizeGapGapMatches = TRUE,
                       correction = "none",
                       processors = 1,
                       verbose = FALSE)
        mat <- 1 - mat
        diag(mat) <- 0
        if(mat[max.col(mat)[1],1] > baseIdentity){
            # Instead of picking one, pick all 
            matZ <- mat
            matZ <- matZ[order(matZ[,1],decreasing=TRUE),]
            rownames(matZ)[matZ[,1]>baseIdentity]
            MSA <- c(ConsensusSequence[ConsensusI], FastaSubset)[as.character(rownames(matZ)[matZ[,1]>baseIdentity])]
            # MSA <- c(ConsensusSequence[ConsensusI], FastaSubset)[as.character(rownames(mat)[max.col(mat)[1]])]
        }
        else{
            MSA <- c()
        }
        # Remove the temp files
        unlink(tempFasta)
        if (length(MSA)!=0){
        # Append the results to MSA
        MSAset <- append(MSAset, MSA)
        }
    }
    if(length(MSAset)>1){
        MSAset <- MSAset[2:length(MSAset)]                                     
        writeXStringSet(MSAset, outputNames, append=FALSE, compress=FALSE, format="fasta")   
    }
}

# Create a condensed MSA for all TE families
fastaFiles <- list.files("./Fasta/")
fastaFaFiles <- grep(".fa.msa.msa.cor$", fastaFiles, value=TRUE)
ConsensusNames <- names(ConsensusSequence[shortListedTEs])
for(Consensus in (ConsensusNames)){
    outputNames <- paste("./Results/", Consensus, ".msa",sep="")
    MSAset <- c(ConsensusSequence[Consensus])
    for (fasta in fastaFaFiles){
        speciesName <- str_remove(fasta, ".fasta.fa.msa.msa.cor")
        fastaLocation = paste("./Fasta/", fasta, sep="")
        FastaSequence <- readDNAStringSet(fastaLocation)
        isMatching <- grep(Consensus, names(FastaSequence), value=TRUE)
        if(length(isMatching) == 0){
            next
        }
        # Grep for the consensus
        FastaSubset <- FastaSequence[sapply(Consensus, function(x) grep(x, names(FastaSequence)))]
        subsetNames <- names(FastaSubset)
        names(FastaSubset) <- paste(speciesName, subsetNames, sep="_")
        MSAset <- append(MSAset, FastaSubset)                                    
    }
    writeXStringSet(MSAset, outputNames, append=FALSE, compress=FALSE, format="fasta")
}


# Create bowtie indexes and start aligning
dir.create(file.path(getwd(), "Alignment"), showWarnings = FALSE)
fastaFiles <- list.files("./Fasta/")
fastaFiles <- grep("*.fa.msa.msa.cor$", fastaFiles, value=T)

# Create bowtie indexes and start aligning
for (i in fastaFiles){
    cmd = "rm ./temp/*"
    system(cmd)
    cmd = paste("bowtie-build -f ./Fasta/", i ," ./temp/TE", sep="")
    system(cmd)
    cmd = paste("bowtie -a -v 0 -p 8 --sam --no-unal -q -x ./temp/TE", "./Reads/reads_filtered.fq", paste("> ./Alignment/", i, ".sam", sep=""), sep=" ")
    system(cmd)
}

# for (i in fastaFiles){
#     # Converting SAM to BAM and indexing
#     cmd = paste("samtools sort", paste(" ./Alignment/", i, ".sam", sep=""), paste(" > ./Alignment/", i, ".bam", sep=""), sep="")
#     system(cmd)
#     cmd = paste("samtools index", paste(" ./Alignment/", i, ".bam", sep=""), sep="")
#     system(cmd)
# }

# Do further analysis in Python
cmd = "python ./PythonTEAnalysis.py"
system(cmd)

# Read the results into table
piRNAstats <- read.table("./Results/CoreSet.tsv")
# Work with TE family names
colnames(piRNAstats) <- gsub("_.*","",colnames(piRNAstats))
# Make it tidy
piRNATable <- piRNAstats %>% 
  as.data.frame() %>%
  rownames_to_column("Species") %>%
  pivot_longer(-c(Species), names_to = "TE", values_to = "piRNA") 
# Remove the file extension
# piRNATable$Species <- tools::file_path_sans_ext(piRNATable$Species)
piRNATable <- as.data.frame(piRNATable)

# For all species with no piRNA hits - find if there exists a potentially active copy or not
TEfamilyNames <- colnames(piRNAstats)
SpeciesNames <- rownames(piRNAstats)
TEMSAfileNames <- list.files("./Results/")
TEMSAfileNames <- TEMSAfileNames[!grepl('All_Tables_|All_TEs_', TEMSAfileNames)]
for(species in SpeciesNames){
    for(family in TEfamilyNames){
        piRNAcount <- piRNAstats[species, family]
        if(piRNAcount == 0){
            TEmsaFile <- grep(paste("^", family, "_", sep=""), TEMSAfileNames, value = TRUE, ignore.case = FALSE)
            tempTEmsaFile <- readChar(paste("./Results/", TEmsaFile, sep=""), file.info(paste("./Results/", TEmsaFile, sep=""))$size)
            if(grepl(species,tempTEmsaFile)[1]==FALSE){
                piRNAstats[species, family] = NA
            }
        }   
    }
}

# Read the information about TE presence and absence in table
TEstats <- read.table("./Logical.tsv", header = TRUE, row.names = 1)
TEstats <- t(TEstats)
colnames(TEstats) <- gsub("_.*","",colnames(TEstats))
TEstats <- as.data.frame(TEstats)
# Subset the dataframe
TEstats <- TEstats[, colnames(piRNAstats)]
TEstats <- TEstats[rownames(piRNAstats),]
# Set the colors
# TEstats[TEstats != 0] <- "black"
# TEstats[TEstats == 0] <- "white"

# Save the tables
write.table(TEstats, "TEstats.txt", sep = ",", quote = FALSE)
write.table(piRNAstats, "piRNATable.txt", sep = ",", quote = FALSE)

piRNA <- piRNAstats

# Read the tree
tree <- ape::read.tree("./BigGene_.tre")# Remove labels not present
#tree <- keep.tip(tree, c(row.names(piRNAstats)))
species.names.fix <- c(`Drosophila.willistoni.17`="D.willistoni")
outgroup <- c("C.costata", "L.varia")
tree$tip.label[match(names(species.names.fix), tree$tip.label)] <- species.names.fix
tree <- root(tree, outgroup)
rtree <- rotateConstr(tree, c("D.melanogaster", tree$tip[!tree$tip %in% c("D.melanogaster",outgroup)], outgroup))
#plot(rtree, align.tip.label=TRUE, cex=0.6, y.lim=c(0+4.5, length(rtree$tip)+1-4.5), x.lim=c(0, 0.35))
tips.plotordered <- rtree$tip[rtree$edge[rtree$edge[,2] <= length(rtree$tip),2]]

#par(mar=c(2,1,6,1))||

piRNA <- read.table("piRNATable.txt", sep=",")

rownames(piRNA)[match(names(species.names.fix), rownames(piRNA))] <- species.names.fix
stopifnot(all(rownames(piRNA) %in% tips.plotordered))

missing.species <- tips.plotordered[!tips.plotordered %in% rownames(piRNA)]
piRNA <- rbind(piRNA, as.data.frame(matrix(NA, ncol=ncol(piRNA), nrow=length(missing.species), dimnames=list(missing.species, colnames(piRNA)))))

# piRNA <- piRNA[tips.plotordered,]
piRNA <- piRNA[rev(tips.plotordered),]

piRNA <- as.matrix(piRNA)
piRNA <- t( t(piRNA) / piRNA["D.melanogaster",])

#image(x=seq_along(colnames(piRNA)), z=t(piRNA), axes=FALSE, xlab="", ylab="", col=hcl.colors(1024, "Purples2", rev = TRUE))
#axis(3, at=seq_along(colnames(piRNA)), tick=FALSE, label=colnames(piRNA), las=2)
# MIDPOINT ROOTING - to fix
# SpeciesTree <- midpoint.root(SpeciesTree)

piRNA <- as.data.frame(piRNA) %>%
  #Select column whose max value is greater than equal to 2 
  select_if(~max(., na.rm = TRUE) < 1.2)
piRNA <- as.matrix(piRNA)
# Also drop transib and transib3
piRNA <- piRNA[, colnames(piRNA) != "Transib1"]
piRNA <- piRNA[, colnames(piRNA) != "TRANSIB3"]
piRNA <- piRNA[, colnames(piRNA) != "HELITRON1"]
# Make > 1 to 1
piRNA[piRNA > 1] <- 1.0

# pdf("./AllSpeciesPlot.pdf", width=20, height=20)
# phylo.heatmap(rtree, piRNA, color=brewer.pal(n = 9, name = "YlOrRd"))
# heatmapRedone(rtree, piRNA, color=brewer.pal(n = 9, name = "Purples"))
# dev.off()
