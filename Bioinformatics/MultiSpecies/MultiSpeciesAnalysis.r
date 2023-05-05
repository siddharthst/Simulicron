# Load the required libs
suppressPackageStartupMessages({
require(tidyverse)
require(ape)
require(phytools)
require(Biostrings)
library(seqinr)
library(DescTools)}
)

# Make sure blast, seqtk, bedtools, bedops and mafft are installed and accessible via path
# Set additional parameters below
genomes        <- "~/Documents/Projects/Droso101/"
TEdbOrig       <- "./TEDB.fa"
SmallRNAReads  <- "~/Documents/Projects/SequencingLib/SRR14569563/reads_trimmed.fq"

# Distance threshold for merging TE fragments
TEfragmentDistance <- 250

# Species tree
speciesTree <- "./BigGene_.tre"

# Minimum consensus length
minimumConLength <- 1500

# Minumum number of piRNA reads in DM required to keep the TE family
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
dir.create(file.path(getwd(), "Sequences"), showWarnings = FALSE)
dir.create(file.path(getwd(), "Alignments"), showWarnings = FALSE)

# Threshold for minimum length
minTElen <- 0.7

# Generate list of TEs 
TEsToAnalyze <- names(read.fasta(file = TEdbFile, seqtype = "DNA"))

# Create a list of all the species present
listSpecies <- grep(".fasta$", list.files(path=genomes), value=TRUE) %>% str_replace(".fasta", "")

# Load consensus sequence
ConsensusSequence <- readDNAStringSet(TEdbFile)

# Create directories for the TEs
for (TE in TEsToAnalyze){
    dir.create(paste("./Results/",TE, sep=""), showWarnings = FALSE)
}

# Create dataframe to store the alignment information
ResultDataFrame <- data.frame(matrix(nrow = length(listSpecies), ncol = length(TEsToAnalyze)))
colnames(ResultDataFrame) <- TEsToAnalyze
rownames(ResultDataFrame) <- listSpecies
ResultDataFrame[is.na(ResultDataFrame)] <- 0

# Lists for storing stats
speciesToRemove <- c()

# Record for consensus sequences
ConsensusSequence <- readDNAStringSet(TEdbFile)

for(species in listSpecies){
    # Specify locations
    genomeLocation  <- paste(genomes, species, ".fasta", sep="")
    indexLocation   <- paste("./genomeDB/", species, sep="")
    NCBIResults     <- paste("./NCBIBLResults/", species, sep="")
    bedFileLocation <- paste("./BEDFiles/", species, ".bed", sep="")
    fastaLocation   <- paste("./Sequences/", species, ".fa", sep="")
    
    # First build a blast database
    cmd = paste("makeblastdb -in ", genomeLocation, " -dbtype nucl -out ", indexLocation, sep="")
    
    system(cmd)
    
    # For each species, identify the TE insertions and the amount of piRNA for each family
    # Run blast 
    cmd = paste("blastn -db ", indexLocation, " -query ", TEdbFile, " -outfmt 6 -evalue 20 -out ./NCBIBLResults/" ,species, sep="")
    system(cmd)
    
    # Check if there are any empty files, if yes, then do not process that species
    fileSize = file.info(NCBIResults)$size
    if(fileSize == 0){
        speciesToRemove <- c(speciesToRemove, species)
        next
    }
    
    # If not, then isolate the TEs, and append the name of the species and #
    NCBIResultFile    <- read.csv(NCBIResults, header=FALSE, sep="\t")
    NCBIResultFile$V1 <- paste(species, NCBIResultFile$V1, sep="_")
    write.table(NCBIResultFile, file=NCBIResults, quote=FALSE, sep='\t', row.names=FALSE, col.names=FALSE)
    cmd = paste("bash blast2bed.sh ", NCBIResults, sep="")
    system(cmd)
    
    # Work on the bed file
    cmd = paste("sort-bed ", NCBIResults, ".bed > ", bedFileLocation, sep="")
    system(cmd)
    unlink(paste(NCBIResults, ".bed", sep=""))
    
    # Merge the bed files
    cmd = paste("bedtools merge -s -d ", TEfragmentDistance, " -c 4,5,6 -o distinct,mean,distinct -i ", bedFileLocation, " > ", bedFileLocation, ".merged.bed", sep="")
    system(cmd)
    
    # Add unique number to each TE copy
    MergedResultFile      <- read.csv(paste(bedFileLocation, ".merged.bed", sep=""), header=FALSE, sep="\t")
    MergedResultFile$V4   <- paste(MergedResultFile$V4,seq.int(nrow(MergedResultFile)), sep="_")
    write.table(MergedResultFile, file=paste(bedFileLocation, ".merged.bed", sep=""), quote=FALSE, sep='\t', row.names=FALSE, col.names=FALSE)
    
    # Extract the sequences
    cmd = paste("bedtools getfasta -nameOnly -s -fi ", genomeLocation, " -bed ", bedFileLocation, ".merged.bed > ", fastaLocation, sep="")
    system(cmd)
    
    # Analyse the piRNA content
    for(TE in TEsToAnalyze){
        # Temporary files
        tempFasta = tempfile(pattern = TE, tmpdir = "./temp")
        tempMSA   = tempfile(pattern = TE, tmpdir = "./temp")
        tempIndex = tempfile(pattern = TE, tmpdir = "./temp")
        
        # Output names
        outputNames    <- paste("./Results/", TE,  "/", species, "_", TE, ".msa",sep="")
        bestTELocation <- paste("./Results/", TE,  "/", species, "_", TE, ".best.fa",sep="")
        
        # Isolate records
        cmd = paste("grep -A1 ", TE, " ", fastaLocation, " > ./temp/", TE, species, ".fa", sep="")
        system(cmd)
        fileSize = file.info(paste("./temp/", TE, species, ".fa", sep=""))$size
        
        # If no TE copies are present for the family
        if(fileSize == 0){
            ResultDataFrame[species, TE] = NA
            next
        }
        
        ConsensusLength <- nchar(toString(ConsensusSequence[TE]))
        minReqLength <- minTElen * ConsensusLength
        TempfastaLocation <- paste("./temp/", TE, species, ".fa", sep="")
        FastaSubset <- readDNAStringSet(TempfastaLocation)
        
        # Perform the MSA for all identified copies
        MSAset <- c(ConsensusSequence[TE])
        for (i in 1:length(FastaSubset)){
            # Write the object into fasta file
            writeXStringSet(c(ConsensusSequence[TE],FastaSubset[i]), tempFasta, append=FALSE, compress=FALSE, format="fasta")
            # Perform MSA using mafft
            cmd = paste("mafft --auto --preservecase ",  tempFasta," > ", tempMSA, sep="")
            system(cmd)
            # Read the results into another DNAStringSet
            MSA <- readDNAStringSet(tempMSA)
            unlink(tempFasta)
            unlink(tempMSA)
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
                    MSAset <- append(MSAset, MSA[2])
                }
            }
            # Check if there are any full length copies left
            if (length(MSAset) > 2)
                {
                MSAset <- MSAset[2:length(MSAset)]
                # Save the alignment and align reads to it
                writeXStringSet(MSAset, outputNames, append=FALSE, compress=FALSE, format="fasta")
                cmd=paste("bowtie-build -f ", outputNames, " ", tempIndex, sep="")
                system(cmd)
                pathToAlignment    = paste("./temp/", species, "_", TE, ".sam", sep="")
                pathToBAM          = paste("./temp/", species, "_", TE, ".bam", sep="")
                pathToAlignmentBED = paste("./temp/", species, "_", TE, ".bed", sep="")
                cmd = paste("bowtie -a -v 0 -p 8 --sam --no-unal -q -x ", tempIndex, " ", SmallRNAReads, " > ", pathToAlignment, sep="")
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
                # Check for any valid alignments
                if (file.exists(pathToAlignmentBED)){
                    fileSize = file.info(pathToAlignmentBED)$size
                    if (fileSize > 0){
                        # Load the alignment and find the the copy with most hits
                        AlignmentDataFrame <- as.data.frame(read.table(pathToAlignmentBED, header = FALSE, sep="\t",stringsAsFactors=FALSE, quote=""))
                        bestTE <- Mode(AlignmentDataFrame$V1)[1]
                        # bestTE <- gsub("[()+-]", "", bestTE)
                        ResultDataFrame[species, TE] = length(which(AlignmentDataFrame$V1==bestTE))
                        # Isolate the best TE for each species for tree construction
                        cmd = paste("grep -A1 '", bestTE, "' ", outputNames, " > ", bestTELocation, sep="")
                        system(cmd)
                    }
                }
                else{
                    ResultDataFrame[species, TE] = 0
                }
            
            }
        }
    cmd = "rm -rf ./temp/*"
    system(cmd)
}

save.image(file = "1.RData")

write.table(ResultDataFrame, file='./Results/Stats.tsv', quote=FALSE, sep='\t')


