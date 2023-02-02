# Make sure blast, seqtk, bedtools and mafft are installed and accessible via path
# set additional parameters below
genome       <- "~/Documents/Projects/Genomes/Drosophila_melanogaster.BDGP6.32.dna.toplevel.fa"
TEsToAnalyze <- c("IDM", "GYPSY", "COPIA", "BATUMI")
TEdbFile     <- "./TEDB.fa"
# Distance threshold for merging TE fragments
TEfragmentDistance <- 250
# Create dirs
dir.create(file.path(getwd(), "temp"), showWarnings = FALSE)
# Threshold for minimum length
minTElen <- 0.7

# Load the required libs
require(tidyverse)
require(ape)
require(phytools)
require(Biostrings)

# First build a blast database
cmd = paste("makeblastdb -in ", genome, " -dbtype nucl -out ./genomeDB", sep="")
system(cmd)

# Now run BLAST on the sequences
cmd = paste("blastn -db ./genomeDB -query ", TEdbFile, " -outfmt 6 -evalue 20 -out results.tab", sep="")
system(cmd)

# Convert the output to a BED file
# https://github.com/nterhoeven/blast2bed
cmd = "bash blast2bed.sh ./results.tab"
system(cmd)

# Working on each TE individually 
ConsensusSequence <- readDNAStringSet("./TEDB.fa")
for(TE in TEsToAnalyze){
    cmd = paste("grep \t", TE, "\t ./results.tab.bed | bedtools sort > ", TE, ".bed", sep="")
    system(cmd)
    # Merge TEs
    cmd = paste("bedtools merge -s -d 250 -c 4,5,6 -o distinct,mean,distinct -i ", TE, ".bed > ./", TE, ".merged.bed", sep="")
    system(cmd)
    # Extract the sequences
    cmd = paste("bedtools getfasta -s -fi ", genome, " -bed ", TE, ".merged.bed >", TE, ".fa", sep="")
    system(cmd)
    # Perform MSA
    ConsensusLength <- nchar(toString(ConsensusSequence[TE]))
    minReqLength <- minTElen * ConsensusLength
    fastaLocation = paste("./", TE, ".fa", sep="")
    FastaSubset <- readDNAStringSet(fastaLocation) 
    MSAset <- c(ConsensusSequence[TE])
    Alignmentset <- c(ConsensusSequence[TE])
    outputNames <- paste(TE, ".msa",sep="")
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
                MSAset <- append(MSAset, MSA[2])
                # Append the results to Alignment file by replacing - with N
                MSA[[2]] <- gsub("-", "N", MSA[[2]])
                Alignmentset <- append(Alignmentset, MSA[2])
        }
        # Remove the temp files
        unlink(tempFasta)
        unlink(tempMSA)
    }  
    MSAset <- MSAset[2:length(MSAset)]                                        
    writeXStringSet(MSAset, outputNames, append=FALSE, compress=FALSE, format="fasta")
    # Before generating trees, we need to rename the fasta headers as fasttree doesn't like :
    cmd = paste("seqtk rename ", TE, ".msa > ", TE, ".msa.fixed", sep="")
    system(cmd)
    # Generate the tree
    cmd = paste("fasttree -nt -gtr < ", TE, ".msa.fixed > ", TE, ".tree", sep="")
    system(cmd)
    # Load the tree and plot it
    tree <- read.tree(paste(TE, ".tree", sep=""))
    tree <- midpoint.root(tree)
    pdf(paste(TE, ".pdf", sep=""))
    plot(tree)
    add.scale.bar()
    dev.off()
}
