#!/bin/bash

### Siddharth S. Tomar
### Aurélie Hua-Van
### Arnaud Le Rouzic


################ CHANGE PARAMETER HERE #############################
Project=~/CrossRegulation
method=gff # for RM analysis (with -gff option), otherwise choose blastn
#sra=SRR11846566
#sra=SRR5687217 # the file and path must be in the form $Project/Data/$sra.trimmed.collapse.fq
#sra=SRR25922470
sra=SRR14569563
#mismatch=0
mismatch=0
flavor=All
TEs=Data/RepBaseDmel.fa
Species=''
echo $1 $2 $3 $4 $5 $6 $7

if [ -n "$1" ]; then
    Project=$1
fi
if [ -n "$2" ]; then
    method=$2
fi
if [ -n "$3" ]; then
    sra=$3
fi
if [ -n "$4" ]; then
    mismatch=$4
fi
if [ -n "$5" ]; then
    flavor=$5
fi
if [ -n "$6" ]; then
    TEs=$6
fi
if [ -n "$7" ]; then
    Species=$7
fi

# flavor=Cp  do only the Cp count
# flavor=keep alternate=keep , for not redoing the first part which is independent of the mapping,
# flavor=All redo everything
cd $Project

echo $Project $method $sra $mismatch $flavor $TEs

######################## NEED ###################
# Genomes (Folder with genomes)
# RepeatMasker (Folder with RepeatMasker output (if already done)
# Data (with the SRA dataset)

######################### Clean before relaunching ($flavor=All)
if [ $flavor == All ];then
    rm -r AllSpecies/StatsSeq
    rm -r AllSpecies/Bed
    rm -r AllSpecies/Fasta
    rm -r AllSpecies/Bowtie
    rm -r AllSpecies/Bowtie.$sra.$mismatch
    rm -r AllSpecies/Stats
    rm -r AllSpecies/Index
    rm -r AllSpecies/Tree
    #list of TEs
    TElist=$(cat $TEs | grep '>' |cut -f 1 | cut -c2- |uniq)
    for te in $TElist
       do
         TEfile=TEcons/$te.cons.fa
        if [ ! -s ${TEfile} ];then
            grep -A1 -P '(?<![\w-])'$te'(?![\w-])'  $TEs > ${TEfile}
        fi
       done
fi

if [ $flavor == keep ];then
    rm -r AllSpecies/Bowtie.$sra.$mismatch
fi
################ For each genome
# will run the RM extraction of sequence (if $ flavor All), and the mapping (All and keep)
for file in $(ls Genomes/$Species*.fasta); do
    SpeciesName=$(basename "${file%.*}")
    echo ""
    echo $SpeciesName
    echo ""
    # Add more commands as needed
    bash Scripts/AnalyseOneSpecies.sh $Project $method $sra $mismatch $flavor $TEs $SpeciesName.fasta
done
cat AllSpecies/StatsSeq/*.CpNb.txt > AllSpecies/StatsSeq/AllSpecies.CpNb.txt
