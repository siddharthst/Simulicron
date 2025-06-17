#!/bin/bash

### Siddharth S. Tomar
### Aurélie Hua-Van
### Arnaud Le Rouzic

#### Please ensure that you have installed:
# blast+
# bedtools 2.30
# python 3
# mafft
# FastTree
# bowtie
# R and R packages ape and phytools

### bash Scripts/AnalyseDmel.sh erase iso1 r3 SRR11846566 3

#bash Scripts/AnalyseDmel.sh $Project $strain  $sra $pi $mismatch $flavor $TEs 

# defaults parameters
Project=~/CrossRegulation
strain=iso1
pi=r0 # Pi annot
sra=SRR11846566
#sra=SRR5687217 # the file and path must be in the form $Project/Data/$sra.trimmed.collapse.fq
#sra=SRR25922470
#sra=SRR14569563
#mismatch=0
mismatch=3
flavor=All
TEs=Data/RepBaseDmel.fa

# parameters new

if [ -n "$1" ]; then
    Project=$1
fi
if [ -n "$2" ]; then
    strain=$2
fi
if [ -n "$3" ]; then
    sra=$3
fi
if [ -n "$4" ]; then
    pi=$4
fi
if [ -n "$5" ]; then
    mismatch=$5
fi
if [ -n "$6" ]; then
    flavor=$6
fi
if [ -n "$7" ]; then
    TEs=$7
fi

echo $1 $2 $3 $4 $5 $6 $7

#########################################################################################
# Initialize

cd $Project

#### run

### bash Scripts/AnalyseDmel.sh erase iso1 r3 SRR11846566 3
### bash Scripts/AnalyseDmel.sh keep iso1 r1 SRR11846566 3
### bash Scripts/AnalyseDmel.sh no iso1
### bash Scripts/AnalyseDmel.sh erase OreR
### bash Scripts/AnalyseDmel.sh keep OreR
### bash Scripts/AnalyseDmel.sh no OreR
### bash Scripts/AnalyseDmel.sh erase CanS
### bash Scripts/AnalyseDmel.sh keep CanS
### bash Scripts/AnalyseDmel.sh no CanS
#########################################################################################
# Argument 2

#### Choice of the pi annotation
#### if you change and do not want to erase, choose "no"
#### Please comment the ones you do not want, by default then r0
#pi=r1
#pi=r2
#pi=r0
echo $Project $strain  $sra $pi $mismatch $flavor $TEs # default = 0


#### Choice of the genome
#if [[ $2 == iso1 ]]; then
    #strain=iso1
    #genomeVersion=Dmel.$strain # genome name
    #piC=piCluster.dmelIso1.$pi.bed
#fi

#if [[ $2 == OreR ]]; then
    #strain=OreR
    #genomeVersion=Dmel.$strain # genome name
    #piC=piCluster.dmelOreR.$pi.bed
#fi

#if [[ $2 == CanS ]]; then
    #strain=CanS
    #genomeVersion=Dmel.$strain # genome name
    #piC=piCluster.dmelCanS.$pi.bed
#fi
genomeVersion=Dmel.$strain # genome name
genome=Data/$genomeVersion.fasta # genome fasta file

piC=piCluster.dmel.$strain.$pi.bed



combi=$strain.$pi
echo $strain
echo $combi
echo $genome
echo $piC


#########################################################################################
# Argument 1
# keep or erase
# keep do the pi mapping only, keeping the TE mapping
# erase  remove the folders
# no option should behave keep the folders but redo everything (some old stuff may persit and create errors
if [[ $flavor == keep ]]; then
    echo "keep"
fi
if [[ $flavor == erase ]]; then
    
    rm -r ResultsDmel.$strain/Tree
    rm -r ResultsDmel.$strain/AllTEs
    rm -r ResultsDmel.$strain/Bed
    rm -r ResultsDmel.$strain/Stats
    rm -r ResultsDmel.$strain/Figures
    rm -r ResultsDmel.$strain/Bowtie
    rm -r Index
    #exit 0

fi

if [[ $flavor == mapPi ]]; then
    
    rm   ResultsDmel.$strain/AllTEs/$genomeVersion.$pi.$sra.$distance.$mismatch.Cp.summary.txt
    #rm -r ResultsDmel.$strain/Bowtie
fi

mkdir -p ResultsDmel.$strain
mkdir -p ResultsDmel.$strain/AllTEs
mkdir -p ResultsDmel.$strain/Bowtie
mkdir -p ResultsDmel.$strain/Bed
mkdir -p ResultsDmel.$strain/Tree
mkdir -p ResultsDmel.$strain/Figures
mkdir -p ResultsDmel.$strain/Stats
mkdir -p Index
#########################################################################################
# The variables . Don't forget to start defining them each time you reconnect to the serveur

# Variable to test
#sra=SRR14569563 #piRNA dataset default Gebert 2021
#sra=SRR25922470 # piRNA from Andy Clarck
#sra=SRR5687217
#sra=SRR11846566 #piRNA dataset iso1 from Kofler

distance=250 	# distance for merging

threshold1=0.7     # length threshold to keep copies not in piCluster
N=50           # Number of cpu to use





#List of TEs
TElist=$(cat $TEs | grep '>' |cut -f 1 | cut -c2- |uniq)
#echo $TElist
filename=$(basename -- "$TEs")

TEfile=${filename%.*}
echo $TEfile

smallRNAreads=Data/$sra.trimmed.collapse.fq # small RNA dataset

# Total nb of reads in the not collapsed piRNA dataset
nbreads=$(echo $(cat Data/$sra.trimmed.fq|wc -l)/4|bc)

# We expect that the piRNA reads (decollapsed) for one TE is above 1/10000 of the total nb of reads
# Minimum piNb per pi copy - should depend of the size of the dataset
# To be uses when creating the matrix of sahre reads
threshold2=$(echo $nbreads/10000 |bc )   


var=$pi.$sra.$distance.$mismatch


#TElist=Gypsy1_DM
statfile=ResultsDmel.$strain/Stats/Stat.$genomeVersion.$var.txt


echo $(date) > $statfile

#########################################################################################
# Check
echo $(which blastn) >> $statfile
echo $(bedtools --version) >> $statfile
echo $(python3 --version) >> $statfile
echo $(which mafft) >> $statfile
echo $(which FastTree ) >> $statfile
echo $(which bowtie) >> $statfile
echo -e '\n' >> $statfile

echo -e $distance >> $statfile
echo -e $sra >> $statfile
echo -e $genomeVersion >> $statfile
echo -e '\n' >> $statfile

##################################################
################# FUNCTIONS
##################################################

##################################################
############ build the genome db, blast
############ Run the global analysis : transform to bed, tag copies
##################################################
Blastgenome(){
    echo '>>> blastn TEs on genome'

    makeblastdb -in $genome -dbtype nucl -out Index/$genomeVersion

    blastn -db Index/$genomeVersion -query $TEs -outfmt "6 std qlen slen" -evalue 20 -num_threads $N -out ResultsDmel.$strain/AllTEs/$genomeVersion.$TEfile.blastn

    #Check you have permission
    sort -k2,2 -k9n,9 -k10n,10 ResultsDmel.$strain/AllTEs/$genomeVersion.$TEfile.blastn -o ResultsDmel.$strain/AllTEs/$genomeVersion.$TEfile.sorted.blastn

    Scripts/blast2bed.sh ResultsDmel.$strain/AllTEs/$genomeVersion.$TEfile.blastn
    # output ResultsDmel.$strain/AllTEs/$genomeVersion.$TEfile.blastn.bed
    
    # Merge the result for altogether analysis
    bedtools sort -i ResultsDmel.$strain/AllTEs/$genomeVersion.$TEfile.blastn.bed  | mergeBed -s -d $distance -c 4,5,6,7,8,9,9 -o distinct,mean,distinct,min,max,min,max > ResultsDmel.$strain/AllTEs/$genomeVersion.$TEfile.$distance.merged.bed

    sort -k1,1 -k2n,2 -i ResultsDmel.$strain/AllTEs/$genomeVersion.$TEfile.$distance.merged.bed -o ResultsDmel.$strain/AllTEs/$genomeVersion.$TEfile.$distance.merged.sorted.bed

    # Mark the merged.bed file for merged hits (when several names in col 4)
    # Output an overlap file used in IdentifyPiCopies. Overlapping copies are tagged "_r"
    python3 Scripts/TagRed.py ResultsDmel.$strain/AllTEs/$genomeVersion.$TEfile.$distance.merged.bed

}
##################################################
#### Run TE by TE analysis
#### transform blastn to bed, sort, merge
#### Extract sequence
##################################################
ExtractNAlign(){

    
    echo '>>> Extract sequences for '  $te
    # Save the consensus sequence of the TE (seq must be on one line) Done before
    #grep -A1 -w $te $TEs > ResultsDmel.$strain/Bed/$te.cons.seq
    #grep -A1 -P '(?<![\w-])'$te'(?![\w-])'  $TEs > ResultsDmel.$strain/Bed/$te.cons.seq
    
    # Filter by TE and save the hits from the TE
    cat ResultsDmel.$strain/AllTEs/$genomeVersion.$TEfile.blastn | grep -P '(?<![\w-])'$te'(?![\w-])' > ResultsDmel.$strain/Bed/$te.blastn
    # Transform blast into bed file
    Scripts/blast2bed.sh ResultsDmel.$strain/Bed/$te.blastn
    # Sort and merged the bed file (TE by TE)
    bedtools sort -i ResultsDmel.$strain/Bed/$te.blastn.bed | bedtools merge -s -d $distance -c 4,5,6,7,8,9,10 -o distinct,mean,distinct,min,max,collapse,collapse > ResultsDmel.$strain/Bed/$te.$distance.merged.bed
    # Get the sequences corresponding to the merged file
    bedtools getfasta -s -fi Data/$genomeVersion.fasta -bed ResultsDmel.$strain/Bed/$te.$distance.merged.bed > ResultsDmel.$strain/Bed/$te.$distance.fa
}

##################################################
#### Run TE by TE analysis
#### Identify copies into pCluster (depend on the piAnnotation $pi)
#### Filter for size and align
######################################
IdentifyPiCopies(){
    # Initialize TE stat file
    telogfile=ResultsDmel.$strain/Stats/$te.$distance.stats.txt
    echo $(date) > $telogfile
    # Detect TE in piCluster, and merged TE
    bedtools intersect -a ResultsDmel.$strain/Bed/$te.$distance.merged.bed -b Data/$piC > ResultsDmel.$strain/Bed/$te.$distance.$pi.intersect.bed
    bedtools intersect -a ResultsDmel.$strain/Bed/$te.$distance.merged.bed -b ResultsDmel.$strain/AllTEs/$genomeVersion.$TEfile.$distance.overlap.bed > ResultsDmel.$strain/Bed/$te.$distance.intersect.overlap.bed

    # Change name and filter for size (0.7 - threshold1) for the euchromatic copies only
    # Produce a renamed.filtered file (pi and _r)
    ######## Concat result summary
    python3 Scripts/ChangeName.py  ResultsDmel.$strain/Bed/$te.$distance $pi $threshold1 >> ResultsDmel.$strain/AllTEs/$genomeVersion.$TEfile.$pi.$distance.filter1.summary.txt
    
    # Align sequences by removing all nt insertions compared to consensus
    FILENAME=ResultsDmel.$strain/Bed/$te.$distance.$pi.renamed.filtered.fa
    if [ ! -s "${FILENAME}" ]; then
        rm $FILENAME
        exit 0
    else
        mafft --quiet --reorder --keeplength --addfragments  ResultsDmel.$strain/Bed/$te.$distance.$pi.renamed.filtered.fa ResultsDmel.$strain/Bed/$te.cons.seq > ResultsDmel.$strain/Tree/$te.$distance.$pi.cons.afa 2>> $telogfile
    fi
    }
    
##################################################
#### Run TE by TE analysis
#### Map reads on te, te by te (depend on $sra)
#### Count reads
#### Construct the tree
######################################    
MapNAnalyse(){
        # Build bowtie index (on full sequence, before removing insertion)
    bowtie-build --quiet -f ResultsDmel.$strain/Bed/$te.$distance.$pi.renamed.filtered.fa Index/$te.$distance.$pi.renamed.filtered
    # Map the small RNA
    bowtie -a -v $mismatch -p 8 --sam --no-unal -q -x Index/$te.$distance.$pi.renamed.filtered $smallRNAreads 2>> ResultsDmel.$strain/Stats/$te.$pi.$distance.stats.txt | samtools sort | bamToBed > ResultsDmel.$strain/Bowtie/$genomeVersion.$te.$var.bed
    
    
    # Count reads # save matrix and summaries in Bowtie folder
    ######## Concat result summary for TE
    ######## Concat reads summary for TE

    totCpNb=$(cat ResultsDmel.$strain/Bed/$te.$distance.merged.bed| wc -l )
    ## Analyse
    python3 Scripts/AnalyseBed.py ResultsDmel.$strain/Bowtie/$genomeVersion.$te.$var.bed $te ResultsDmel.$strain/Bed/$te.$distance.$pi.summary.txt $totCpNb $threshold2 >> ResultsDmel.$strain/AllTEs/$genomeVersion.$var.TE.summary.txt

    ######## Concat result summary for Copies and piRNA (all TEs)
    cat ResultsDmel.$strain/Bowtie/$genomeVersion.$te.$var.bed.bowtie.summary.txt >> ResultsDmel.$strain/AllTEs/$genomeVersion.$var.Cp.summary.txt
    cat ResultsDmel.$strain/Bowtie/$genomeVersion.$te.$var.bed.bowtie.histo.txt >> ResultsDmel.$strain/AllTEs/$genomeVersion.$TEfile.$var.histo.txt
    
    # To do only if the TE has been kept (based on  minimum number of sequences, files created )
    FILE=ResultsDmel.$strain/Bowtie/$genomeVersion.$te.$var.bed.kept.txt
    FILEDISCARD=ResultsDmel.$strain/Bowtie/$genomeVersion.$te.$var.bed.discarded.txt
    if [ -f "$FILE" ]; then
        echo $te 'Analyse and Construct tree'
        # Add content of file in TE stat file
        cat $FILE >> ResultsDmel.$strain/Stats/$te.$pi.$distance.stats.txt
        # Remove discarded copies (< 200 pi reads) and cons from alignment
        # Create a TE.afa file in Tree folder
        python3 Scripts/DropCons.py ResultsDmel.$strain/Tree/$te.$distance.$pi.cons.afa $FILE $te
        echo 'Kept for analysis' >> ResultsDmel.$strain/Stats/$te.$distance.stats.txt
        # Construct tree, save stderr in stat file
        FastTree -quiet -pseudo -gtr -nt ResultsDmel.$strain/Tree/$te.$var.afa 2>> ResultsDmel.$strain/Stats/$te.$distance.stats.txt 1> ResultsDmel.$strain/Tree/$te.$var.tre
        else
        echo 'Discarded from analysis' >> ResultsDmel.$strain/Stats/$te.$distance.stats.txt
        python3 Scripts/DropCons.py ResultsDmel.$strain/Tree/$te.$distance.$pi.cons.afa $FILEDISCARD $te
        FastTree -quiet -pseudo -gtr -nt ResultsDmel.$strain/Tree/$te.$var.afa 2>> ResultsDmel.$strain/Stats/$te.$distance.stats.txt 1> ResultsDmel.$strain/Tree/$te.$var.Discarded.tre
    fi
}



if [[ $flavor == Cp ]]
    then
        Blastgenome
elif [[ $flavor == CpPi ]]
    then
        Blastgenome
elif [[ $flavor == erase ]]
    then
        Blastgenome
else
	echo "Skipping blast"    
    #ExtractNAlign
fi
##########################################
##########################################
################ Saving stat ################
echo -e $distance 'nt for merging' >> $statfile
echo -e 'Genome' $(LL_ALL=C grep '>' $genome|wc -l) 'sequences' >> $statfile
echo -e 'TEfile' $(LL_ALL=C grep '>' $TEs|wc -l) 'sequences' >> $statfile
echo -e 'SRA' $(echo $(LL_ALL=C cat $smallRNAreads|wc -l) /4|bc) 'sequences' >> $statfile
echo -e '\n' >> $statfile
echo -e 'Total blastn hits:' $(wc -l ResultsDmel.$strain/AllTEs/$genomeVersion.$TEfile.blastn.bed) >> $statfile
echo -e 'Merged hits: ' $(wc -l ResultsDmel.$strain/AllTEs/$genomeVersion.$TEfile.$distance.merged.bed) >> $statfile
wc -l ResultsDmel.$strain/AllTEs/$genomeVersion.$TEfile.$distance.merged.bed
# 33269 250
# 13430  50

## TODO
#echo $(date) > Stats/Stat.txt
##########################################################
######### FOR EACH TE
##########################################################

## save summary in AllTEs folder
## for each TE, after filter 1 (python ChangeName)
echo -e 'TE\tTotCpNb\tFiltCpNb\tEucCpNb\tPiCpNb\tExcNonPi' > ResultsDmel.$strain/AllTEs/$genomeVersion.$TEfile.$distance.filter1.summary.txt

## save the bowtie reads
#> ResultsDmel.$strain/AllTEs/0_AllTE.individual.reads.bed

## for each TE, after read mapping and analysis (python AnalyseBed)
echo -e 'TE\tTotCpNb\tFiltCpNb\tKeptCpNb\tPiCpNb\tEucCpNb\tExcPi\tNbCollapsedReads\tTotalNbReads\tFreqMinus\tNbReadsKept\tReadsInNonPi\tReadsInPi\tShared\tMax\tmaxpi\tmaxte\tStatus' > ResultsDmel.$strain/AllTEs/$genomeVersion.$var.TE.summary.txt

## for each Copy, after read mapping and analysis
echo -e 'Chr\tStart\tEnd\tTE\tzero\tstrand\tqstart\tqend\tevalue\tdiv\tpercentStart\tlength\tConsLength\tCoordName\tCopyName\tMeanDiv\tstatus\tNbPi\tFinalStatus' > ResultsDmel.$strain/AllTEs/$genomeVersion.$var.Cp.summary.txt
##########################################################################
############## ANALYSE ###################################################
for te in $TElist; do
    (
       # .. do your stuff here
        grep -A1 -P '(?<![\w-])'$te'(?![\w-])'  $TEs > ResultsDmel.$strain/Bed/$te.cons.seq

        if [[ $flavor == keep ]]
            then 
            echo "Skipping sequence extraction, Identifying pi and align"
            IdentifyPiCopies $te
        elif [[ $flavor == mapPi ]]
            then echo "Skipping sequence extraction and pi Identification"
        else
            echo $te 'Extract, align,identify Pi Copies'
            ExtractNAlign $te
            IdentifyPiCopies $te
        fi
        
        FILENAME=ResultsDmel.$strain/Bed/$te.$distance.$pi.renamed.filtered.fa
        if [ ! -s "${FILENAME}" ] && [ -e "${FILENAME}" ]; then
            rm $FILENAME
            exit 0
        elif [ ! -e "${FILENAME}" ]; then
            exit 0
        elif [[ $flavor == Cp ]];then
        	echo "Skipping Mapping"
        elif  [[ $flavor == CpPi ]];then
        	echo "Skipping Mapping"
        else
            MapNAnalyse $te
        fi
    ) &

   # allow to execute up to $N jobs in parallel
    if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
       # now there are $N jobs already running, so wait here for any job
       # to be finished so there is a place to start next one.
        wait -n
    fi

done

# no more jobs to be started but wait for pending jobs
# (all need to be finished)
wait

echo "All done"


######################################################################
# Concat and sort all Cp information
cat ResultsDmel.$strain/Bed/*.$distance.$pi.summary.txt |sort -k1,1 -k2n,2 > ResultsDmel.$strain/AllTEs/$genomeVersion.$pi.$sra.$distance.Cp.filter1.bed

######################################################################
##### Mapping statistics 
#> ResultsDmel.$strain/AllTEs/$genomeVersion.$sra.$mismatch.bowtie.uniq.count

# Create file
> ResultsDmel.$strain/AllTEs/$genomeVersion.$var.bowtie.AllReads.uniq

Analyse(){
#    echo $te 'Compute Stats'
#    # Count and save to Stat file
#    iCpNb=$(grep '>' ResultsDmel.$strain/Bed/$te.$distance.fa |wc -l| cut -d ' ' -f 1)
#    iNtNb=$(cat ResultsDmel.$strain/Bed/$te.$distance.fa | paste - - | cut -f 2 | tr -d '\n' | wc -c)
#    fCpNb=$(grep '>' ResultsDmel.$strain/Bed/$te.$distance.renamed.filtered.fa |wc -l | cut -d ' ' -f 1)
#    fNtNb=$(cat ResultsDmel.$strain/Bed/$te.$distance.renamed.filtered.fa | paste - - | cut -f 2 | tr -d '\n' | wc #-c)
#    iPi=$(wc -l ResultsDmel.$strain/Bed/$te.$distance.intersect.pi.bed| cut -d ' ' -f 1)
#    fPi=$(grep 'piRNA' ResultsDmel.$strain/Bed/$te.$distance.renamed.filtered.fa|wc -l| cut -d ' ' -f 1)
#    iOL=$(wc -l ResultsDmel.$strain/Bed/$te.$distance.intersect.overlap.bed| cut -d ' ' -f 1)
#    fOL=$(grep '*' ResultsDmel.$strain/Bed/$te.$distance.renamed.filtered.fa|wc -l | cut -d ' ' -f 1)
#    echo -e $te $iCpNb $iNtNb $iPi $iOL >> $statfile
#    echo -e $te $fCpNb $fNtNb $fPi $fOL >> $statfile
#    echo -e "" "CpNb" "NtNb" "Pi" "OL" >> ResultsDmel.$strain/Stats/$te.$distance.stats.txt
#    echo -e "Before merging :" $iCpNb $iNtNb $iPi $iOL >> ResultsDmel.$strain/Stats/$te.$distance.stats.txt
#    echo -e "After merging :" $fCpNb $fNtNb $fPi $fOL >> ResultsDmel.$strain/Stats/$te.$distance.stats.txt
#    #### write stats
#    echo -e "Initial number of hits in genome" $(wc -l ResultsDmel.$strain/Bed/$te.blastn.bed) >> #ResultsDmel.$strain/Stats/$te.$distance.stats.txt
#    echo -e "Number of hits in genome after merging" '\t' $(wc -l #ResultsDmel.$strain/Bed/$te.$distance.merged.bed) >> ResultsDmel.$strain/Stats/$te.$distance.stats.txt
#    echo -e "Number of hits in genome after filtering" '\t' $(wc -l #ResultsDmel.$strain/Bed/$te.$distance.summary.txt) >> ResultsDmel.$strain/Stats/$te.$distance.stats.txt
#    echo -e "including piCluster copies:" '\t' $(grep 'piRNA' ResultsDmel.$strain/Bed/$te.$distance.summary.txt | #wc -l) >> ResultsDmel.$strain/Stats/$te.$distance.stats.txt
    

    # Add bowtie bed results to AllTEs file
    #cat ResultsDmel.$strain/Bowtie/$genomeVersion.$te.$sra.$mismatch.bed >> ResultsDmel.$strain/AllTEs/0_AllTE.individual.reads.bed
    
    # list of unique reads (still collapsed). One TE per line.
    UniqTE=$(cat ResultsDmel.$strain/Bowtie/$genomeVersion.$te.$var.bed| cut -f 4| sort|uniq )
    # add the number of unique reads per TE to file uniq (and nb of character)
#    echo -e $te $(echo $UniqTE|wc -wc) >> ResultsDmel.$strain/AllTEs/$genomeVersion.$sra.$mismatch.bowtie.uniq.count
    # add the unique reads per TE to file uniq
    echo -e $te $UniqTE >> ResultsDmel.$strain/AllTEs/$genomeVersion.$var.bowtie.AllReads.uniq
}

(
for te in $TElist
    do
        #((i=i%N)); ((i++==0)) && wait
        FILENAME=ResultsDmel.$strain/Bowtie/$genomeVersion.$te.$var.bed
        if [ -f "$FILENAME" ]; then
            Analyse $te
        fi
    done
)
wait

# Concatenate merged blastn bed files to AllTEs file
cat ResultsDmel.$strain/Bed/*.$distance.merged.bed > ResultsDmel.$strain/AllTEs/0_AllTE.individual.$distance.merged.bed

python3 Scripts/CountRead.py ResultsDmel.$strain/AllTEs/$genomeVersion.$var.bowtie.AllReads.uniq
echo 'Computing stats done'
# nb TE, useless
#cat ResultsDmel.$strain/AllTEs/0_AllTE.individual.reads.bed | cut -f 4| sort|uniq |wc -l

#########################################################################
###################### CHECKING redundancy ###############################

### by comparing annotations for all TE directly from the merged blast or from the merged TE separated blast files
GenomeCov(){
    # Calculate coverage nt from initial merged  file merged (no redundancy expected)
    bedtools genomecov  -i  ResultsDmel.$strain/AllTEs/$genomeVersion.$TEfile.$distance.merged.bed -g $genome.length  > ResultsDmel.$strain/AllTEs/0_AllTE.together.$distance.genomecov

    # Calculate coverage nt from TE merged bed file concatenated (redundancy expected)
    ### Concatenate merged TE annotation to estimate redundancy (see table1)

    bedtools sort -i ResultsDmel.$strain/AllTEs/0_AllTE.individual.$distance.merged.bed > ResultsDmel.$strain/AllTEs/0_AllTE.individual.$distance.sorted.bed

    bedtools genomecov  -i ResultsDmel.$strain/AllTEs/0_AllTE.individual.$distance.sorted.bed -g $genome.length > ResultsDmel.$strain/AllTEs/0_AllTE.individual.$distance.genomecov
}
if [[ $1 == keep ]]; then echo "Skipping GenomeCov"; else GenomeCov; fi
wait
# save time
echo $(date) >> $statfile

