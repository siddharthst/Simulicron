#!/bin/bash

### Siddharth S. Tomar
### Aurélie Hua-Van
### Arnaud Le Rouzic

#### Please ensure that you have installed:
# blast+ or RepeatMasker
# bedtools
# python 3
# mafft
# FastTree2
# bowtie
# R and R packages ape, phytools, tidyverse, ggtree



#########################################################################################
# Initialize
Project=$1
#cd /mnt/35To/huavan/CrossRegulation
cd $Project

#### run
### bash Scripts/AnalyseAllSpecies.sh $Project Z.ornatus.fasta r1 SRR5687217
method=$2
#######################################################################################
#Choice of the species
species=$7
# Get read of fasta extension
SpeciesName="${species%.*}"
genome=Genomes/$SpeciesName.fasta # genome fasta file
echo $genome

# $5 flavor
	# flavor=Cp  do only the Cp count (not implemented)
	# flavor=keep, mapping small reads only, for not redoing the first part which is independent of the mapping,
	# flavor=All redo everything (except RM)
	# flavor=RM redo everything including RM or blastn  (not implemented, done if not keep and not All)

#########################################################################################
# The variables . Don't forget to start defining them each time you reconnect to the serveur

#cd /mnt/Pole/Genomes/huavan/CrossRegulation
cd $Project
# Variable to test
#sra=SRR14569563 #piRNA dataset default Gebert 2021
#sra=SRR25922470 # piRNA from Andy Clarck
#sra=SRR5687217
#sra=SRR11846566 #piRNA dataset iso1 from Kofler 2021
sra=$3
#mismatch=0      # mismatch allowed for mapping to copy
mismatch=$4
echo $sra $mismatch

distance=250     # distance for merging
threshold1=0.7     # length threshold to keep copies not in piCluster
N=10           # Number of cpu to use



TEfile=RepBaseDmel #TE dataset
TEs=Data/$TEfile.fa
TElist=$(cat $TEs | grep '>' |cut -f 1 | cut -c2- |uniq)



smallRNAreads=Data/$sra.trimmed.collapse.fq # small RNA dataset
#List of TEs
#echo $TElist
#TElist=Gypsy1_DM
var=$sra.$distance.$mismatch

################################
### Create adhoc folders
################################
mkdir -p AllSpecies
mkdir -p AllSpecies/StatsSeq # Stat for the copy number of each TE/Species
mkdir -p AllSpecies/Bed # merged and size filtered files per te (RM or BLAST)
mkdir -p AllSpecies/Fasta # aligned filtered files with consensus
mkdir -p AllSpecies/Bowtie # species all.filtered.fa (degap, no cons) for mapping
mkdir -p AllSpecies/Bowtie.$sra.$mismatch # mapping results per species
mkdir -p AllSpecies/Tree
mkdir -p AllSpecies/Figures
mkdir -p AllSpecies/Stats
mkdir -p TEcons
mkdir -p AllSpecies/Index
mkdir -p AllSpecies/Table
mkdir -p AllSpecies/Results
mkdir -p AllSpecies/D.melanogaster

####################################
############ STAT file ##############
statfile=AllSpecies/Stats/$SpeciesName.$var.Stats.txt
# save time
echo $(date) >> $statfile
# Checking saved to statfile
echo $(which blastn) >> $statfile
echo $(bedtools --version) >> $statfile
echo $(python3 --version) >> $statfile
echo $(which mafft) >> $statfile
echo $(which FastTree ) >> $statfile
echo $(which bowtie) >> $statfile
echo -e '\n' >> $statfile

echo -e $distance 'nt for merging' >> $statfile
echo -e $sra >> $statfile
echo -e $SpeciesName >> $statfile
echo -e '\n' >> $statfile

#############################################
################# FUNCTIONS
#############################################
BlastNAlign(){
    # Initialize TE stat file
    telogfile=AllSpecies/Stats/$te.$distance.stats.txt
    echo $(date) > $telogfile
    echo $te 'Extract sequences'
    # Save the consensus sequence of the TE (seq must be on one line)
    #grep -A1 -w $te $TEs > AllSpecies/Bed/$te.cons.seq
    grep -A1 -P '(?<![\w-])'$te'(?![\w-])'  $TEs > AllSpecies/Bed/$te.cons.seq
    # Filter and save the hits from the TE
    cat AllSpecies/AllTEs/$SpeciesName.$TEfile.blastn | grep -P '(?<![\w-])'$te'(?![\w-])' > AllSpecies/Bed/$te.blastn
    # Transform blast into bed file
    Scripts/blast2bed.sh AllSpecies/Bed/$te.blastn
    # Sort and merged the bed file
    bedtools sort -i AllSpecies/Bed/$te.blastn.bed | bedtools merge -s -d $distance -c 4,5,6,7,8,9,10 -o distinct,mean,distinct,min,max,collapse,collapse > AllSpecies/Bed/$te.$distance.merged.bed
    # Get the sequences corresponding to the merged file
    bedtools getfasta -s -fi Data/$SpeciesName.fasta -bed AllSpecies/Bed/$te.$distance.merged.bed > AllSpecies/Bed/$te.$distance.fa

    # Detect TE in piCluster, and merged TE
    bedtools intersect -a AllSpecies/Bed/$te.$distance.merged.bed -b Data/$piC > AllSpecies/Bed/$te.$distance.intersect.pi.bed
    bedtools intersect -a AllSpecies/Bed/$te.$distance.merged.bed -b AllSpecies/AllTEs/$SpeciesName.$TEfile.$distance.overlap.bed > AllSpecies/Bed/$te.$distance.intersect.overlap.bed

    # Change name and filter for size (0.7 - threshold1)
    # Produce a renamed.filtered file
    ######## Concat result summary
    python Scripts/ChangeName.py  AllSpecies/Bed/$te.$distance $threshold1 >> AllSpecies/AllTEs/$SpeciesName.$TEfile.$distance.filter1.summary.txt
    
    # Align sequences by removing all insertions compared to consensus
    FILENAME=AllSpecies/Bed/$te.$distance.renamed.filtered.fa
    if [ ! -s "${FILENAME}" ]; then
        rm $FILENAME
        exit 0
    else
        mafft --quiet --reorder --keeplength --addfragments  AllSpecies/Bed/$te.$distance.renamed.filtered.fa AllSpecies/Bed/$te.cons.seq > AllSpecies/Tree/$te.$distance.cons.afa 2>> $telogfile
    fi
    }
    

##################################################
############ build the genome db blast , transform to bed, sort, merge tag redundancy(python) ########
##################################################
Blastgenome(){
#makeblastdb -in $genome -dbtype nucl -out Index/$SpeciesName

#blastn -db Index/$SpeciesName -query $TEs -outfmt "6 std qlen slen" -evalue 20 -num_threads $N -out AllSpecies/AllTEs/$SpeciesName.$TEfile.blastn
echo 'Blast done'
}


RMgenome(){
OUTFILE=RepeatMasker/$SpeciesName".out.gff"
if [ ! -f "$FILE2" ]; then
    echo OUTFILE" does not exist"
    RepeatMasker -pa 8 -s -gff -no_is -nolow -norna -div 40 -dir RepeatMasker -lib $TEs Genomes/$SpeciesName.fasta
fi
echo 'RM done'
}

MakeBed(){
    echo '>>> Make Bed'
    python3 Scripts/Gff2Bed.py RepeatMasker/$SpeciesName.fasta.out.gff AllSpecies/Bed/
    # output AllSpecies/AllTEs/$SpeciesName.$TEfile.blastn.bed
    echo '>>> Sort'
    ~/bin/bedtools sort -i AllSpecies/Bed/$SpeciesName.fasta.out.gff.bed > AllSpecies/Bed/$SpeciesName.fasta.out.gff.bed.sorted
}


#if [[ $1 == keep ]]; then echo "Skipping RM"; else RMgenome; fi

ExtractSeq()
{
echo ">>> Extract Sequences TE by TE"
> AllSpecies/StatsSeq/$SpeciesName.CpNb.txt
for te in $TElist; do
    (
    #size=$(grep -w $te Data/RepBaseDmel.fa.length |cut -f2)
    size=$(awk -v target="$te" -F'\t' '$1 == target {print $2}' Data/RepBaseDmel.fa.length)
    cutof=0.7
    t=$(echo $size*$cutof | bc)
    ### This separate te, then merge, then filter for size
    ### Bedtools merge give an error when nothing is provided by the pipe
    ### Supposedly due to the use of Bedtool 2.30.0 and the argument -c
    filteredfile=AllSpecies/Bed/$te.$SpeciesName.fasta.out.gff.bed.filtered
    mergefile=AllSpecies/Bed/$te.$SpeciesName.fasta.out.gff.bed.merged
    filteredfasta=AllSpecies/Fasta/$te.$SpeciesName.filtered.fa
    
    # filter for the te
    awk -v te="$te" -v size="$t" 'BEGIN{OFS="\t"}{if($4==te)print $0}' AllSpecies/Bed/$SpeciesName.fasta.out.gff.bed.sorted > $filteredfile
    
    if [ ! -s "${filteredfile}" ];then
        rm AllSpecies/Bed/$te.$SpeciesName.fasta.out.gff.bed.filtered
    else
        # merge bed and filter for size
        ~/bin/bedtools merge  -s -d 250 -c 4,5,6,7,8 -o distinct,mean,distinct,min,max -i $filteredfile | awk -v te="$te" -v size="$t" 'BEGIN{OFS="\t"}{if($3-$2>size)print $0}'> AllSpecies/Bed/$te.$SpeciesName.fasta.out.gff.bed.merged
        # if copies, count copies and extract seq
        if [ -e "${mergefile}" ] && [ -s "${mergefile}" ];then
            Initialcount=$(wc -l "${mergefile}" )
            
            # Get fasta
            ~/bin/bedtools getfasta -s -name -fi /mnt/Poles/Genomes/huavan/Genomes/$SpeciesName.fasta -bed $mergefile > AllSpecies/Fasta/$te.$SpeciesName.fasta

            # Aligning
            mafft --quiet --thread 10 --keeplength  --addfragments AllSpecies/Fasta/$te.$SpeciesName.fasta TEcons/$te.cons.fa > AllSpecies/Fasta/$te.$SpeciesName.mafft
            
            # Get rid of the consensus, Degap, Filter again
            awk '/^>/ {if (++n == 2) p=1} p' AllSpecies/Fasta/$te.$SpeciesName.mafft | bioawk -c fastx '{gsub(/-/, "", $seq); print ">"$name"\n"$seq}' | bioawk -c fastx -v threshold="$t" '{if (length($seq) >= threshold) print ">"$name"\n"$seq}'  > AllSpecies/Fasta/$te.$SpeciesName.filtered.fa
            # write stats (Initial number of sequences per te, and after filtering)
            # Should be in StatsSeq
            if  [ -e "${filteredfasta}" ];then
                count=$(grep -c '^>' "$filteredfasta" |awk 'BEGIN{FS=OFS="\t"} {sub(/[[:space:]].*/, "", $1)} 1')
                echo -e "$SpeciesName\t${te}\t${Initialcount}\t${count}"  >> AllSpecies/StatsSeq/$SpeciesName.CpNb.txt
            fi
        else
            rm AllSpecies/Bed/$te.$SpeciesName.fasta.out.gff.bed.merged
        fi
    fi
    ) &

   # allow to execute up to $N jobs in parallel
    if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
        #now there are $N jobs already running, so wait here for any job
        #to be finished so there is a place to start next one.
        wait -n
    fi

    done

# no more jobs to be started but wait for pending jobs
# (all need to be finished)
wait

echo "All done"
}


MapNAnalyse(){
    echo '>>> Map reads'
    # Build bowtie Index
    bowtie-build --quiet -f AllSpecies/Bowtie/$SpeciesName.All.filtered.fa AllSpecies/Index/$SpeciesName.AllTEs
    # Map the small RNA
    bowtie -a -v $mismatch -p 50 --sam --no-unal -q -x AllSpecies/Index/$SpeciesName.AllTEs $smallRNAreads 2>> AllSpecies/Stats/$SpeciesName.$sra.$mismatch.stats.txt | samtools sort | bamToBed > AllSpecies/Bowtie.$sra.$mismatch/$SpeciesName.sam.bed
    echo 'Mapping done' #OK
    # Compute the number of reads for the sequence with the most read for each TE
    awk 'BEGIN{OFS="\t"} {split($1, arr1, "::"); split($4, arr2, "_x"); print arr1[1], arr1[2], arr2[2], $2, arr2[1], $5, $6}' AllSpecies/Bowtie.$sra.$mismatch/$SpeciesName.sam.bed | awk 'BEGIN{OFS="\t"} {sum[$1"\t"$2] += $3} END{for (key in sum) {split(key, arr, "\t"); print arr[1], arr[2], sum[key]}}' |awk 'BEGIN{OFS="\t"} {sum[$1"\t"$2] += $3; if ($3 > max[$1]) {max[$1] = $3; line[$1] = $0}} END{for (key in line) print line[key]}' | sort > AllSpecies/Bowtie.$sra.$mismatch/$SpeciesName.readcount.txt
    
    echo 'Read analysis done' #OK
}


##### Pipeline starts here

######################################################
########### Step1: BLAST or RM
########### Step2: Transform gff or blast into bed file
########### Step3: Sequence extraction from bed
######################################################

if [[ $5 == keep ]]
then
    echo "Skipping first part"
elif [[ $5 == All ]]
then
    #RMgenome 
    MakeBed
    ExtractSeq #te by te
    # Merge bed per species
    cat AllSpecies/Bed/*.$SpeciesName.fasta.out.gff.bed.merged > AllSpecies/StatsSeq/$SpeciesName.fasta.out.gff.bed.merged
    # The All.filtered.afa used for mapping
    cat AllSpecies/Fasta/*.$SpeciesName.filtered.fa > AllSpecies/Bowtie/$SpeciesName.All.filtered.fa
else
    RMgenome 
    MakeBed
    ExtractSeq #te by te
    # Merge bed per species
    cat AllSpecies/Bed/*.$SpeciesName.fasta.out.gff.bed.merged > AllSpecies/StatsSeq/$SpeciesName.fasta.out.gff.bed.merged
    # The All.filtered.afa used for mapping
    cat AllSpecies/Fasta/*.$SpeciesName.filtered.fa > AllSpecies/Bowtie/$SpeciesName.All.filtered.fa
    
fi

##########################################
########### Step4: Map the piRNA
##########################################
MapNAnalyse

##########################################
################ Saving stat ################
##########################################
echo -e 'Genome' $(LL_ALL=C grep '>' $genome|wc -l) 'sequences' >> $statfile
echo -e 'TEfile' $(LL_ALL=C grep '>' $TEs|wc -l) 'sequences' >> $statfile
echo -e 'SRA' $(echo $(LL_ALL=C cat $smallRNAreads|wc -l) /4|bc) 'sequences' >> $statfile
echo -e '\n' >> $statfile
# StatsSeq
echo -e 'Total blastn/RM hits:' $(wc -l AllSpecies/Bed/$SpeciesName.fasta.out.gff.bed) >> $statfile
# StatsSeq
echo -e 'Merged hits: ' $(wc -l AllSpecies/StatsSeq/$SpeciesName.fasta.out.gff.bed.merged) >> $statfile

##### Mapping statistics
Analyse(){
    # list all the reads mapped per species (non decollapsed)
    # should be stored into a folder StatsMap!
    cut -f 4 AllSpecies/Bowtie.$sra.$mismatch/$SpeciesName.sam.bed | sort|uniq  > AllSpecies/Stats/$SpeciesName.bowtie.$sra.$mismatch.AllReads.uniq
    
    # print the min and max divergence found with RM (5th column of bed.merged)
    # should be stored in a file stored into StatsSeq!
    cut -f 5 AllSpecies/Bed/*.$SpeciesName.fasta.out.gff.bed.merged | sort -n | tee >(echo "min=$(head -1)") \ > >(echo "max=$(tail -1)")
    
}

Analyse



