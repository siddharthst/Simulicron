#!/bin/bash

### Siddharth S. Tomar
### Aurélie Hua-Van
### Arnaud Le Rouzic

### cd to Project
Project=~/CrossRegulation
####################################
## Parameters
###################################
Project=$1
sra=$2
mismatch=$3
method=$4
keep=$5

cd $Project

####################################
# Get the sequences from Dmel genome
####################################
if [[ ! $5 == keep ]]; then
# Make Bed
python3 Scripts/Gff2Bed.py RepeatMasker/D.melanogaster.fasta.out.$method AllSpecies/D.melanogaster/

bedtools sort -i AllSpecies/D.melanogaster/D.melanogaster.fasta.out.$method.bed > AllSpecies/D.melanogaster/D.melanogaster.$method.sorted.bed

# Extract seq TE by TE
TEs=Data/RepBaseDmel.fa
TElist=$(cat $TEs | grep '>' |cut -f 1 | cut -c2- |uniq)

for te in $TElist; do

awk -v te="$te" -v size="$t" 'BEGIN{OFS="\t"}{if($4==te)print $0}' AllSpecies/D.melanogaster/D.melanogaster.$method.sorted.bed | bedtools merge -s -d 100 -c 4,5,6,7,8 -o distinct,mean,distinct,min,max |awk -v te="$te" -v size="$t" 'BEGIN{OFS="\t"}{if($3-$2>size)print $0}'> AllSpecies/D.melanogaster/D.melanogaster.$method.$te.merged.bed

bedtools getfasta -s -name -fi /mnt/Poles/Genomes/huavan/Genomes/D.melanogaster.fasta -bed AllSpecies/D.melanogaster/D.melanogaster.$method.$te.merged.bed > AllSpecies/D.melanogaster/TE.D.melanogaster.$method.$te.fasta

#mafft --quiet --keeplength  --addfragments AllSpecies/D.melanogaster/TE.D.melanogaster.$method.$te.fasta temp/$te.cons.fa > AllSpecies/D.melanogaster/$te.Dmel.Ref.mafft

#awk '/^>/ {if (++n == 2) p=1} p' AllSpecies/D.melanogaster/$te.Dmel.Ref.mafft | bioawk -c fastx '{gsub(/-/, "", $seq); print ">"$name"\n"$seq}' > AllSpecies/D.melanogaster/TE.D.melanogaster.$method.$te.degap.fasta
done

# Concatenate into one file
cat AllSpecies/D.melanogaster/D.melanogaster.$method.*.merged.bed  | bedtools sort > AllSpecies/D.melanogaster/Dmel.Ref.$method.merged.bed
cat AllSpecies/D.melanogaster/TE.D.melanogaster.$method.*.fasta > AllSpecies/D.melanogaster/Dmel.Ref.$method.TE.fasta

# build the index
bowtie-build --quiet -f AllSpecies/D.melanogaster/Dmel.Ref.$method.TE.fasta AllSpecies/D.melanogaster/Dmel.Ref.$method.TE

fi


####################################
# Map piRNA with Bowtie 
####################################
echo ">>> mapping Reads"
bowtie -a -v $mismatch -p 70 --sam --no-unal -q -x AllSpecies/D.melanogaster/Dmel.Ref.$method.TE Data/$sra.trimmed.collapse.fq  2>> AllSpecies/D.melanogaster/Dmel.Ref.$method.TE.$sra.$mismatch.stats.txt | samtools sort | bamToBed > AllSpecies/D.melanogaster/Dmel.Ref.$method.TE.$sra.$mismatch.sam.bed

awk 'BEGIN{OFS="\t"} {split($1, arr1, "::"); split($4, arr2, "_x"); print arr1[1], arr2[1], arr2[2]}' AllSpecies/D.melanogaster/Dmel.Ref.gff.TE.$sra.$mismatch.sam.bed |sort |uniq | awk 'BEGIN{OFS="\t"} {sum[$1] += $3} END{for (key in sum) {print key, "AllCopies", sum[key]}}'  > AllSpecies/D.melanogaster/Dmel.Ref.gff.TE.$sra.$mismatch.readcount.txt 

# copy to bowtie directory
cp AllSpecies/D.melanogaster/Dmel.Ref.gff.TE.$sra.$mismatch.readcount.txt AllSpecies/Bowtie.$sra.$mismatch/

mv AllSpecies/Bowtie.$sra.$mismatch/Dmel.Ref.gff.TE.$sra.$mismatch.readcount.txt AllSpecies/Bowtie.$sra.$mismatch/Dmel.Ref.readcount.txt



####################################
# Get statistics
####################################
# from all copies
echo Number of sequences
grep -c '>' AllSpecies/D.melanogaster/Dmel.Ref.$method.TE.fasta

echo "3- Total number of mapping"
wc -l AllSpecies/D.melanogaster/Dmel.Ref.$method.TE.$sra.$mismatch.sam.bed

echo "4-Total number of Reads collapsed"
cut -f4 AllSpecies/D.melanogaster/Dmel.Ref.$method.TE.$sra.$mismatch.sam.bed | sort |uniq |wc -l

echo "5- Unique decollapsed reads on all Copies"
cut -f4 AllSpecies/D.melanogaster/Dmel.Ref.$method.TE.$sra.$mismatch.sam.bed | sort |uniq|awk '{split($1, arr1, "_x"); print arr1[2]}' | paste -sd+ | bc

# Statistics from D.mel in the whole Drosophila species analysis
echo "6- Unique decollapsed reads on 70 percent"
cd $Project/AllSpecies/Bowtie.$sra.$mismatch/
cut -f4 D.melanogaster.sam.bed | sort | uniq| awk '{split($1, arr1, "_x"); print arr1[2]}' | paste -sd+ | bc
cd $Project

echo "7- Non unique decollapsed reads on best copy"
cd $Project/AllSpecies/Bowtie.$sra.$mismatch/
cut -f3 D.melanogaster.readcount.txt | paste -sd+ | bc
wc -l  D.melanogaster.readcount.txt
cd $Project

echo "8- Unique decollapsed reads on the best copy"
cd $Project/AllSpecies/Bowtie.$sra.$mismatch/
awk '{print $1 "::" $2}' D.melanogaster.readcount.txt > /tmp/contig.txt
LL_C=C grep -wFf /tmp/contig.txt D.melanogaster.sam.bed | cut -f4 |sort| uniq|awk '{split($1, arr1, "_x"); print arr1[2]}'  | awk '{s+=$1} END {print s}'
cd $Project



echo "Done"
