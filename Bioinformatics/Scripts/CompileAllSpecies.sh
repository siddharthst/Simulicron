
#!/bin/bash

### Siddharth S. Tomar
### Aurélie Hua-Van
### Arnaud Le Rouzic


################ DEFAULT PARAMETER HERE #############################
Project=~/CrossRegulation
method=gff # for RM analysis (with -gff option), otherwise choose blastn
#sra=SRR11846566
#sra=SRR5687217 # the file and path must be in the form $Project/Data/$sra.trimmed.collapse.fq
#sra=SRR25922470
sra=SRR14569563
mismatch=0
#mismatch=3
flavor=keep
TEs=Data/RepBaseDmel.fa
# flavor=Cp  do only the Cp count
# flavor=keep alternate=keep , for not redoing the first part which is independent of the mapping,
# flavor=All redo everything

# command line
#bash Scripts/CompileAllSpecies.sh $Project $method $sra $mismatch $flavor $TEs 


if [ -n "$1" ]; then
    echo var exists
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


cd $Project


######################## Get reference data (Dmel all copies)
###### need D.melanogaster to be done
bash Scripts/Dmel4AllSpecies.sh $Project $sra $mismatch $method $flavor


################ Phylogenetic tree of TE copie (depends on the sra and mismatch)
# Save the best mapped copy sequence from each genome per te
# Get TE fasta file and save it for phylogeny
python3 Scripts/MakeTree.py AllSpecies/Bowtie.$sra.$mismatch/*.readcount.txt

# Do the tree with best copies
listfile=$(ls AllSpecies/Tree/*.$sra.$mismatch.fasta)
for file in ${listfile[@]}
do
fasttree -quiet -nt -gtr $file > $file.tre
done

############ Save a file with the number of copy per TE for extra analysis
if [ $flavor == All ];then
cd AllSpecies/Bowtie
grep '>' *.filtered.fa | awk '{split($1, arr1, ".All.filtered.fa:>"); print arr1[1], arr1[2]}'| awk '{split($2, arr2, "::"); print $1,arr2[1]}' | awk 'BEGIN{OFS="\t"} {sum[$1"\t"$2] += 1} END{for (key in sum) { print key, sum[key]}}' | sort > ../StatsSeq/All.CpNb.filtered.txt
cd $Project
fi

############## Draw the figure
# Rscript Scripts/AllspeciesHeatMap.R $sra $mismatch


