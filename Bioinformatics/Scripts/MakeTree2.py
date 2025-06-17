#!/usr/bin/env python

### Siddharth S. Tomar
### Aurélie Hua-Van
### Arnaud Le Rouzic


# This script extracts the best TE copy for each TE for each species
# It gather into one fasta file all the copies for the same TE
# Input: $sra.$mismatch.crossRegulated.txt (all species-TE pairs)
# Input: $Species.readcount.txt file for each species (to get best copy's name)
# Input: $TE.$Species.mafft (all seq for this TE in this species) (to get best copy's seq)
# Output: $TE.$sra.$mismatch.2.fasta (aligned seq per TE for making a tree)

import glob
from collections import defaultdict
import os, sys

dictCross=defaultdict(list) # dictCross[Species]=list(TE)

BestCopiesDict=defaultdict(list)

# For each species, get sequence of best copy saved in BestCopiesDict[TE]=[fasta of best seq for each species]
def GetBestCopySeq(Speciesfile):
    
    Species=os.path.basename(Speciesfile).split('.readcount.txt')[0]
    readcount=open(Speciesfile).readlines()
    
    # Extract the copy on the form TE::contig:posX-PosY+)
    bestcopyList=[line.split('\t')[0]+'::'+line.split('\t')[1] for line in readcount]
    
    # USELESS
    # Open the fasta file (filtered)
    #fasta=open('AllSpecies/Bowtie/'+Species+".All.filtered.fa").read().split('>')[1:]

    # Open the fasta aligned file corresponding to the species and the TE    
    # Put the fasta sequence into the BestCopiesDict[TE]=[fasta of best seq for each species]
    for copie in bestcopyList:
        TE=copie.split('::')[0]
        fasta=open('AllSpecies/Fasta/'+TE+'.'+Species+".mafft").read().split('>')[1:]
        for seq in fasta:
            if copie in seq:
                BestCopiesDict[copie.split('::')[0]].append('>'+Species+'\n'+''.join(seq.split('\n')[1:]))

print(sys.argv[1])

# Extract parameters
mismatch=os.path.basename(sys.argv[1]).split('.')[-3]
sra=os.path.basename(sys.argv[1]).split('.')[-4]
print(sra,mismatch)

# open the file with the pairs species-TE (noref)
CrossReg=open(sys.argv[1]).readlines()[1:]
# Fill the dict[Species]=[TE1,TE2,...]
for line in CrossReg:
    line=line.rstrip().split('\t')
    dictCross[line[1]].append(line[0])

#recountFolder="AllSpecies/Bowtie."+sra+'.'+mismatch+'/'
#for Speciesfile in glob.glob(eval("'"+recountFolder*.readcount.txt"'")):

# For each species, open the readcount.txt file that contains the piRNA count for the best copy per TE
for Species in dictCross.keys():
    GetBestCopySeq("AllSpecies/Bowtie."+sra+'.'+mismatch+'/'+Species+".readcount.txt")

# for each TE family, write into a new file in Tree folder
for TE in BestCopiesDict.keys():
    print(TE)
    with open('AllSpecies/Tree/'+TE+'.'+sra+'.'+mismatch+'.2.fasta','w') as TEfile:
        TEfile.write('\n'.join(BestCopiesDict[TE]))
    
# python3 Scripts/MakeTree.py AllSpecies/Bowtie.$sra.$mismatch/*.readcount.txt

# python3 Scripts//MakeTree2.py AllSpecies/Results/Species.TE.$sra.$mismatch.CrossRegulated.txt 