#!/usr/bin/env python

import glob
from collections import defaultdict
import os, sys


BestCopiesDict=defaultdict(list)
def GetBestCopy(Speciesfile):
    
    Species=os.path.basename(Speciesfile).split('.readcount.txt')[0]
    readcount=open(Speciesfile).readlines()
    bestcopyList=[line.split('\t')[0]+'::'+line.split('\t')[1] for line in readcount]

    fasta=open('AllSpecies/Bowtie/'+Species+".All.filtered.fa").read().split('>')[1:]
    for copie in bestcopyList:
            TE=copie.split('::')[0]
            fasta=open('AllSpecies/Fasta/'+TE+'.'+Species+".mafft").read().split('>')[1:]
            for seq in fasta:
                if copie in seq:
                    BestCopiesDict[copie.split('::')[0]].append('>'+Species+'\n'+''.join(seq.split('\n')[1:]))

print(sys.argv[1])
for Speciesfile in sys.argv[1:]:
#for Speciesfile in glob.glob(eval("'"+sys.argv[1]+"'")):
    print(Speciesfile)
    if 'Dmel.Ref' in Speciesfile: pass
    else:
        GetBestCopy(Speciesfile)
mismatch=os.path.dirname(sys.argv[1]).split('.')[-1]
sra=os.path.dirname(sys.argv[1]).split('.')[-2]
print(sra,mismatch)

for TE in BestCopiesDict.keys():

    with open('AllSpecies/Tree/'+TE+'.'+sra+'.'+mismatch+'.fasta','w') as TEfile:
        TEfile.write('\n'.join(BestCopiesDict[TE]))
    
# python3 Scripts/MakeTree.py AllSpecies/Bowtie.$sra.$mismatch/*.readcount.txt