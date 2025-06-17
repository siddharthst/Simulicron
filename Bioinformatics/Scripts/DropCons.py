#!/usr/bin/env python

import sys, os

""" 
Clean the aligned msa by removing cons
and every non kept pi sequences (not enough piRNA)
before tree construction
"""

##### TO DO: recalculate distance [DONE]
### Put it on coord.txt

def GetGap(name,seq,cons):
    # name of the sequence, start, end, div with consensus (0-1)
    # one line per segment
    seq=seq.upper()
    div=GetPercentID(name,seq,cons)
    SegStart=[]
    SegEnd=[]
    retline=''
    if seq[0]!='-':SegStart.append(0)
    for (i, nt) in enumerate(seq[:-1]):
        if nt=='-' and seq[i+1] != '-': SegStart.append(i+1) # '-N'
        elif nt!='-' and seq[i+1] == '-': SegEnd.append(i) # 'N-'
    if seq[-1]!='-':SegEnd.append(len(seq))
    for i,seg in enumerate(SegStart):
        retline=retline+name+'\t'+str(seg)+'\t'+str(SegEnd[i])+'\t'+div+'\n'
    return(retline)
    
def GetPercentID(name,seq,cons):
    
    seq=seq.upper()
    cons=cons.upper()
    seqwoN=seq.replace('-','')
    n=0
    for i,j in enumerate(seq):
        if seq[i]==cons[i] and seq[i]!='-': n+=1
    return(str(round(n/len(seqwoN),2)))
    

def DropCons(fastafile, keptfile, TE):
 
    # Drop the cons and kept seq as list
    fasta = open(fastafile).read().split('>')[1:]
    cons=''.join(fasta[0].split('\n')[1:])
    fasta=fasta[1:] # List all seqs except cons
    fastaname=[seq.split('\n')[0] for seq in fasta]# List names of all seqs except cons
    fastaseq=[''.join(seq.split('\n')[1:]) for seq in fasta]
         
    # list of seq kept after bowtie (>200 reads)
    if "kept" in keptfile:
    	# line 4 contains all names of sequences to keep (exclude discarded copies)
        #keptSeq  = open(keptfile).readlines()[4].split('\t')
        # comment the previous and uncomment the following if we do not want to exclude pi copies not expressed.
        keptSeq = fastaname
    elif "discarded" in keptfile:
    	# We keep all copies 
        keptSeq = fastaname
    # var= sra.dist.mm
    var = os.path.basename(keptfile).split('.bed')[0].split('.'+TE+'.')[1]

    # Output
    newfasta = open(os.path.dirname(fastafile)+'/'+TE+'.'+var+'.afa','w')
    seq_coord = open(os.path.dirname(fastafile)+'/'+TE+'.'+var+'.coord.txt','w')
    
    # kept seq 2 keep
    for i,seq in enumerate(fastaname):
        if seq in keptSeq:
            newfasta.write('>'+seq+'\n'+fastaseq[i]+'\n')            
            seq_coord.write(GetGap(seq,fastaseq[i],cons))
    newfasta.close()
    seq_coord.close()
    
DropCons(sys.argv[1], sys.argv[2], sys.argv[3])

# input1 .cons.afa (after mafft on renamed.filtered)
# input2 /D.melanogaster.r6.TE.sra.dist.mm.bed.{kept/discarded}.txt(after mafft on renamed.filtered)
# output1 .afa to be used for tree
# output2 .coord.txt with all info for all insertions in the sequences
    