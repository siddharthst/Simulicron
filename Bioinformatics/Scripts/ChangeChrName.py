#!/usr/bin/env python3


import sys, os

def ChangeBedChrName(bedfilename, dicoName):
    listeChr=['2L','2R','3L','3R','4','Y','X','rDNA','mitochondrion']
    
    # Open bedfile and find name to search and replace
    bedfile=open(bedfilename).readlines()
    clusterchr=[line.split('\t')[0] for line in bedfile if line.split('\t')[0] not in listeChr]
    clusterchr=list(set([line.split('v1')[0].split('_')[1] for line in clusterchr]))
    bedfile2=[]
    for i,line in enumerate(bedfile):
        line=line.split('\t')
        if line[0] not in listeChr:
            line[0]=line[0].split('v1')[0].split('_')[1]
        bedfile2.append('\t'.join(line))
    
    
    Savefile=open(bedfilename+'.mod2','w')
    for i,line in enumerate(bedfile2):
         for code in dicoName.keys():
             if code in line:
                 bedfile2[i]=line.replace(code,dicoName[code])
         Savefile.write(bedfile2[i])
    Savefile.close()
    
def ChangeGenomeChrName(genomefile):
    genome = open(os.path.expanduser(genomefile)).read()
    chromosome=genome.split('>')[:]
    genomeSave = open(os.path.expanduser(genomefile)+'.mod','w')
    Chrname=[chromo.split('\n')[0] for chromo in chromosome[1:]]
    ChrSeq=[''.join(chromo.split('\n')[1:]) for chromo in chromosome[1:]]    
    #print(Chrname)
    dicoName={}
    for i,Chr in enumerate(Chrname):
        modified_name=''.join(Chr.split()[-2:]).replace('chromosome','').replace('sequence','').replace('completegenome','mitochondrion').strip()
        dicoName[Chr.split(' ')[0].split('.')[0]]=modified_name
        genomeSave.write('>'+modified_name+'\n'+ChrSeq[i]+'\n')        
    return(dicoName)

if len(sys.argv)==2:
    ChangeGenomeChrName(os.path.expanduser(sys.argv[1]))
else:
    dico=ChangeGenomeChrName(os.path.expanduser(sys.argv[2]))
	ChangeBedChrName(os.path.expanduser(sys.argv[1]),dico)


#python Scripts/ChangeBedChrName.py Data/piCluster.dmelr6.r2.bed Data/dmel-all-chromosome-r6.38.fasta