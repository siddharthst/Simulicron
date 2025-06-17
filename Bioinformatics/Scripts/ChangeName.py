#!/usr/bin/env python


import sys, os

""" 
After identification of copies into piCluster (bedtools intersect),
The names of copies in the fasta file is changed by a number with ou without the mention
whether it is in piCluster or not (renamed)
Another fastafile is created (renamed.filtered.) containing size-filtered copies
A Bed/te.$distance.$pi.summary.txt file is created, with information for each copy
"""


def ChangeName(fileroot,piAnnot,prop):
    
    # Path to files (Bed folder)
    # fileroot=$te.$distance
    piInsertionFile=fileroot+'.'+piAnnot+'.intersect'+'.bed'
    overlapfile=fileroot+'.intersect.overlap.bed'
    bedfile=fileroot+'.merged.bed'

    fastafile=fileroot+'.fa'
    # $te.cons.seq
    consfile=os.path.splitext(fileroot)[0]+'.cons.seq'
    te=os.path.basename(fileroot)
    
    # Output file (in Bed folder)
    #fastafile = Bed/#te.$distance.fa
    fastacor=open(fileroot+'.'+piAnnot+'.renamed.fa','w')
    #fastacor=open(os.path.splitext(fastafile)[0]+'.'+piAnnot+'.renamed.fa','w')
    fastafin=open(os.path.splitext(fastafile)[0]+'.'+piAnnot+'.renamed.filtered.fa','w')
    summary=open(os.path.splitext(fastafile)[0]+'.'+piAnnot+'.summary.txt','w')
    
    fasta=open(fastafile).read().split('>')[1:]
    piInsertion=open(piInsertionFile).readlines()
    consseq=open(consfile).read().split('\n')[1].replace('-','')
    overlap=open(overlapfile).readlines()
    bed=open(bedfile).readlines()
    
    # for all copies after merge (bed not filtered)
    dicoBed={}
    listDiv=[]

    for j in bed:
        i=j.rstrip().split('\t')
        div=list(map(float,i[8].split(',')))
        evalue=list(map(float,i[9].split(',')))
        ###### Need to recalculate pairwise identity with the consensus (have to be done before the tree)
        divkept=min([div[i] for i,j in enumerate(evalue) if j== min(evalue)])
        #Take the first div of the min evalue
        #divkept=div[evalue.index(min(evalue))]
        percentstart=round(int(i[6])/len(consseq),2)

        name=i[0]+':'+i[1]+'-'+i[2]
        if name in dicoBed: continue # avoid
        #dicoBed[name]=i[:8]+[str(min(evalue)),str(divkept),str(percentstart)]
        dicoBed[name]=i[:8]+[str(min(evalue)),str(divkept),str(i[6])]
        listDiv.append(divkept)
        
        
    meanDiv=round(sum(listDiv)/len(dicoBed),2)
    piName=[i.split('\t')[0]+':'+i.split('\t')[1]+'-'+i.split('\t')[2] for i in piInsertion]
    RedName=[i.split('\t')[0]+':'+i.split('\t')[1]+'-'+i.split('\t')[2] for i in overlap]
   
    # for all copies after merge
    # Change name and filter for size for euchromatic copies
    ListName=[]
    n,m,i,e=0,0,0,0
    for copie in fasta:
        name=copie.split('(')[0]
        if name in ListName: continue
        else: ListName.append(name)
        seq=''.join(copie.split('\n')[1:])
        newname='Cp.'+str(n)
        tab=dicoBed[name]
        if name in piName:
            newname=newname+'_piRNA'
        if name in RedName:
             newname=newname+'_r'
        tab.extend([str(len(seq)), str(len(consseq)),name, newname ])
        fastacor.write('>'+ newname+'\n'+seq+'\n')

        if 'piRNA' in newname:        
            fastafin.write('>'+ newname+'\n'+seq+'\n')
            summary.write('\t'.join(tab)+'\t'+str(meanDiv)+'\tKept\n')
            i+=1
        elif len(seq)>(len(consseq) * float(prop)):
            fastafin.write('>'+ newname+'\n'+seq+'\n')
            summary.write('\t'.join(tab)+'\t'+str(meanDiv)+'\tKept\n')
            m+=1
        else:
            summary.write('\t'.join(tab)+'\t'+str(meanDiv)+'\tDiscarded\n')
            e+=1

        n+=1
    print('\t'.join(map(str,[te,n,m+i,m,i,e])))
        
    fastacor.close()
    summary.close()
    
    
ChangeName(sys.argv[1],sys.argv[2],sys.argv[3])
# python Scripts/ChangeName.py  ResultsDmel/Bed/$te.$distance $threshold1 
# output to ResultsDmel/AllTEs/$genomeVersion.$TEfile.$distance.filter1.summary.txt