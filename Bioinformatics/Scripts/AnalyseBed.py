#!/usr/bin/env python


import sys, os
from collections import defaultdict
from collections import Counter


""" 
After mapping with bowtie and transformation into bed
Generate the matrix to plot the figure
Update the summaries file
## NO MORE: At this step, pi Copies are excluded if they map less than 1/10000 of the total numer of reads (2240 reads)
output:
List2print

"""

def CountReadByLength(listofhit,TE):
    #### We should filter for name and strand
    #### then unique, if both direction, take the most frequent
    #### get the size from the sequence or bowtie without indel
    #### Dereplicate the read for counting
    
    ### or 
    ### do a tuple [(name, length, strand)]
    ### Count the tuples (Counter)
    ### Make unique
    ### 


    ## read 1 an 2 + coordinates
    ## read 3 name of the read
    ## read 5 strand
    
    # One read SRR.123_x3 is 25 nt mapping + [mapping (4 28) (14 38) (1 25) (300 324)]
    # One read SRR.456_x2 is 20 nt mapping + [mapping (20 40)]
    # One read SRR.789_x1 is 25 nt mapping + [mapping (104 128)]
    # One read SRR.357_x4 is 20 nt mapping - [mapping (104 128)]
    
    # transform read into more compact: SizeStrand.read_x3
    ReadTuple2calculate=[ str(int(read[2])-int(read[1])+1)+read[5]+'.'+read[3] for read in listofhit]  
    # Make unique
    resunique=list(set(ReadTuple2calculate))  
    # Rereplicate
    res=Rereplicate(resunique)
    # keep SizeStrand
    res2= [r.split('.')[0] for r in res]
    
    
#     # Previously
#     ReadList2calculate=[ str(int(read[2])-int(read[1])+1)+'.'+read[5]+'_x'+read[3].split('_x')[1] for read in listofhit]   
#     res=Rereplicate(ReadList2calculate)
#     # We send to Rereplicate:
#     # 25.+_x3
#     # 20.+_x2
#     # 25.+_x1
#     # 20.-_x1
#     # Dereplication should give 25.+_x3 25.+_x3 25.+_x3 20.+_x2 20.+_x2 25.+_x1 20.-_x4 20.-_x4 20.-_x4 20.-_x4
#     res2= [r.split('_x')[0] for r in res]
#     # 25+ 25+ 25+ 20+ 20+ 25+ 20- 20- 20- 20-
#     # We cut the last part (.split('_x')[0]) which has became useless
#     # Then it is:  25 25 25  20 20 25
    
    
    
    # We count with Counter and make it as a dict
    D= dict(Counter(res2))
    # dict of the form : {25+ : 5, 20+ : 1, 20- : 4}
    # We transform the dict into a list of sorted tuples
    # ("25", "+", 3568, "TE")
    listoftuples=[(r[:-1], r[-1] ,str(D[r]),TE) for r in D]
    # We sort by strand then size
    sorted_tuples= sorted(listoftuples,key=lambda x: (x[1],x[0]))
     # we make a big string
    stuff2save= '\n'.join(map(str,['\t'.join(tup) for tup in sorted_tuples]))+'\n'
    
    # We calculate the sum of (-) and the total of decollapsed read
    NbReadsTotS = sum([int(item[2]) for item in sorted_tuples])
    Nbminus = sum([int(item[2]) for item in sorted_tuples if item[1]=='-'])
    if NbReadsTotS !=0 :
        FreqMinus=round(Nbminus/NbReadsTotS,4)
    
    # We calculate the frequence of '-'
    #CounterOri=Counter(elem[1] for elem in listoftuples)
    #if len(listoftuples) !=0 and '-' in CounterOri.keys(): FreqMinus= round(CounterOri['-'] /len(listoftuples),4)
    else: FreqMinus= 'NA'
    return(stuff2save,FreqMinus)
    

#Function to decollapse reads. Need to be used everytime at the last moment
def Rereplicate(listeRead):
    #SRAXXXXXXX.12345_x5
    #size.strand._x3
    if len(listeRead)==0:
        ReadListDerep=listeRead       
    elif '_x' in listeRead[0]:
        DerepList=[[read] * int(read.split('_x')[1]) for read in listeRead]
        ReadListDerep=[num for sublist in DerepList for num in sublist]
    else : ReadListDerep=listeRead
    return(ReadListDerep)
    
def AnalyseBed(bedfile, TE, summaryfile,TotCpNb,threshold):
    ### FOR EACH TE
    # Create the  output files (bed.mat and summary.txt)
    matrix=open(bedfile+'.mat','w')
    summary2 = open(bedfile+'.bowtie.summary.txt','w')
    List2print=[TE,TotCpNb]

    
    ####### Open Bedfile bowtie (reads are collapsed and redundant, should be 23-30)   
    bed=open(bedfile).readlines()  
    # List of reads (collapsed but redundant) come with different columns (
    # [0]: copy name
    # [1]: start
    # [2]: end
    # [3]: read name (_x1)
    # [4]: score (255)
    # [5]: strand
    
    listofhit=[hit.rstrip().split('\t') for hit in bed]

    #### Prepare read length for histogram with orientation
    stuff2save,FreqMinus = CountReadByLength(listofhit, TE)
    with open(bedfile+'.bowtie.histo.txt','w') as histo:
        histo.write(stuff2save)
    
    #### Prepare read length for histogram with orientation    
    #ReadsMappingOri = list(set([(read[3],read[5]) for read in listofhit]))
    #NbPlus=len([i for i in ReadsMappingOri if i[1]=='+'])
    #NbMinus=len([i for i in ReadsMappingOri if i[1]=='-'])
    #if NbMinus+NbPlus != 0: FreqMinus=round(NbMinus/(NbMinus+NbPlus),2)
    #else: FreqMinus='NA'
   
    # Remove redundancy ListOfHits --> ReadsMapping (keep only the name of the read)
    ReadsMapping = list(set([read[3] for read in listofhit])) # Read name
    # Number of collapsed reads 
    NbReadsDerep = len(ReadsMapping)
    # Total number of reads (after decollapsing)
    NbReadsTot = len(Rereplicate(ReadsMapping))

    # List the insertions and the insertions into pi
    ListOfInsertions=list(set([line[0] for line in listofhit])) # Copy Name non redundant
    ListOfPi=[name for name in ListOfInsertions if 'piRNA' in name]
            
    ####### Open summary with info on TE copies, put it into dicoSummary
    ####### If the copy is not one of the mapping --> Discarded
    # name is col 11
    summary= open(summaryfile).readlines()
    dicoSummary=defaultdict(list)    
    for line in summary: # each copy (All copies)
        insertion= line.split('\t')[14] # Copy name Cp.xxxx
        # Transform the line into a dict[name]
        dicoSummary[insertion]=line.rstrip().split('\t')
        if insertion not in ListOfInsertions: # euchromatic short?
            dicoSummary[insertion].extend(["NA","Discarded"])
            summary2.write('\t'.join(dicoSummary[insertion])+ '\n') 


    # open dictionaries (decollapsed)
    dicoExcluded=defaultdict(list)    
    dicoInsertion=defaultdict(list)
    dicoPi=defaultdict(list)
    
    # List of all reads (decollapsed, redundant)
    ListRead=[]
    # List of all reads (decollapsed, but redundant if on several copies) mapping piRNA sequence
    ListReadPi=[]

    
    
    nFL,neuk,npi,nex=0,0,0,0
    flag=''
    for insertion in ListOfInsertions:
    
        # List of unique collapsed reads
        ReadListTE=list(set([read[3] for read in listofhit if read[0]==insertion]))
        # List of decollapsed reads (contains collapsed names in duplicates)
        ReadListTEDecol=Rereplicate(ReadListTE)
        
        ## We keep insertion if it is into piCluster and mapped with "threshold" reads or more
        # piRNA above threshold (here =0)
        if 'piRNA' in insertion and len(ReadListTEDecol)>=0:  # No threshold
            # list of reads not decollapsed, mapped to all the insertion in (pi)
            ListReadPi.extend(ReadListTE)
            # Decollapsed reads are listed in the dict with pi, and the dict with tot
            dicoPi[insertion]=ReadListTEDecol
            npi+=1
            dicoInsertion[insertion]=ReadListTEDecol
            nFL+=1
            dicoSummary[insertion].append(str(len(ReadListTEDecol)))
            dicoSummary[insertion].append('Kept')
            summary2.write('\t'.join(dicoSummary[insertion])+ '\n') 
            
        ## We also keep insertion if it is not into piCluster
        ## EUCHROMATIC
        elif 'piRNA' not in insertion :
            # list of reads not decollapsed, mapped to all the insertion in pi or not
            ListRead.extend(ReadListTE)
            # Decollapsed reads are listed in the dict with tot only
            dicoInsertion[insertion]=ReadListTEDecol
            neuk+=1
            nFL+=1
            dicoSummary[insertion].append(str(len(ReadListTEDecol)))
            dicoSummary[insertion].append('Kept')
            summary2.write('\t'.join(dicoSummary[insertion])+ '\n') 
            
        ## This is the piRNA insertion with less than threshold reads
        else: 
            dicoSummary[insertion].append(str(len(ReadListTEDecol)))
            dicoSummary[insertion].append('Discarded')
            summary2.write('\t'.join(dicoSummary[insertion])+ '\n') 
            nex+=1   
    summary2.close()
    
    # Remove redundancy due to sharing between copies
    UniqRead=list(set(ListRead))
    UniqPiRead=list(set(ListReadPi))
    
    
    #### Dereplicate the list of reads and Decollapse it for having total
    DecolUniqReads=Rereplicate(UniqRead)
    DecolUniqPiReads=Rereplicate(UniqPiRead)
    
    # Calculate the length
    NbReadsFromKeptNonPi=len(DecolUniqReads)
    NbReadsFromKeptPi=len(DecolUniqPiReads)
    
    # Total decollapsed reads mapping kept copies
    NbReadsKept=len(Rereplicate(list(set(ListReadPi+ListRead))))
    
    #### SHARED: Decollapsed, unique
    # List of reads (collapsed) mapping at least once both pi and non pi copies
    Common=list(set(ListReadPi).intersection(set(ListRead)))
    # Total decollapsed reads mapping both pi and nonpi copies (shared)
    CommonRerep=len(Rereplicate(Common))
    
    #### for the copies
    FiltCpNb=len(ListOfInsertions) # mapped with bowtie
    KeptCpNb=len(dicoInsertion) # kept (excluding piRNA with less than 200 decollapsed pi ) [should equal nFL]
    KeptPiCpNb=len(dicoPi) #[should equal npi]
    KeptEucCpNb=KeptCpNb-KeptPiCpNb
    ExcPi=str(nex) 
    
    List2print.extend([FiltCpNb,KeptCpNb,KeptPiCpNb,KeptEucCpNb,ExcPi])
    List2print.extend([NbReadsDerep,NbReadsTot,FreqMinus,NbReadsKept,NbReadsFromKeptNonPi,NbReadsFromKeptPi,CommonRerep])
    
    ### Filtering TE families with at least 6 copies including 2 piClusters
    ### A file "Kept" or "Discarded" is created. 
    ### Trees will be constructed only if the "Kept" file exist for the TE
    ### R will analyse only TEs for which a tree exist
    #if nFL >= 6 and npi >= 2 : flag='kept' # 7 FL copies not pi  and 3 pi copies
    if nFL >= 4 and npi >= 1 and neuk >= 1 : flag='kept' # 7 FL copies not pi  and 3 pi copies
    else: flag='discarded' 
    
    with open(bedfile+'.'+flag+'.txt','w') as TE2write:
        TE2write.write('Total insertions: ' + str(len(ListOfInsertions))+'\n')
        TE2write.write('Kept insertions: '  + str(len(dicoInsertion))+'\n')
        TE2write.write('Kept pi insertion: '+ str(len(dicoPi))+'\n')
        TE2write.write('Excluded pi insertions: ' + str(nex)+'\n')
        TE2write.write('\t'.join(dicoInsertion.keys()))
            
    
    ### Create matrix with all insertions
    maxPi=0
    maxte=0
    cellpi='NA'
    cellte='NA'
    matrix.write("\t"+'\t'.join(sorted(dicoPi.keys()))+'\n')
    for te in dicoInsertion:
        line2write=[te]
        for pi in sorted(dicoPi.keys()):
            sharedReads=[read for read in set(dicoInsertion[te]).intersection(dicoPi[pi])]
            shareReadsDerep=Rereplicate(sharedReads)
            line2write.append(str(len(shareReadsDerep)))
            maxPi=max(maxPi,len(shareReadsDerep))
            if max(maxPi,len(shareReadsDerep)) == len(shareReadsDerep): cellpi=pi
            if 'piRNA' not in te: 
                maxte=max(maxte,len(shareReadsDerep))
                if max(maxte,len(shareReadsDerep)) == len(shareReadsDerep): cellte=pi
        matrix.write('\t'.join(line2write)+ '\n')
    matrix.close()
    List2print.extend([maxPi,cellpi,cellte,flag])
    
    # To be save in a summary file (bash) (TE.summary.txt)
    print('\t'.join(map(str,List2print)))
    
AnalyseBed(sys.argv[1], sys.argv[2],sys.argv[3], sys.argv[4], sys.argv[5])
# 1: bedfile for the te
# 2 te name
# 3 te summary.txt (genomic copies)
# 4 TETotCpNb (before filtration)
# 5 threshold (1/100000 of the total decollapsed reads)

#print(res)
#AnalyseBed('Results/Bowtie/D.melanogaster.r6.Gypsy1_DM.SRR14569563.0.bed')
#python Scripts/AnalyseBed.py Results/Bowtie/D.melanogaster.r6.Gypsy1_DM.SRR14569563.0.bed Gypsy1_DM  Results/Bed/Gypsy1_DM.summary.txt TETotCpNb threshold