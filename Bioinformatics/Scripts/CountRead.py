#!/usr/bin/env python

import sys, os
from collections import defaultdict

""" 
Open the file with all concatenated bowtie bed results
For each read, report the TE families mapped

"""


def Rereplicate(listeRead):
    if len(listeRead)==0:
        ReadListDerep=listeRead       
    elif '_x' in listeRead[0]:
        DerepList=[[read] * int(read.split('_x')[1]) for read in listeRead]
        ReadListDerep=[num for sublist in DerepList for num in sublist]
    else : ReadListDerep=listeRead
    print('Decollapsed:',len(ReadListDerep))
    return(ReadListDerep)
    
    
def CountRead(ConcatBowtieFile):
 
    # open the file as list of lines
    bowtie = open(ConcatBowtieFile).readlines()
    AllReads=[]
    print(len(bowtie), 'TEs families were mapped by small RNAs')
    dictRead=defaultdict(list)
    for line in bowtie:
        line=line.rstrip().split(' ')
        AllReads.extend(line[1:])
        for read in line[1:] : 
            dictRead[read].append(line[0])
    print('initial collapsed: ', len(AllReads))   
    print('unique collapsed: ', len(list(set(AllReads)) )  )
       
    AllReadsDR=Rereplicate(list(set(AllReads)))
    print('Decollapsed from non unique:', len(Rereplicate(AllReads)))
            
    with open(ConcatBowtieFile+'.txt','w') as output:
        for read in dictRead:
            output.write('\t'.join(map(str,[read, len(dictRead[read]), ','.join(dictRead[read])]))+'\n')
    
    
CountRead(sys.argv[1])

# 