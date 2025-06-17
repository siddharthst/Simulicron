#!/usr/bin/env/ python

import re,os
import sys
from glob import iglob

nom=''
dom=[]
adel=[]

##################### PROCESSING DU FICHIER D'ACCESSION ########################
def ask():
    '''
    Ce programme mesure la longueur des sequences entre deux points.
    '''
    fastafile=raw_input('length >> Fasta file ? : ').rstrip()
    save=raw_input('length >> Output file ? : [file.length]').rstrip()
    if save=='':save=fastafile+'.length'

    List(fastafile,save)


################## OUVERTURE DU FICHIER D"ACCESSION ############################
def List(fastafile,save) :
    print ("Running length...")
   
    fichier=open(fastafile) 
    file=fichier.read().split('>')

    ############ TRAITEMENT DE CHAQUE ACCESSION ################################
    somme=0
    outfile=open(save,'w') 
    for x in file[1:]: # for each sequence
        outfile.write(x.split('\n')[0].split('\t')[0]+'\t')
        s=''.join(x.split('\n')[1:])
        t=s.replace('-','')
        somme+=len(t)
        outfile.write(str(len(t))+'\n')
    outfile.close()
    print( somme)

if __name__ == '__main__' : 
    if len(sys.argv)==3:List(sys.argv[1],sys.argv[2])
    elif '*' in sys.argv[1]:
        for file in iglob(sys.argv[1]):
            List(file,file+'.length')
    elif len(sys.argv)==2:List(sys.argv[1],sys.argv[1]+'.length')
    else: ask()
else: 
    pass