#!/usr/bin/env python


import sys, os

""" 
After blasting and merging, the general bed file indicate which families have been merged
The script uses the number of different families to tag the sequence
Input: The merged bedfile from the general blast (all TEs together)
Output: The filtered merged bedfile from the general blast keeping only lines having map on several TEs
"""


def TagRedundancy(mergedbedfile):

    overlap=open(mergedbedfile.split('.merged.bed')[0]+'.overlap.bed','w')
    
    bed=open(mergedbedfile).readlines()
    for line in bed:
        columns=line.split('\t')
        TEs=line.split('\t')[3].split(',')
        if ',' in columns[3]: 
            columns[3]='*'
            overlap.write('\t'.join(columns))
    overlap.close()
        
TagRedundancy(sys.argv[1])
        
