#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys,os, glob

def GFF2bed(gfffile, out):
	print(">>> Gff2Bed")
	with open(gfffile ,"r") as fh:
		tab=[]
		for line in fh:
			if '#' in line: continue
			TE=line.split('Motif:')[1].split('"')[0] 
			qstart=line.split(" ")[-2]
			qend=line.split(" ")[-1].rstrip()
			line=line.rstrip().split("\t")
			tab.append('\t'.join([line[0],line[3],line[4],TE,line[5].strip(),line[6], qstart,qend]))
			#tab.append(line[0]+"\t"+line[3]+"\t"+line[4]+"\t"+TE+"\t"+line[5].strip()+"\t"+line[6])
	with open(out+os.path.basename(gfffile)+'.bed' ,"w") as bedfile:
		#print(os.path.expanduser('ResultsDro/Bed/'+os.path.basename(gfffile)+'.bed'))
		bedfile.write('\n'.join(tab))

#print(sys.argv[1])
for file in sys.argv[1:-1]:
    print(file)
    GFF2bed(file,sys.argv[-1])
