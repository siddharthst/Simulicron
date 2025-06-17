#!/bin/bash

### Siddharth S. Tomar
### Aurélie Hua-Van
### Arnaud Le Rouzic

####   The parameters

#1- Set here the project directory
#2- Set the method (gff/blastn) for copy identification (gff: extension for RM) (blast not really implemented)
#3- sra
	#sra=SRR11846566
	#sra=SRR5687217 # the file and path must be in the form $Project/Data/$sra.trimmed.collapse.fq
	#sra=SRR25922470
	#sra=SRR14569563
#4- mismatch
	#mismatch=0
	#mismatch=3
#5- flavor
	# flavor=Cp  do only the Cp count (not implemented)
	# flavor=keep alternate=keep , for not redoing the first part which is independent of the mapping,
	# flavor=All redo everything (except RM)
	# flavor=RM redo everything including RM or blastn (if not keep and not All)
#6- TE database
#7- Species 
    # '' or nothing : All species
    # Z. only Zaprionus (Z.*)
    
####   The scripts

# 1- Scripts/LaunchAllSpecies.sh
	# it starts the analysis for each genome
	# AnalyseOneSpecies
		# (RM)
		# Extracting copies
		# Mapping reads
		
# 2- Scripts/CompileAllSpecies.sh
	#When all genomes are done
		# Run all copies extraction for Dmel
		# Make the phylogenetic trees with the best copies
		
# 3- R script AllSpeciesHeatmap.R 
	# Make some figures

# 4- GetSpecificReads.py 
	# Calculate shared reads between species
	
# 5- AllSpecies_DrawSpecificReadBarplot.R
	# Draw the supplementary figure
				

############## RUN
Project=~/CrossRegulation

method=gff 
sra=SRR11846566
mismatch=3
flavor=All
TEs=Data/RepBaseDmel.fa


SECONDS=0
cd $Project
bash Scripts/LaunchAllSpecies.sh $Project $method $sra $mismatch $flavor $TEs $Species
bash Scripts/CompileAllSpecies.sh $Project $method $sra $mismatch $flavor $TEs 
Rscript Scripts/AllspeciesHeatMap.R $Project $method $sra $mismatch
python3 Scripts/GetSpecificReads.py AllSpecies/Table/Species.TE.$sra.$mismatch.CrossRegulated.txt
Rscript Scripts/AllSpecies_DrawSpecificReadBarplot.R $Project $method $sra $mismatch

ELAPSED="Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo `date` >> Runtime.txt
echo $Project >> Runtime.txt
echo $method >> Runtime.txt
echo $sra >> Runtime.txt
echo $mismatch >> Runtime.txt
echo $flavor >> Runtime.txt
echo $TEs >> Runtime.txt
echo $Species >> Runtime.txt

echo $ELAPSED >> Runtime.txt
echo '#######' >> Runtime.txt

