#!/bin/bash

### Siddharth S. Tomar
### Aurélie Hua-Van
### Arnaud Le Rouzic

###########################################################################################
####   The parameters
###########################################################################################

#1- Set here the project directory
#2- Set the strain (iso1, CanS, OreR)
#3- sra
	#sra=SRR11846566
	#sra=SRR5687217 # the file and path must be in the form $Project/Data/$sra.trimmed.collapse.fq
	#sra=SRR25922470
	#sra=SRR14569563
#4- Set the pi annotation (r0 to r4) 
#5- mismatch
	#mismatch=0
	#mismatch=3
#6- flavor
	# flavor=Cp,  do only the Cp count
	# flavor=CpPi,  do only the Cp count and the Cp annotation
	# flavor=keep , skip the blast and sequence extraction, redo the copy annotation and the mapping
	# flavor=All , redo everything
	# flavor=mapPi , do the mapping only
	# flavor=erase , erase everything
#7- TE database 
	# TEs=Data/RepBaseDmel.fa
    
###########################################################################################
####   The scripts
###########################################################################################

# 1- Scripts/AnalyseDmel.sh
	# Blast the genome
	# For each TE
		# Extract copies 
		# Identify pi copies
		# Filter out euchromatic copies <70% length of the consensus
		# Align copies to the consensus, deleting any insertions
		# Build the tree
		# Map the reads on the copies
		# Calculate the number of read for each copy
		# Save a matrix
	# Concat all results into two tables
		# 1- table with some data for each TE family (1 family per row)
		# 2- table with some data for each copy (1 copy per row)

		
# 2- Scripts/MakeUltrametric.3.R
	# Create Ultrametric tree from the tree of each TE
	# Determine the status of each TE family regarding presence absence of old pi copy and their regulatory potential
	# Determine if recent burst or not (old burst or not burst)
	# Determine if no/one/several waves of amplification (not really used)
	
# 3- Scripts/Fig4.R 
	# Draw the Figure 5 for each TE (should be renamed!)
	
# Create the log file computing time (uncomment if you want to erase previous data)
#> Dmel.Runtime.txt

###########################################################################################
### 	The different runs (the 3 first correspond to the 3 genomes and corresponding pi annotation, and mismatch 0)
###########################################################################################

############# RUN iso1
# iso1 r0 0 mismatch the shortest, on the server
Project=/home/huavan/CrossRegulation
strain=iso1
pi=r0
sra=SRR11846566
mismatch=0
flavor=mapPi
TEs=Data/RepBaseDmel.fa

SECONDS=0
cd $Project
#bash Scripts/AnalyseDmel.sh $Project $strain $sra $pi $mismatch $flavor $TEs 
Rscript Scripts/DefineCategory.R  $sra $mismatch $pi $strain
#Rscript Scripts/Fig4_Combined.R  $sra $mismatch $pi $strain

#wait
ELAPSED="Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo `date`>> Dmel.Runtime.txt
echo Run1_0 >> Dmel.Runtime.txt
echo $Project $strain $sra $pi $mismatch $flavor $TEs >> Dmel.Runtime.txt
echo $ELAPSED >> Dmel.Runtime.txt 
echo '#######'>> Dmel.Runtime.txt


############# Run OreR 
Project=/home/huavan/CrossRegulation
strain=OreR
pi=r1
sra=SRR25922470
mismatch=0
flavor=mapPi
TEs=Data/RepBaseDmel.fa

SECONDS=0
cd $Project
#bash Scripts/AnalyseDmel.sh $Project $strain $sra $pi $mismatch $flavor $TEs 
Rscript Scripts/DefineCategory.R  $sra $mismatch $pi $strain
#Rscript Scripts/Fig4_Combined.R  $sra $mismatch $pi $strain

ELAPSED="Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"

echo `date`>> Dmel.Runtime.txt
echo Run_OreR_0 >> Dmel.Runtime.txt
echo $Project $strain $sra $pi $mismatch $flavor $TEs >> Dmel.Runtime.txt
echo $ELAPSED >> Dmel.Runtime.txt 
echo '#######'>> Dmel.Runtime.txt


############# Run CanS the shortest, on the server
Project=/home/huavan/CrossRegulation
strain=CanS
pi=r1
sra=SRR5687217
mismatch=0
flavor=mapPi
TEs=Data/RepBaseDmel.fa

SECONDS=0
cd $Project
#bash Scripts/AnalyseDmel.sh $Project $strain $sra $pi $mismatch $flavor $TEs 
Rscript Scripts/DefineCategory.R  $sra $mismatch $pi $strain
#Rscript Scripts/Fig4_Combined.R  $sra $mismatch $pi $strain
wait
ELAPSED="Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
# 
echo `date`>> Dmel.Runtime.txt
echo Run_CanS_0 >> Dmel.Runtime.txt
echo $Project $strain $sra $pi $mismatch $flavor $TEs >> Dmel.Runtime.txt
echo $ELAPSED >> Dmel.Runtime.txt 
echo '#######'>> Dmel.Runtime.txt

Rscript Scripts/Fig5_Categories.R

# ##### Run iso1 with canS pi the shortest, on the server
# Project=/mnt/Poles/Genomes/huavan/CrossRegulation
# strain=iso1
# pi=r1
# sra=SRR25922470
# mismatch=0
# flavor=keep
# TEs=Data/RepBaseDmel.fa
# 
# SECONDS=0
# cd $Project
# #bash Scripts/AnalyseDmel.sh $Project $strain $sra $pi $mismatch $flavor $TEs 
# #Rscript Scripts/MakeUltrametric.R  $sra $mismatch $pi $strain
# #Rscript Scripts/Fig4.R  $sra $mismatch $pi $strain
# 
# ELAPSED="Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
# 
# echo `date`>> Dmel.Runtime.txt
# echo Run_Iso1_OreR_0 >> Dmel.Runtime.txt
# echo $Project $strain $sra $pi $mismatch $flavor $TEs >> Dmel.Runtime.txt
# 
# echo $ELAPSED >> Dmel.Runtime.txt 
# echo '#######'>> Dmel.Runtime.txt
# 
# 
# # # Run iso1 with OreR pi the shortest, on the server
# Project=/mnt/Poles/Genomes/huavan/CrossRegulation
# strain=iso1
# pi=r1
# sra=SRR5687217
# mismatch=0
# flavor=keepPi
# TEs=Data/RepBaseDmel.fa
# 
# SECONDS=0
# cd $Project
# #bash Scripts/AnalyseDmel.sh $Project $strain $sra $pi $mismatch $flavor $TEs 
# #Rscript Scripts/MakeUltrametric.R  $sra $mismatch $pi $strain
# #Rscript Scripts/Fig4.R  $sra $mismatch $pi $strain
# 
# ELAPSED="Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
# 
# echo `date`>> Dmel.Runtime.txt
# echo Run_Iso1_CanS_0 >> Dmel.Runtime.txt
# echo $Project $strain $sra $pi $mismatch $flavor $TEs >> Dmel.Runtime.txt
# 
# 
# echo $ELAPSED >> Dmel.Runtime.txt 
# echo '#######'>> Dmel.Runtime.txt
# # # iso1 r0 3 mismatch the shortest, on the server
# Project=/mnt/Poles/Genomes/huavan/CrossRegulation
# strain=iso1
# pi=r0
# sra=SRR11846566
# mismatch=3
# flavor=keepPi
# TEs=Data/RepBaseDmel.fa
# 
# SECONDS=0
# cd $Project
# #bash Scripts/AnalyseDmel.sh $Project $strain $sra $pi $mismatch $flavor $TEs 
# Rscript Scripts/MakeUltrametric.R  $sra $mismatch $pi $strain
# #Rscript Scripts/Fig4.R  $sra $mismatch $pi $strain
# 
# ELAPSED="Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
# 
# echo `date`>> Dmel.Runtime.txt
# echo Run2_0 >> Dmel.Runtime.txt
# echo $Project $strain $sra $pi $mismatch $flavor $TEs >> Dmel.Runtime.txt
# 
# 
# echo $ELAPSED >> Dmel.Runtime.txt 
# echo '#######'>> Dmel.Runtime.txt
# 
# # # Run3 iso1 r1 0 mismatch the shortest, on the server
# Project=/mnt/Poles/Genomes/huavan/CrossRegulation
# strain=iso1
# pi=r1
# sra=SRR11846566
# mismatch=0
# flavor=keep
# TEs=Data/RepBaseDmel.fa
# 
# SECONDS=0
# cd $Project
# #bash Scripts/AnalyseDmel.sh $Project $strain $sra $pi $mismatch $flavor $TEs 
# #Rscript Scripts/MakeUltrametric.3.R  $sra $mismatch $pi $strain
# #Rscript Scripts/Fig4.R  $sra $mismatch $pi $strain
# 
# ELAPSED="Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
# 
# echo `date`>> Dmel.Runtime.txt
# echo Run3_0 >> Dmel.Runtime.txt
# echo $Project $strain $sra $pi $mismatch $flavor $TEs >> Dmel.Runtime.txt
# 
# echo $ELAPSED >> Dmel.Runtime.txt 
# echo '#######'>> Dmel.Runtime.txt
# 
# # # Run4 iso1 r1 3 mismatch the shortest, on the server
# Project=/mnt/Poles/Genomes/huavan/CrossRegulation
# strain=iso1
# pi=r1
# sra=SRR11846566
# mismatch=3
# flavor=keepPi
# TEs=Data/RepBaseDmel.fa
# 
# SECONDS=0
# cd $Project
# #bash Scripts/AnalyseDmel.sh $Project $strain $sra $pi $mismatch $flavor $TEs 
# #Rscript Scripts/MakeUltrametric.R  $sra $mismatch $pi $strain
# Rscript Scripts/Fig4.R  $sra $mismatch $pi $strain
# 
# ELAPSED="Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
# 
# echo `date`>> Dmel.Runtime.txt
# echo Run4_0 >> Dmel.Runtime.txt
# echo $Project $strain $sra $pi $mismatch $flavor $TEs >> Dmel.Runtime.txt
# 
# 
# echo $ELAPSED >> Dmel.Runtime.txt 
# echo '#######'>> Dmel.Runtime.txt
# 
# # # Run5 iso1 r3 0 mismatchthe shortest, on the server
# Project=/mnt/Poles/Genomes/huavan/CrossRegulation
# strain=iso1
# pi=r3
# sra=SRR11846566
# mismatch=0
# flavor=keep
# TEs=Data/RepBaseDmel.fa
# 
# SECONDS=0
# cd $Project
# #bash Scripts/AnalyseDmel.sh $Project $strain $sra $pi $mismatch $flavor $TEs 
# #Rscript Scripts/MakeUltrametric.R  $sra $mismatch $pi $strain
# #Rscript Scripts/Fig4.R  $sra $mismatch $pi $strain
# 
# ELAPSED="Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
# 
# echo `date`>> Dmel.Runtime.txt
# echo Run5_0 >> Dmel.Runtime.txt
# echo $Project $strain $sra $pi $mismatch $flavor $TEs >> Dmel.Runtime.txt
# 
# echo $ELAPSED >> Dmel.Runtime.txt 
# echo '#######'>> Dmel.Runtime.txt
# 
# # # Run5 iso1 r3 3 mismatch the shortest, on the server
# Project=/mnt/Poles/Genomes/huavan/CrossRegulation
# strain=iso1
# pi=r3
# sra=SRR11846566
# mismatch=3
# flavor=keepPi
# TEs=Data/RepBaseDmel.fa
# 
# SECONDS=0
# cd $Project
# #bash Scripts/AnalyseDmel.sh $Project $strain $sra $pi $mismatch $flavor $TEs 
# #Rscript Scripts/MakeUltrametric.R  $sra $mismatch $pi $strain
# Rscript Scripts/Fig4.R  $sra $mismatch $pi $strain
# 
# ELAPSED="Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
# 
# echo `date`>> Dmel.Runtime.txt
# echo Run6_0 >> Dmel.Runtime.txt
# echo $Project $strain $sra $pi $mismatch $flavor $TEs >> Dmel.Runtime.txt
# 
# echo $ELAPSED >> Dmel.Runtime.txt 
# echo '#######'>> Dmel.Runtime.txt