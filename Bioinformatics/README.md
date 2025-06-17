

# Readme
(C) Siddharth S. Tomar, Arnaud Le Rouzic, Aurélie Hua-Van

## Setting up the development environment for Bioinformatics analysis

Requirements:

- Python (3.10.12)
- R (4.4.3) (ape, phytools ggtree gridExtra ggplot2 ggpubr adephylo castor viridis ggfortify multcomp cluster)
- cutadapt (3.5)
- fqtrim (v0.9.7) [https://ccb.jhu.edu/software/fqtrim/](https://ccb.jhu.edu/software/fqtrim/)
- mafft (v7.490)
- blast+ (2.12.0+)
- RepeatMasker
- fasttree (2.1) 
- bedtools (2.31.1)

Optional:

- sra-toolkit for downloading RNA reads from NCBI
- datasets ([https://github.com/huggingface/datasets]()) for downloading genome assemblies from NCBI

## Copy the Bioinformatics directory

This directory will be the ***Project*** directory.
It contains:

- **Scripts** directory that contains all the scripts to run the analyses and the main figures (R, Bash, Python)

- **Data** directory
	- The phylogenetic tree of the Drosophilini species (BigGene_.tre)
	- The TE classification file (RepBaseDmel.Classif.txt)
	- The multifasta file with the TE consensus (RepBaseDmel.fasta)
	- The file with the taxonomy of Drosophilini species (TaxoDroso.txt)
	- The file with accessions for drosophilidae assemblies (94_Drosophilidae_Accessions.txt)
	- The picluster annotation files in bed format (piCluster.dmel.Strain.r1.bed, etc ...)
			**Sources of annotations**:
		- piCluster.dmel.CanS.r1.bed (Srivastav et al. 2024, restrictive replicate 1)
		- piCluster.dmel.CanS.r2.bed (Srivastav et al. 2024, restrictive replicate 2)
		- piCluster.dmel.OreR.r1.bed (Srivastav et al. 2024, restrictive replicate 1)
		- piCluster.dmel.OreR.r2.bed (Srivastav et al. 2024, restrictive replicate 2)
		- piCluster.dmel.iso1.r0.bed (Brennecke et al. 2007)
		- piCluster.dmel.iso1.r1.bed (Srivastav et al. 2024, restrictive replicate 1)
		- piCluster.dmel.iso1.r2.bed (Srivastav et al. 2024, restrictive replicate 2)

- **Genomes** directory (empty) for the species analysis
- **RepeatMasker** directory (empty) for the species analysis


## Data preparation:

`cd Project`


### small RNAreads: (download, select 23-29 nt-long reads and collapse identical reads)
	- download the fq file for iso1 small RNA reads from ovaries 
		`fasterq-dump SRR11846566 -o Project/Data/SRR11846566.fq`
		
	- download the fq file for CanS small RNA reads from ovaries
		`fasterq-dump SRR5687217 -o Project/Data/SRR5687217.fq`
		
	- download the fq file for OreR small RNA reads from ovaries
		`fasterq-dump SRR25922470 -o Project/Data/SRR25922470.fq`
		
	For each accession: 

		cd Project/Data
		cutadapt -m 23 -M 29 -e 0.1 -q 20 -O 1 -a TGGAATTCTCGG ACCESSION.fq -o ACCESSION.trimmed.fq
		fqtrim -D -C ACCESSION.trimmed.fq -o collapse.fq

	In the Data directory, a file with extension .trimmed.collapse.fq is created
	
###*D. melanogaster* genome assemblies
	 

		cd Project/Data
		for genome in GCA_000001215.4 GCA_003401735.1 GCA_003402015.1
		do
			/opt/datasets download genome accession $genome
			unzip -o ncbi_dataset.zip
			rm ncbi_dataset.zip
			mv ncbi_dataset/data/GCA_000001215.4/* Data
			rm -r ncbi_dataset
		done
		mv Data/GCA_000001215.4 Data/Dmel.iso1.fasta
		mv Data/GCA_003401735.1 Data/Dmel.CanS.fasta
		mv Data/GCA_003402015.1 Data/Dmel.OreR.fasta
		
		
		
		

		



###iso1 
For iso1, the chromosome names are different in the annotation files (Srivastav et al. 2024) and the assemblies. Names have been modified both in the assembly and in the annotation files using a script. The modified annotations are provided. Run the script to change the chromosome names in the newly downloaded genome.


		cd Project
		python Scripts/ChangeChrName.py Data/Dmel.iso1.fasta
		mv Data/Dmel.iso1.fasta.mod Data/Dmel.iso1.fasta

### Generate files containing the lengths of the sequences###
		python3 Scripts/length.py Data/Dmel.iso1.fasta
		python3 Scripts/length.py Data/Dmel.CanS.fasta
		python3 Scripts/length.py Data/Dmel.OreR.fasta

## Run the analyses for *D. melanogaster*

- open the configuration bash script (ConfigDmel2.sh) to set up parameters

- run the script

`cd Project`

`python Scripts/ConfigDmel2.sh`

The script can run successively the analysis on the three assemblies, with different parameters. You can comment or uncomment if you want to rerun one part only

**Parameters**:

- project 
- strain (iso1/CanS/OreR)
- sra (SRR11846566/SRR5687217/SRR25922470)
- pi annotation (r0/r1/r2)
- mismatch
- flavor (erase/keep/All/mapPi/CpPi/...) skip some steps (not fully implemented)
- TE database [Data/RepBaseDmel.fa]

For each assembly, it creates a ResultsDmel.Strain folder and various subfolders.

For each, it calls successively 3 scripts:


### AnalyseDmel.sh (run the analysis) 
- identification of copies within assemblies (independant of the settings)
- determination of copies status based on the piCluster annotation
- mapping of the small RNA reads and analysis 
	
### DefineCategory.R 
(for each TE, analyse the copies tree to classify the TE family according to :)

- recent transposition of young copies
- presence of sequences within piCluster
- presence of old sequences within piCluster
- number of piRNA reads shared between piRNA Cluster copies and euchromatic copies

### Fig4_combined.R 
- create figures 4 for each TE with the copies tree, the shared reads matrix, the divergence of each copy with the consensus


### Run the summary figure (Figure5). 

After analysis of the 3 assemblies with the same parameters is done (Beware, parameters are hard-encoded)

`Rscript Scripts/Fig5_Categories.R`

### List of scripts and embedded scripts:
* ConfigDmel.sh (Launch with the desired parameters)
* AnalyseDmel.sh (run the assembly/smallRNA analysis TE by TE)
* blast2bed.sh (Transform the blast result into bed file)
* TagRed.py (Identify copies within piRNA cluster)
* ChangeName.py (Filter and change the name of the copies)
* AnalyseBed.py (Count the reads per copy after bowtie)
* DropCons.py (Filter and discard the consensus from the alignment)
* CountRead.py (Create the matrix with the read count)

### List of scripts for supplementary figures:
* Fig5_Categories.m3.R (Supp Fig S3, using mismatch parameter = 3)

`$ Rscript Scripts/Fig5_Categories.m3.R`

### list of output files


## Run the analysis for Drosophilini species

### Download the genomes


- run the script: 

`$ python Scripts/ConfigAllSpecies.sh`

The script runs successively the analysis on the all genomes in the Genomes directory.
It runs an analysis on *D. melanogaster* genome as a reference
It creates the main figure (Figure 6) and some supplementary figures on the cross-regulated TEs


### List of scripts and embedded scripts:
* LaunchAllSpecies.sh (launch the analysis for the different genomes)
* AnalyseOneSpecies.sh (Analyse one genome)
* Gff2bed.py (Transform the RM gff file into bed)
* CompileAllSpecies.sh (Gather all results)
* Dmel4AllSpecies.sh (Analyse the reference genome for all genomes)
* MakeTree.py (select the best copies and align for phylogenetic tree for all TEs)
* GetSpecificReads.py (For cross-regulated TEs, analyse the reads compared to the reference)
* AllSpecies_DrawSpecificReadBarplot.R (SupFig S5)
* AllspeciesHeatMap.R (Draw the heatmap figure Figure 6)
* AllSpecies.Data.R (helper script once the analyse is done)
* AllSpecies.functions.R (helper script once the analyse is done)

### Other scripts (Supplementary Figures)
* AllSpecies.Correlation.SupFig.R (SupFig S4 A)
* AllSpecies_SupFigure_SpeciesTree_TEnb.R (SupFig S4 B and C)
* AllSpecies.DrawTEtree.R (SupFig 6 Color the TE trees' tips according to clade)
