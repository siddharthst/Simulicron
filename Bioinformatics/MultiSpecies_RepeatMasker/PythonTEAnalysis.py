import os
import pandas as pd
from functools import reduce

# Create list of dataFrames
SAMlist                 = []
# Create a list of result files
resultDir               = './Alignment//'        
pattern                 = '.sam'
Result_files            = [f for f in os.listdir(resultDir) if pattern in f]
# Open a sam file 
for SAMfile in Result_files:
    SpeciesName                 = SAMfile.split('/')[-1].split('.fasta.fa.msa.msa.cor.sam')[0]
    try:
        SAMframe                = pd.read_csv('./Alignment/' + SAMfile, sep="\t", comment='@', header=None)
        SAMDict                 = SAMframe.groupby(2)[0].apply(list).to_dict()
        SAMframe                = pd.DataFrame({"TE":list(SAMDict.keys()), "Reads":list(SAMDict.values())})
        SAMframe['TE']          = SAMframe['TE'].str.split('::').str[0]
        SAMframe[SpeciesName]   = SAMframe['Reads'].str.len()
        SAMframe                = SAMframe.drop('Reads', axis=1)
        # Choose the best TE copy in terms of hits
        SAMframe                = SAMframe.sort_values(SpeciesName, ascending=False).drop_duplicates(['TE'])
        SAMframe                = SAMframe.set_index('TE')
        SAMlist.append(SAMframe)
    except:
        print(SpeciesName, " No SAM records and has been skipped.")


# Merge dataframes on index
MergedDataFrame = SAMlist[0]
for df in SAMlist[1:]:
    MergedDataFrame = MergedDataFrame.join(df, how='outer')
# Transpose dataframe
MergedDataFrame = MergedDataFrame.fillna(0)
MergedDataFrame = MergedDataFrame.transpose()
# Save the dataframe
MergedDataFrame.to_csv("./Results/CoreSet.tsv", sep="\t")




