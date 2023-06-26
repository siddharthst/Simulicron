# coding: utf-8

# In[1]:


import os
import pandas as pd
from functools import reduce


# In[2]:


resultDir = './Results/'        
pattern = 'All_Tables'
Result_files = [f for f in os.listdir(resultDir) if pattern in f]


# In[3]:


ResultDict = {}
for files in Result_files:
    with open(resultDir + files) as testFile:
        if 'Superfamily' in testFile.read():            
            species = files.replace('All_Tables_', '').replace('.fasta.out.bed','')
            dataFrame = pd.read_csv(resultDir + files, sep="\t", index_col=0, usecols=[0,1])
            dataFrame.columns = [species]
            ResultDict[species] = dataFrame    


# In[4]:


AllDataFrames = list(ResultDict.values())
# df_merge = reduce(lambda df_x, df_y: pd.join(df_x, df_y, how="outer"), AllDataFrames)
MergedDataFrame = AllDataFrames[0]
for df in AllDataFrames[1:]:
    MergedDataFrame = MergedDataFrame.join(df, how='outer')
MergedDataFrame = MergedDataFrame.fillna(0)
LogicalDataframe = MergedDataFrame.copy(deep=True)
LogicalDataframe[LogicalDataframe != 0] = 1
LogicalDataframe['Total'] = LogicalDataframe.sum(axis=1)
LogicalDataframe = LogicalDataframe.sort_values('Total', ascending=False)
ShortListedDataFrame = LogicalDataframe.head(50)
ShortListedDataFrame = ShortListedDataFrame[['Total']]


# In[11]:


# Save the shortlisted TEs
with open('shortlisted_TEs.txt', mode='wt', encoding='utf-8') as myfile:
    myfile.write('\n'.join(list(ShortListedDataFrame.index)
))
# Save all data
LogicalDataframe.to_csv('Extended.tsv', sep="\t")
MergedDataFrame.to_csv('Logical.tsv', sep="\t")


# In[ ]:




