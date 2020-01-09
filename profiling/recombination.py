#!/usr/bin/env python
# coding: utf-8

# In[1]:


from cProfile import Profile
prof = Profile()
prof.disable()
import os
import sys
module_path = os.path.abspath(os.path.join('../src/simulicronalpha/'))
if module_path not in sys.path:
    sys.path.append(module_path)
# Imports
import inspect
import random
import tempfile
import uuid
from functools import wraps
import collections
import seaborn as sns
import matplotlib.pyplot as plt
import coloredlogs
import logging
import numpy as np
from numpy import random_intel
import pandas as pd
from coolname import generate_slug
import warnings
import cProfile

warnings.filterwarnings("ignore")

# Logger
logger = logging.getLogger(__name__)
coloredlogs.install(level="DEBUG")


# In[2]:


# Generate sample dataset for the function:
nSites = 1000000
nTransposons = 10
TE = [uuid.uuid4().hex for _ in range(nTransposons)]
TEinsertionSites = random.sample(range(nSites), nTransposons)
TEFrame = pd.DataFrame()
TEFrame['TID'] = TE
TEFrame['InsertionSite'] = TEinsertionSites
GenFrame = pd.DataFrame()
GenFrame['InsertionSiteID'] = list(range(nSites))
GenFrame['RecombinationRate'] = np.random.uniform(0.01, 0.99, nSites)
genomeOrg = pd.DataFrame()
genomeOrg['TEfather'] = [TE] * 10
genomeOrg['TEmother'] = [[0] * nTransposons]*10


# In[3]:


#prof.enable()
@profile
def recombination(TranspFrame, GenFrame, genomeOrg=None):
    # Create insertion list for the progeny
    TEid = []
    # Create a copy of genomeFrame
    genomeCopy = GenFrame.copy(deep=True)
    # Added lookup to transposon frame
    # Instead of looking into population database the insertion
    # sites are now accessed directly from the transposon data-
    # base
    TEid_Father = genomeOrg["TEfather"]
    TEid_Mother = genomeOrg["TEmother"]

    if TEid_Father[0] == 0 and TEid_Mother[0] == 0:
        return [0]
    else:
        initParent = random.choice(["M", "F"])
        surrogateParent = "M" if ("M" != initParent) else "F"
        RandArray = np.random_intel.uniform(0, 1.0, GenFrame.shape[0])
        switch = genomeCopy["RecombinationRate"] > RandArray
        genomeCopy["Progeny"] = np.where(
            switch.cumsum() % 2 == 0, initParent, surrogateParent
        )
        # Blank statement
        z = np.where(
            switch.cumsum() % 2 == 0, initParent, surrogateParent
        )
        if all(v == 0 for v in TEid_Father):
            pass
        else:
            for i in TEid_Father:
                insertionSite = TranspFrame.loc[TranspFrame["TID"] == i][
                    "InsertionSite"
                ].values[0]
                se = genomeCopy[genomeCopy["InsertionSiteID"] == insertionSite][
                    "Progeny"
                ].values[0]
                if se == "M":
                    TEid.append(i)
                    
        if all(v == 0 for v in TEid_Mother):
            pass
        else:
            for i in TEid_Mother:
                insertionSite = TranspFrame.loc[TranspFrame["TID"] == i][
                    "InsertionSite"
                ].values[0]
                se = genomeCopy[genomeCopy["InsertionSiteID"] == insertionSite][
                    "Progeny"
                ].values[0]
                if se == "F":
                    TEid.append(i)

    if not TEid:
        TEid.append(0)
    return TEid

#-------------------------#
for i in list(range(100)):
    for k in list(range(100)):
        z = recombination(TEFrame, GenFrame, genomeOrg.iloc[0])
#prof.disable() 
#prof.dump_stats('recombination.stats')


# In[ ]:




