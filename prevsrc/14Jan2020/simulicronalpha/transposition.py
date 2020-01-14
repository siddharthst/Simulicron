import numpy as np
from numpy import random, random_intel
from coolname import generate_slug
import pandas as pd
import uuid


def transposition(self, genomeFrame=None, transposonFrame=None, orgFrame=None):
    # Skip the check for empty transposon content. Do it in the parent
    TEfather = 
    TEmother = 
    filledSites = list(
        filter((0).__ne__, orgFrame.iloc[0]["TEfather"] + orgFrame.iloc[0]["TEmother"])
    )
    GenomicSites = pd.Series(
        genomeFrame.InsertionProbability.values, index=genomeFrame.InsertionSiteID,
    )    
    GenomicSites.drop(filledSites, axis=1, inplace=True)
    for i in filledSites:
        transpositionRate = transposonFrame.loc[transposonFrame["TID"] == i, "TraRate"]
        if transpositionRate > np.random_intel.uniform(0, 1.0):
            insertionSite = np.random_intel.choice(
                GenomicSites.index.values, 1, p=GenomicSites.InsertionProbability.values
            )
            newParent = random.choice(["M", "F"])
            GenomicSites.drop(insertionSite, axis=1, inplace=True)
            TID = uuid.uuid4().hex
            row = pd.Series(
                {
                    "TID": TID,
                    "InsertionSite": insertionSite,
                    "TraRate": transpositionRate,
                    "Name": generate_slug(),
                    "Class": transposonFrame.loc[transposonFrame["TID"] == i, "Class"],
                    "Traceback": [i],
                    "Generation": 1,
                    "Parent": newParent,
                }
            )
            if newParent == "M":
                TEfather.append(TID)
            else:
                TEmother.append(TID)
            transposonFrame = transposonFrame.append(row, ignore_index=True)

    return 0
