import numpy as np
from numpy import random, random_intel
from coolname import generate_slug
import pandas as pd
import uuid


def transposition(TEfather, TEmother, genomeFrame, transposonFrame, generation):
    # Skip the check for empty transposon content. Do it in the parent
    Insertion_Father = []
    Insertion_Mother = []
    GenomicSites = pd.Series(
        genomeFrame.InsertionProbability.values, index=genomeFrame.InsertionSiteID,
    )
    # GenomicSites = GenomicSites[~GenomicSites.index.isin(TEfather + TEmother)]
    GenomeSites.drop(TEfather + TEmother, axis=1, inplace=True)
    for i in TEfather + TEmother:
        transpositionRate = transposonFrame.loc[transposonFrame["TID"] == i, "TraRate"]
        if transpositionRate > np.random_intel.uniform(0, 1.0):
            insertionSite = np.random_intel.choice(
                GenomicSites.index.values, 1, p=GenomicSites.InsertionProbability.values
            )
            newParent = random.choice(["M", "F"])
            GenomeSites.drop(insertionSite, axis=1, inplace=True)
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
            if (newParent == 'M'):
                TEfather.append(TID)
            else:
                TEmother.append(TID)
            transposonFrame = transposonFrame.append(row, ignore_index=True)
        
    return 0
