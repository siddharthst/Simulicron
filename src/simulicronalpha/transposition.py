import numpy as np
from numpy import random, random_intel
import pandas as pd
import uuid


def transposition(TEfather, TEmother, genomeFrame, transposonFrame, generation):
    # Skip the check for empty transposon content. Do it in the main
    GenomicSites = pd.Series(
        genomeFrame.InsertionProbability.values, index=genomeFrame.InsertionSiteID,
    )
    GenomicSites = GenomicSites[~GenomicSites.index.isin(TEfather + TEmother)]
    for i in TEfather + TEmother:
        transpositionRate = transposonFrame.loc[transposonFrame["TID"] == i, "TraRate"]
        if transpositionRate > np.random_intel.uniform(0, 1.0):
            insertionSite = np.choice(
                list_of_candidates, number_of_items_to_pick, p=probability_distribution
            )
    return 0
