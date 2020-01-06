import numpy as np
from numpy import random, random_intel
import pandas as pd


def recombination(genomeFrame, genomeOrg, transposonFrame):
    # Create insertion list for the progeny
    TEid = []
    # Create a copy of genomeFrame
    genomeCopy = genomeFrame.copy(deep=True)
    # Added lookup to transposon frame
    # Instead of looking into population database the insertion
    # sites are now accessed directly from the transposon data-
    # base
    TEid_Father = genomeOrg["TEfather"].tolist()
    TEid_Mother = genomeOrg["TEmother"].tolist()

    if TEid_Father[0] == 0 and TEid_Mother[0] == 0:
        return [0]
    else:
        initParent = random.choice(["M", "F"])
        surrogateParent = "M" if ("M" != initParent) else "F"
        RandArray = np.random_intel.uniform(0, 1.0, genomeFrame.shape[0])
        switch = genomeCopy["RecombinationRate"] > RandArray
        # counter = collections.Counter(switch)
        # print(counter)
        genomeCopy["Progeny"] = np.where(
            switch.cumsum() % 2 == 0, initParent, surrogateParent
        )
        if all(v == 0 for v in TEid_Father):
            pass
        else:
            for i in TEid_Father:
                insertionSite = transposonFrame.loc[transposonFrame["TID"] == i][
                    "InsertionSite"
                ].values[0]
                se = genomeCopy[genomeCopy["InsertionSite"] == insertionSite][
                    "Progeny"
                ].values[0]
                if se == "M":
                    TEid.append(i)

        if all(v == 0 for v in TEid_Mother):
            pass
        else:
            for i in TEid_Mother:
                insertionSite = transposonFrame.loc[transposonFrame["TID"] == i][
                    "InsertionSite"
                ].values[0]
                se = genomeCopy[genomeCopy["InsertionSite"] == insertionSite][
                    "Progeny"
                ].values[0]
                if se == "F":
                    TEid.append(i)

        # print(genomeCopy['Progeny'].value_counts())
    if not TEid:
        TEid.append(0)
    return (TEid)
