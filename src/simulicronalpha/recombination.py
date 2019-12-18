import numpy as np
from numpy import random, random_intel
import pandas as pd


def recombination(genomeFrame, genomeOrg):
    # Create insertion list for the progeny
    TEprogeny = []
    TEid = []
    # Create a copy of genomeFrame
    genomeCopy = genomeFrame.copy(deep=True)
    TEid_Father = genomeOrg["TEfather"].copy()
    TEid_Mother = genomeOrg["TEmother"].copy()
    Insertion_Father = genomeOrg["Insertion_Father"]
    Insertion_Mother = genomeOrg["Insertion_Mother"]
    if Insertion_Father[0] == 0 and Insertion_Mother[0] == 0:
        return (
            [0],
            random.choice(
                [genomeOrg["Insertion_Father"], genomeOrg["Insertion_Mother"],]
            ),
        )
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
        if all(v == 0 for v in Insertion_Father):
            pass
        else:
            for i in Insertion_Father:
                se = genomeCopy[genomeCopy["InsertionSite"] == i]["Progeny"].values[0]
                if se == "M":
                    TEprogeny.append(i)
                    TEid.append(TEid_Father.pop(0))
        if all(v == 0 for v in Insertion_Mother):
            pass
        else:
            for i in Insertion_Mother:
                se = genomeCopy[genomeCopy["InsertionSite"] == i]["Progeny"].values[0]
                if se == "F":
                    TEprogeny.append(i)
                    TEid.append(TEid_Mother.pop(0))
        # print(genomeCopy['Progeny'].value_counts())
    if not TEprogeny:
        TEprogeny.append(0)
    if not TEid:
        TEid.append(0)
    return (TEprogeny, TEid)
