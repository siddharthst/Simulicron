import numpy as np
from numpy import cumsum
from numpy import concatenate as c
import random


def recombination(rates, transposonMatrix, v1, v2):
    # Empty vectors to store result
    r1 = []
    r2 = []
    # Creating lambda (macro)
    # "Match" does not exist in python
    # match = lambda a, b: [ b.index(x) if x in b else 0 for x in a ]
    # Get the postion of transposons
    positionV1 = transposonMatrix[v1, 1].astype(int).tolist()
    positionV2 = transposonMatrix[v2, 1].astype(int).tolist()
    # This step sorts the locations and adds another location,
    # 0 if not already present
    unqiquePos = list(set(positionV1 + positionV2) | set([0]))
    unqiquePos.sort()
    # Calculate the effective rate from genome map
    effectiveRates = 0.5 * (1 - np.exp(-2 * np.diff(rates[unqiquePos])))
    # print (effectiveRates)
    # Performing "Recombination"
    rec = np.random.uniform(size=len(effectiveRates)) < effectiveRates
    # Select the direction to start from
    start = [(np.random.uniform() < 0.5)]
    # Concat. start and recombination
    # Also remove the added 0 and start from
    # whichhaplo and uniqpos
    whichhaplo = 1 + cumsum(c((start, rec))) % 2
    whichhaplo = np.delete(whichhaplo, 0)
    del unqiquePos[0]
    unqiquePos = np.asarray(unqiquePos)
    # Generating the haplotype
    # Also checking if there is no transposon left
    # In the case above, return a (int) 0
    # Else return the array containing transposons
    if positionV1 == [0]:
        r1 = []
    else:
        if not any(whichhaplo == 1):
            pass
        else:
            pos = set(unqiquePos[whichhaplo == 1])
            for i in v1:
                if (transposonMatrix[i, 1]) in pos:
                    r1.append(i)
    if positionV2 == [0]:
        r2 = []
    else:
        if not any(whichhaplo == 2):
            pass
        else:
            pos = set(unqiquePos[whichhaplo == 2])
            for i in v2:
                if (transposonMatrix[i, 1]) in pos:
                    r2.append(i)
    # Merge to create gamate
    r = r1 + r2
    # Return 0 if no transposon remains
    if not r:
        return 0
    else:
        return r
