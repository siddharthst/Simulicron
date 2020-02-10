import numpy as np
from numpy import cumsum
from numpy import concatenate as c
import random

def generateFitness(
    populationMatrix, transposonMatrix, fitnessFunction=1
):
    allele1Index = np.nonzero(populationMatrix[0:, 0])[0]
    allele2Index = np.nonzero(populationMatrix[0:, 1])[0]
    for i in np.concatenate((allele1Index, allele2Index), axis=0):
        transposonContent = np.trim_zeros(
            np.hstack(
                [populationMatrix[i][0]] + [populationMatrix[i][1]]
            )
        )
        SelectionCoef = transposonMatrix[transposonContent, 2]
        if fitnessFunction == 1:
            populationMatrix[i][2] = np.exp(sum(SelectionCoef))
    return populationMatrix[:, 2]


def calculateFitness(
    transposonMatrix, v1, v2, fitnessFunction=1,
):
    cV1 = v1
    cV2 = v2
    if cV1 == 0:
        cV1 = np.asarray([])
    else:
        cV1 = np.asarray(v1)
    if cV2 == 0:
        cV2 = np.asarray([])
    else:
        cV2 = np.asarray(v2)
    teContent = c([cV1, cV2]).astype(int)
    penalties = transposonMatrix[teContent, 2]
    if fitnessFunction == 1:
        return np.exp(sum(penalties))
