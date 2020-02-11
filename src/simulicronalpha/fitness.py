import numpy as np
from numpy import cumsum
from numpy import concatenate as c
import random

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
