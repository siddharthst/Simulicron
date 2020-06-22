import numpy as np
from numpy import cumsum
from numpy import concatenate as c
from regulation import regulation
import random


def transposition(
    transposonMatrix, genomeMatrix, NumberOfTransposonInsertions, TEset, v1, v2,
):
            # Code to track the genealogy of TE
            for k in range(NumberOfTransposonInsertions):
                if transposonsToTranspose[i] in TEset[k + 1]:
                    TEset[k + 1].add(len(transposonMatrix) - 1)
                    pass

    return (allele1Index, allele2Index, transposonMatrix, TEset)
