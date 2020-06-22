import numpy as np
from numpy import cumsum
from numpy import concatenate as c
import random

def recombination(rates, transposonMatrix, v1, v2):
    rec = np.random.uniform(size=len(rates)) < rates
    start = [0 if (np.random.uniform() < 0.5) else 1]
    whichcol = 1 + cumsum(c((start, rec))) % 2
    allele1 = set(np.where(whichcol == 1)[0])
    allele2 = set(np.where(whichcol == 2)[0])
    transposon1 = transposonMatrix[[v1], 1].tolist()[0]
    transposon2 = transposonMatrix[[v2], 1].tolist()[0]
    haplotype1 = [x in allele1 for x in transposon1]
    haplotype2 = [x in allele2 for x in transposon2]
    r1 = v1[haplotype1]
    r2 = v2[haplotype2]
    if r1.size == 1:
        if r1 == [0]:
            r1 = np.asarray([])
    if r2.size == 1:
        if r2 == [0]:
            r2 = np.asarray([])
    if r1.size == 0 and r2.size == 0:
        return 0
    elif r1.size == 0 and r2.size != 0:
        return list(r2)
    elif r1.size != 0 and r2.size == 0:
        return list(r1)
    else:
        return list(c([r1, r2]))

def recombination(rates, transposonMatrix, v1, v2)