from cProfile import Profile
prof = Profile()
prof.disable()
import numpy as np
from numpy.random import uniform
from numpy import cumsum
from numpy import concatenate as c
import pandas as pd

NumberOfSites = 1000
v1 = np.full(NumberOfSites, "M")
v2 = np.full(NumberOfSites, "F")
rates = uniform(size=(NumberOfSites - 1))

def genTup(seq, n):
    return [tuple(seq[max(i, 0):i + n]) for i in range(-n+1, len(seq)-1)]

def recombination1(v1, v2, rates):
    assert len(v1) == len(v2) and len(v1) == (len(rates) + 1), "Length mismatch"
    rec = uniform(size=len(rates)) < rates
    start = [2 if (uniform() < 0.5) else 1]
    RecombPoints = np.unique(c([[0], np.where(rec == True)[0], [len(rates)]]))
    RecombValues = np.empty((len(RecombPoints),))
    RecombValues[::2] = start[0]
    RecombValues[1::2] = [1 if start[0] == 2 else 2]
    it = iter(RecombPoints)
    RecombRanges = list(zip(it, it))
    #RecombRanges = genTup(RecombPoints,2)
    return (dict(zip(RecombRanges, RecombValues)))

for i in list(range(10000)):
    recombination1(v1, v2, rates)
