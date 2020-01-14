#!/usr/bin/env python
# coding: utf-8

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

prof.enable()

def recombination(v1, v2, rates):
    assert len(v1) == len(v2) and len(v1) == (
        len(rates) + 1), "Length mismatch"
    rec = uniform(size=len(rates)) < rates
    start = [0 if (uniform() < 0.5) else 1]
    whichcol = 1 + cumsum(c((start, rec))) % 2
    return (np.where(whichcol == 1, v1, v2))


for i in list(range(10000)):
    recombination(v1, v2, rates)
prof.disable()
prof.dump_stats('recombination_numpy_only.stats')



