#!/usr/bin/env python
# coding: utf-8

from cProfile import Profile
prof = Profile()
prof.disable()
import numpy as np
from _Cyrecombination_ import recombination


NumberOfSites = 1000
v1 = np.full(NumberOfSites, "M")
v2 = np.full(NumberOfSites, "F")
rates = np.random.uniform(size=(NumberOfSites - 1))

prof.enable()
for i in list(range(10000)):
    recombination(v1, v2, rates)
prof.disable()
prof.dump_stats('recombination_cython.stats')



