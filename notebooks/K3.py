# Init
import os
import sys

module_path = os.path.abspath(os.path.join("../src/simulicronalpha/"))
if module_path not in sys.path:
    sys.path.append(module_path)

# Imports
import random
import numpy as np
import pandas as pd
import warnings
import pickle

from numpy import concatenate as c

from popSim import runBatch


# 100
sim   = runBatch(    
    numberOfSimulations=8,
    baseSelection=0,
    baseInsertionProb=1,
    numberOfInsertionSites=10000,
    numberOfChromosomes=6,
    baseRecombinationRate=0.01,
    baseTau=1,
    numberOfPiRNA=6,
    piPercentage=3,
    enablePiRecombination=False,
    NumberOfIndividual=1000,
    NumberOfTransposonTypes=2,
    NumberOfInsertionsPerType=[1,1],
    FrequencyOfInsertions=[0.5,0.5],
    ExcisionRates=[0.2,0.5],
    RepairRates=[1,1],
    InsertionRates=[1,1],
    HardyWeinberg=False,
    NumberOfGenerations=10000,
    numberOfThreads=8,
)
with open("K3U0205.pickle", "wb") as f:
    pickle.dump((sim), f)
