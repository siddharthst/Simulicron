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

'''
# 100
sim   = runBatch(    
    numberOfSimulations=20,
    baseSelection=0,
    baseInsertionProb=1,
    numberOfInsertionSites=10000,
    numberOfChromosomes=6,
    baseRecombinationRate=0.01,
    baseTau=1,
    numberOfPiRNA=6,
    piPercentage=3,
    enablePiRecombination=False,
    NumberOfIndividual=100,
    NumberOfTransposonTypes=1,
    NumberOfInsertionsPerType=[1],
    FrequencyOfInsertions=[1.0],
    ExcisionRates=[1.0],
    RepairRates=[1],
    InsertionRates=[1],
    HardyWeinberg=False,
    NumberOfGenerations=1000,
    numberOfThreads=8,
)
with open("K1.pickle", "wb") as f:
    pickle.dump((sim), f)
'''

# 100
sim   = runBatch(    
    numberOfSimulations=20,
    baseSelection=0,
    baseInsertionProb=1,
    numberOfInsertionSites=10000,
    numberOfChromosomes=6,
    baseRecombinationRate=0.01,
    baseTau=1,
    numberOfPiRNA=6,
    piPercentage=3,
    enablePiRecombination=False,
    NumberOfIndividual=100,
    NumberOfTransposonTypes=1,
    NumberOfInsertionsPerType=[1],
    FrequencyOfInsertions=[1.0],
    ExcisionRates=[0.01],
    RepairRates=[1],
    InsertionRates=[1],
    HardyWeinberg=False,
    NumberOfGenerations=1000,
    numberOfThreads=8,
)
with open("K1U001.pickle", "wb") as f:
    pickle.dump((sim), f)

# 100
sim   = runBatch(    
    numberOfSimulations=20,
    baseSelection=0,
    baseInsertionProb=1,
    numberOfInsertionSites=10000,
    numberOfChromosomes=6,
    baseRecombinationRate=0.01,
    baseTau=1,
    numberOfPiRNA=6,
    piPercentage=3,
    enablePiRecombination=False,
    NumberOfIndividual=100,
    NumberOfTransposonTypes=1,
    NumberOfInsertionsPerType=[1],
    FrequencyOfInsertions=[1.0],
    ExcisionRates=[0.10],
    RepairRates=[1],
    InsertionRates=[1],
    HardyWeinberg=False,
    NumberOfGenerations=1000,
    numberOfThreads=8,
)
with open("K1U010.pickle", "wb") as f:
    pickle.dump((sim), f)

