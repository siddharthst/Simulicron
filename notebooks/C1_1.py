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

from popSim import (
    generateTransposon,
    generatePopulation,
    calculateFitness,
    generateGenome,
    recombination,
    createData,
    transposition,
    runSim,
    runBatch,
)

# 30
sim   = runBatch(numberOfSimulations=100000, NumberOfTransposonInsertions=1, baseSelection=0.01, 
                    baseTransposition=0, numberOfChromosomes=1, baseRecombinationRate=0.0, 
                    NumberOfIndividual=30, insertionFrequency=0.2, HardyWeinberg=True)
with open("c1P30C0-01.pickle", "wb") as f:
    pickle.dump((sim), f)

sim   = runBatch(numberOfSimulations=100000, NumberOfTransposonInsertions=1, baseSelection=0.05, 
                    baseTransposition=0, numberOfChromosomes=1, baseRecombinationRate=0.0, 
                    NumberOfIndividual=30, insertionFrequency=0.2, HardyWeinberg=True)
with open("c1P30C0-05.pickle", "wb") as f:
    pickle.dump((sim), f)
    
sim   = runBatch(numberOfSimulations=100000, NumberOfTransposonInsertions=1, baseSelection=0.10, 
                    baseTransposition=0, numberOfChromosomes=1, baseRecombinationRate=0.0, 
                    NumberOfIndividual=30, insertionFrequency=0.2, HardyWeinberg=True)
with open("c1P30C0-10.pickle", "wb") as f:
    pickle.dump((sim), f)
    
sim   = runBatch(numberOfSimulations=100000, NumberOfTransposonInsertions=1, baseSelection=0.00, 
                    baseTransposition=0, numberOfChromosomes=1, baseRecombinationRate=0.0, 
                    NumberOfIndividual=30, insertionFrequency=0.2, HardyWeinberg=True)
with open("c1P30C0-00.pickle", "wb") as f:
    pickle.dump((sim), f)

    
sim   = runBatch(numberOfSimulations=100000, NumberOfTransposonInsertions=1, baseSelection=-0.01, 
                    baseTransposition=0, numberOfChromosomes=1, baseRecombinationRate=0.0, 
                    NumberOfIndividual=30, insertionFrequency=0.2, HardyWeinberg=True)
with open("c1P30C1-01.pickle", "wb") as f:
    pickle.dump((sim), f)
