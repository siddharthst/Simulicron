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

sim   = runBatch(numberOfSimulations=100000, NumberOfTransposonInsertions=1, baseSelection=0.0, 
                    baseTransposition=0, numberOfChromosomes=1, baseRecombinationRate=0.0, 
                    NumberOfIndividual=12, insertionFrequency=0.5, HardyWeinberg=True)
with open("a2P12C0-0.pickle", "wb") as f:
    pickle.dump((sim), f)
    
sim   = runBatch(numberOfSimulations=100000, NumberOfTransposonInsertions=1, baseSelection=0.0, 
                    baseTransposition=0, numberOfChromosomes=1, baseRecombinationRate=0.0, 
                    NumberOfIndividual=32, insertionFrequency=0.5, HardyWeinberg=True)
with open("a2P32C0-0.pickle", "wb") as f:
    pickle.dump((sim), f)

sim   = runBatch(numberOfSimulations=100000, NumberOfTransposonInsertions=1, baseSelection=0.0, 
                    baseTransposition=0, numberOfChromosomes=1, baseRecombinationRate=0.0, 
                    NumberOfIndividual=64, insertionFrequency=0.5, HardyWeinberg=True)
with open("a2P64C0-0.pickle", "wb") as f:
    pickle.dump((sim), f)
    
sim   = runBatch(numberOfSimulations=100000, NumberOfTransposonInsertions=1, baseSelection=0.0, 
                    baseTransposition=0, numberOfChromosomes=1, baseRecombinationRate=0.0, 
                    NumberOfIndividual=128, insertionFrequency=0.5, HardyWeinberg=True)
with open("a2P128C0-0.pickle", "wb") as f:
    pickle.dump((sim), f)

sim   = runBatch(numberOfSimulations=100000, NumberOfTransposonInsertions=1, baseSelection=0.0, 
                    baseTransposition=0, numberOfChromosomes=1, baseRecombinationRate=0.0, 
                    NumberOfIndividual=192, insertionFrequency=0.5, HardyWeinberg=True)
with open("a2P192C0-0.pickle", "wb") as f:
    pickle.dump((sim), f)
