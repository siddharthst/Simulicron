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
                    NumberOfIndividual=10, insertionFrequency=0.5, HardyWeinberg=True)
with open("a2P10C0-0.pickle", "wb") as f:
    pickle.dump((sim), f)
    
sim   = runBatch(numberOfSimulations=100000, NumberOfTransposonInsertions=1, baseSelection=0.0, 
                    baseTransposition=0, numberOfChromosomes=1, baseRecombinationRate=0.0, 
                    NumberOfIndividual=30, insertionFrequency=0.5, HardyWeinberg=True)
with open("a2P30C0-0.pickle", "wb") as f:
    pickle.dump((sim), f)

sim   = runBatch(numberOfSimulations=100000, NumberOfTransposonInsertions=1, baseSelection=0.0, 
                    baseTransposition=0, numberOfChromosomes=1, baseRecombinationRate=0.0, 
                    NumberOfIndividual=50, insertionFrequency=0.5, HardyWeinberg=True)
with open("a2P50C0-0.pickle", "wb") as f:
    pickle.dump((sim), f)
    
sim   = runBatch(numberOfSimulations=100000, NumberOfTransposonInsertions=1, baseSelection=0.0, 
                    baseTransposition=0, numberOfChromosomes=1, baseRecombinationRate=0.0, 
                    NumberOfIndividual=70, insertionFrequency=0.5, HardyWeinberg=True)
with open("a2P70C0-0.pickle", "wb") as f:
    pickle.dump((sim), f)

sim   = runBatch(numberOfSimulations=100000, NumberOfTransposonInsertions=1, baseSelection=0.0, 
                    baseTransposition=0, numberOfChromosomes=1, baseRecombinationRate=0.0, 
                    NumberOfIndividual=100, insertionFrequency=0.5, HardyWeinberg=True)
with open("a2P100C0-0.pickle", "wb") as f:
    pickle.dump((sim), f)

sim   = runBatch(numberOfSimulations=100000, NumberOfTransposonInsertions=1, baseSelection=0.0, 
                    baseTransposition=0, numberOfChromosomes=1, baseRecombinationRate=0.0, 
                    NumberOfIndividual=200, insertionFrequency=0.5, HardyWeinberg=True)
with open("a2P200C0-0.pickle", "wb") as f:
    pickle.dump((sim), f)


sim   = runBatch(numberOfSimulations=100000, NumberOfTransposonInsertions=1, baseSelection=0.0, 
                    baseTransposition=0, numberOfChromosomes=1, baseRecombinationRate=0.0, 
                    NumberOfIndividual=300, insertionFrequency=0.5, HardyWeinberg=True)
with open("a2P300C0-0.pickle", "wb") as f:
    pickle.dump((sim), f)
    

sim   = runBatch(numberOfSimulations=100000, NumberOfTransposonInsertions=1, baseSelection=0.0, 
                    baseTransposition=0, numberOfChromosomes=1, baseRecombinationRate=0.0, 
                    NumberOfIndividual=400, insertionFrequency=0.5, HardyWeinberg=True)
with open("a2P400C0-0.pickle", "wb") as f:
    pickle.dump((sim), f)
    

sim   = runBatch(numberOfSimulations=100000, NumberOfTransposonInsertions=1, baseSelection=0.0, 
                    baseTransposition=0, numberOfChromosomes=1, baseRecombinationRate=0.0, 
                    NumberOfIndividual=500, insertionFrequency=0.5, HardyWeinberg=True)
with open("a2P500C0-0.pickle", "wb") as f:
    pickle.dump((sim), f)


sim   = runBatch(numberOfSimulations=100000, NumberOfTransposonInsertions=1, baseSelection=0.0, 
                    baseTransposition=0, numberOfChromosomes=1, baseRecombinationRate=0.0, 
                    NumberOfIndividual=600, insertionFrequency=0.5, HardyWeinberg=True)
with open("a2P600C0-0.pickle", "wb") as f:
    pickle.dump((sim), f)
    

sim   = runBatch(numberOfSimulations=100000, NumberOfTransposonInsertions=1, baseSelection=0.0, 
                    baseTransposition=0, numberOfChromosomes=1, baseRecombinationRate=0.0, 
                    NumberOfIndividual=700, insertionFrequency=0.5, HardyWeinberg=True)
with open("a2P700C0-0.pickle", "wb") as f:
    pickle.dump((sim), f)
    

sim   = runBatch(numberOfSimulations=100000, NumberOfTransposonInsertions=1, baseSelection=0.0, 
                    baseTransposition=0, numberOfChromosomes=1, baseRecombinationRate=0.0, 
                    NumberOfIndividual=1000, insertionFrequency=0.5, HardyWeinberg=True)
with open("a2P1000C0-0.pickle", "wb") as f:
    pickle.dump((sim), f)
