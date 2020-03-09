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

# 10
sim   = runBatch(numberOfSimulations=100000, NumberOfTransposonInsertions=1, baseSelection=0.01, 
                    baseTransposition=0, numberOfChromosomes=1, baseRecombinationRate=0.0, 
                    NumberOfIndividual=10, insertionFrequency=0.2, HardyWeinberg=True)
with open("c1P10C0-01.pickle", "wb") as f:
    pickle.dump((sim), f)

sim   = runBatch(numberOfSimulations=100000, NumberOfTransposonInsertions=1, baseSelection=0.1, 
                    baseTransposition=0, numberOfChromosomes=1, baseRecombinationRate=0.0, 
                    NumberOfIndividual=10, insertionFrequency=0.2, HardyWeinberg=True)
with open("c1P10C0-10.pickle", "wb") as f:
    pickle.dump((sim), f)
    
sim   = runBatch(numberOfSimulations=100000, NumberOfTransposonInsertions=1, baseSelection=0.1, 
                    baseTransposition=0, numberOfChromosomes=1, baseRecombinationRate=0.0, 
                    NumberOfIndividual=10, insertionFrequency=0.2, HardyWeinberg=True)
with open("c1P10C0-00.pickle", "wb") as f:
    pickle.dump((sim), f)

# 30
sim   = runBatch(numberOfSimulations=100000, NumberOfTransposonInsertions=1, baseSelection=0.01, 
                    baseTransposition=0, numberOfChromosomes=1, baseRecombinationRate=0.0, 
                    NumberOfIndividual=30, insertionFrequency=0.2, HardyWeinberg=True)
with open("c1P30C0-01.pickle", "wb") as f:
    pickle.dump((sim), f)

sim   = runBatch(numberOfSimulations=100000, NumberOfTransposonInsertions=1, baseSelection=0.1, 
                    baseTransposition=0, numberOfChromosomes=1, baseRecombinationRate=0.0, 
                    NumberOfIndividual=30, insertionFrequency=0.2, HardyWeinberg=True)
with open("c1P30C0-10.pickle", "wb") as f:
    pickle.dump((sim), f)
    
sim   = runBatch(numberOfSimulations=100000, NumberOfTransposonInsertions=1, baseSelection=0.1, 
                    baseTransposition=0, numberOfChromosomes=1, baseRecombinationRate=0.0, 
                    NumberOfIndividual=30, insertionFrequency=0.2, HardyWeinberg=True)
with open("c1P30C0-00.pickle", "wb") as f:
    pickle.dump((sim), f)

# 50
sim   = runBatch(numberOfSimulations=100000, NumberOfTransposonInsertions=1, baseSelection=0.01, 
                    baseTransposition=0, numberOfChromosomes=1, baseRecombinationRate=0.0, 
                    NumberOfIndividual=50, insertionFrequency=0.2, HardyWeinberg=True)
with open("c1P50C0-01.pickle", "wb") as f:
    pickle.dump((sim), f)

sim   = runBatch(numberOfSimulations=100000, NumberOfTransposonInsertions=1, baseSelection=0.1, 
                    baseTransposition=0, numberOfChromosomes=1, baseRecombinationRate=0.0, 
                    NumberOfIndividual=50, insertionFrequency=0.2, HardyWeinberg=True)
with open("c1P50C0-10.pickle", "wb") as f:
    pickle.dump((sim), f)
    
sim   = runBatch(numberOfSimulations=100000, NumberOfTransposonInsertions=1, baseSelection=0.1, 
                    baseTransposition=0, numberOfChromosomes=1, baseRecombinationRate=0.0, 
                    NumberOfIndividual=50, insertionFrequency=0.2, HardyWeinberg=True)
with open("c1P50C0-00.pickle", "wb") as f:
    pickle.dump((sim), f)
    
# 70
sim   = runBatch(numberOfSimulations=100000, NumberOfTransposonInsertions=1, baseSelection=0.01, 
                    baseTransposition=0, numberOfChromosomes=1, baseRecombinationRate=0.0, 
                    NumberOfIndividual=70, insertionFrequency=0.2, HardyWeinberg=True)
with open("c1P70C0-01.pickle", "wb") as f:
    pickle.dump((sim), f)

sim   = runBatch(numberOfSimulations=100000, NumberOfTransposonInsertions=1, baseSelection=0.1, 
                    baseTransposition=0, numberOfChromosomes=1, baseRecombinationRate=0.0, 
                    NumberOfIndividual=70, insertionFrequency=0.2, HardyWeinberg=True)
with open("c1P70C0-10.pickle", "wb") as f:
    pickle.dump((sim), f)
    
sim   = runBatch(numberOfSimulations=100000, NumberOfTransposonInsertions=1, baseSelection=0.1, 
                    baseTransposition=0, numberOfChromosomes=1, baseRecombinationRate=0.0, 
                    NumberOfIndividual=70, insertionFrequency=0.2, HardyWeinberg=True)
with open("c1P70C0-00.pickle", "wb") as f:
    pickle.dump((sim), f)

# 100
sim   = runBatch(numberOfSimulations=100000, NumberOfTransposonInsertions=1, baseSelection=0.01, 
                    baseTransposition=0, numberOfChromosomes=1, baseRecombinationRate=0.0, 
                    NumberOfIndividual=100, insertionFrequency=0.2, HardyWeinberg=True)
with open("c1P100C0-01.pickle", "wb") as f:
    pickle.dump((sim), f)

sim   = runBatch(numberOfSimulations=100000, NumberOfTransposonInsertions=1, baseSelection=0.1, 
                    baseTransposition=0, numberOfChromosomes=1, baseRecombinationRate=0.0, 
                    NumberOfIndividual=100, insertionFrequency=0.2, HardyWeinberg=True)
with open("c1P100C0-10.pickle", "wb") as f:
    pickle.dump((sim), f)
    
sim   = runBatch(numberOfSimulations=100000, NumberOfTransposonInsertions=1, baseSelection=0.1, 
                    baseTransposition=0, numberOfChromosomes=1, baseRecombinationRate=0.0, 
                    NumberOfIndividual=100, insertionFrequency=0.2, HardyWeinberg=True)
with open("c1P100C0-00.pickle", "wb") as f:
    pickle.dump((sim), f)

