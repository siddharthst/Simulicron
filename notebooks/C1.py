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

sim   = runBatch(numberOfSimulations=100000, NumberOfTransposonInsertions=1, baseSelection=0.01, 
                    baseTransposition=0, numberOfChromosomes=1, baseRecombinationRate=0.0, 
                    NumberOfIndividual=10, insertionFrequency=0.4, HardyWeinberg=True)
with open("c1P10C0-0.pickle", "wb") as f:
    pickle.dump((sim), f)

sim   = runBatch(numberOfSimulations=100000, NumberOfTransposonInsertions=1, baseSelection=0.1, 
                    baseTransposition=0, numberOfChromosomes=1, baseRecombinationRate=0.0, 
                    NumberOfIndividual=10, insertionFrequency=0.4, HardyWeinberg=True)
with open("c1P10C0-1.pickle", "wb") as f:
    pickle.dump((sim), f)
    
sim   = runBatch(numberOfSimulations=100000, NumberOfTransposonInsertions=1, baseSelection=-0.01, 
                    baseTransposition=0, numberOfChromosomes=1, baseRecombinationRate=0.0, 
                    NumberOfIndividual=10, insertionFrequency=0.4, HardyWeinberg=True)
with open("c1P10C0-5.pickle", "wb") as f:
    pickle.dump((sim), f)

   
# 30
    
sim   = runBatch(numberOfSimulations=100000, NumberOfTransposonInsertions=1, baseSelection=0.01, 
                    baseTransposition=0, numberOfChromosomes=1, baseRecombinationRate=0.0, 
                    NumberOfIndividual=30, insertionFrequency=0.6, HardyWeinberg=True)
with open("c1P30C0-0.pickle", "wb") as f:
    pickle.dump((sim), f)

sim   = runBatch(numberOfSimulations=100000, NumberOfTransposonInsertions=1, baseSelection=0.1, 
                    baseTransposition=0, numberOfChromosomes=1, baseRecombinationRate=0.0, 
                    NumberOfIndividual=30, insertionFrequency=0.6, HardyWeinberg=True)
with open("c1P30C0-1.pickle", "wb") as f:
    pickle.dump((sim), f)
    
sim   = runBatch(numberOfSimulations=100000, NumberOfTransposonInsertions=1, baseSelection=-0.1, 
                    baseTransposition=0, numberOfChromosomes=1, baseRecombinationRate=0.0, 
                    NumberOfIndividual=30, insertionFrequency=0.6, HardyWeinberg=True)
with open("c1P30C0-5.pickle", "wb") as f:
    pickle.dump((sim), f)
    
#50

sim   = runBatch(numberOfSimulations=100000, NumberOfTransposonInsertions=1, baseSelection=0.01, 
                    baseTransposition=0, numberOfChromosomes=1, baseRecombinationRate=0.0, 
                    NumberOfIndividual=50, insertionFrequency=0.4, HardyWeinberg=True)
with open("c1P50C0-0.pickle", "wb") as f:
    pickle.dump((sim), f)

sim   = runBatch(numberOfSimulations=100000, NumberOfTransposonInsertions=1, baseSelection=0.1, 
                    baseTransposition=0, numberOfChromosomes=1, baseRecombinationRate=0.0, 
                    NumberOfIndividual=50, insertionFrequency=0.4, HardyWeinberg=True)
with open("c1P50C0-1.pickle", "wb") as f:
    pickle.dump((sim), f)
    
sim   = runBatch(numberOfSimulations=100000, NumberOfTransposonInsertions=1, baseSelection=-0.01, 
                    baseTransposition=0, numberOfChromosomes=1, baseRecombinationRate=0.0, 
                    NumberOfIndividual=50, insertionFrequency=0.4, HardyWeinberg=True)
with open("c1P50C0-5.pickle", "wb") as f:
    pickle.dump((sim), f)
    
# 70

sim   = runBatch(numberOfSimulations=100000, NumberOfTransposonInsertions=1, baseSelection=0.01, 
                    baseTransposition=0, numberOfChromosomes=1, baseRecombinationRate=0.0, 
                    NumberOfIndividual=70, insertionFrequency=0.3, HardyWeinberg=True)
with open("c1P70C0-0.pickle", "wb") as f:
    pickle.dump((sim), f)

sim   = runBatch(numberOfSimulations=100000, NumberOfTransposonInsertions=1, baseSelection=0.1, 
                    baseTransposition=0, numberOfChromosomes=1, baseRecombinationRate=0.0, 
                    NumberOfIndividual=70, insertionFrequency=0.3, HardyWeinberg=True)
with open("c1P70C0-1.pickle", "wb") as f:
    pickle.dump((sim), f)
    
sim   = runBatch(numberOfSimulations=100000, NumberOfTransposonInsertions=1, baseSelection=-0.01, 
                    baseTransposition=0, numberOfChromosomes=1, baseRecombinationRate=0.0, 
                    NumberOfIndividual=70, insertionFrequency=0.3, HardyWeinberg=True)
with open("c1P70C0-5.pickle", "wb") as f:
    pickle.dump((sim), f)
    
# 100

sim   = runBatch(numberOfSimulations=100000, NumberOfTransposonInsertions=1, baseSelection=0.01, 
                    baseTransposition=0, numberOfChromosomes=1, baseRecombinationRate=0.0, 
                    NumberOfIndividual=100, insertionFrequency=0.5, HardyWeinberg=True)
with open("c1P100C0-0.pickle", "wb") as f:
    pickle.dump((sim), f)

sim   = runBatch(numberOfSimulations=100000, NumberOfTransposonInsertions=1, baseSelection=0.1, 
                    baseTransposition=0, numberOfChromosomes=1, baseRecombinationRate=0.0, 
                    NumberOfIndividual=100, insertionFrequency=0.5, HardyWeinberg=True)
with open("c1P100C0-1.pickle", "wb") as f:
    pickle.dump((sim), f)
    
sim   = runBatch(numberOfSimulations=100000, NumberOfTransposonInsertions=1, baseSelection=-0.01, 
                    baseTransposition=0, numberOfChromosomes=1, baseRecombinationRate=0.0, 
                    NumberOfIndividual=100, insertionFrequency=0.5, HardyWeinberg=True)
with open("c1P100C0-5.pickle", "wb") as f:
    pickle.dump((sim), f)

# 200

sim   = runBatch(numberOfSimulations=100000, NumberOfTransposonInsertions=1, baseSelection=0.01, 
                    baseTransposition=0, numberOfChromosomes=1, baseRecombinationRate=0.0, 
                    NumberOfIndividual=200, insertionFrequency=0.1, HardyWeinberg=True)
with open("c1P200C0-0.pickle", "wb") as f:
    pickle.dump((sim), f)

sim   = runBatch(numberOfSimulations=100000, NumberOfTransposonInsertions=1, baseSelection=0.1, 
                    baseTransposition=0, numberOfChromosomes=1, baseRecombinationRate=0.0, 
                    NumberOfIndividual=200, insertionFrequency=0.1, HardyWeinberg=True)
with open("c1P200C0-1.pickle", "wb") as f:
    pickle.dump((sim), f)
    
sim   = runBatch(numberOfSimulations=100000, NumberOfTransposonInsertions=1, baseSelection=-0.01, 
                    baseTransposition=0, numberOfChromosomes=1, baseRecombinationRate=0.0, 
                    NumberOfIndividual=200, insertionFrequency=0.1, HardyWeinberg=True)
with open("c1P200C0-5.pickle", "wb") as f:
    pickle.dump((sim), f)

# 300

sim   = runBatch(numberOfSimulations=100000, NumberOfTransposonInsertions=1, baseSelection=0.01, 
                    baseTransposition=0, numberOfChromosomes=1, baseRecombinationRate=0.0, 
                    NumberOfIndividual=300, insertionFrequency=0.1, HardyWeinberg=True)
with open("c1P300C0-0.pickle", "wb") as f:
    pickle.dump((sim), f)

sim   = runBatch(numberOfSimulations=100000, NumberOfTransposonInsertions=1, baseSelection=0.1, 
                    baseTransposition=0, numberOfChromosomes=1, baseRecombinationRate=0.0, 
                    NumberOfIndividual=300, insertionFrequency=0.1, HardyWeinberg=True)
with open("c1P300C0-1.pickle", "wb") as f:
    pickle.dump((sim), f)
    
sim   = runBatch(numberOfSimulations=100000, NumberOfTransposonInsertions=1, baseSelection=-0.01, 
                    baseTransposition=0, numberOfChromosomes=1, baseRecombinationRate=0.0, 
                    NumberOfIndividual=300, insertionFrequency=0.1, HardyWeinberg=True)
with open("c1P300C0-5.pickle", "wb") as f:
    pickle.dump((sim), f)


### selection = 0.0

sim   = runBatch(numberOfSimulations=100000, NumberOfTransposonInsertions=1, baseSelection=0.00, 
                    baseTransposition=0, numberOfChromosomes=1, baseRecombinationRate=0.0, 
                    NumberOfIndividual=10, insertionFrequency=0.4, HardyWeinberg=True)
with open("c1P10C00-00.pickle", "wb") as f:
    pickle.dump((sim), f)
    
sim   = runBatch(numberOfSimulations=100000, NumberOfTransposonInsertions=1, baseSelection=0.00, 
                    baseTransposition=0, numberOfChromosomes=1, baseRecombinationRate=0.0, 
                    NumberOfIndividual=30, insertionFrequency=0.6, HardyWeinberg=True)
with open("c1P30C00-00.pickle", "wb") as f:
    pickle.dump((sim), f)
    
sim   = runBatch(numberOfSimulations=100000, NumberOfTransposonInsertions=1, baseSelection=0.00, 
                    baseTransposition=0, numberOfChromosomes=1, baseRecombinationRate=0.0, 
                    NumberOfIndividual=50, insertionFrequency=0.4, HardyWeinberg=True)
with open("c1P50C00-00.pickle", "wb") as f:
    pickle.dump((sim), f)
    
sim   = runBatch(numberOfSimulations=100000, NumberOfTransposonInsertions=1, baseSelection=0.00, 
                    baseTransposition=0, numberOfChromosomes=1, baseRecombinationRate=0.0, 
                    NumberOfIndividual=70, insertionFrequency=0.3, HardyWeinberg=True)
with open("c1P70C00-00.pickle", "wb") as f:
    pickle.dump((sim), f)
    
sim   = runBatch(numberOfSimulations=100000, NumberOfTransposonInsertions=1, baseSelection=0.00, 
                    baseTransposition=0, numberOfChromosomes=1, baseRecombinationRate=0.0, 
                    NumberOfIndividual=100, insertionFrequency=0.5, HardyWeinberg=True)
with open("c1P100C00-00.pickle", "wb") as f:
    pickle.dump((sim), f)
    
sim   = runBatch(numberOfSimulations=100000, NumberOfTransposonInsertions=1, baseSelection=0.00, 
                    baseTransposition=0, numberOfChromosomes=1, baseRecombinationRate=0.0, 
                    NumberOfIndividual=200, insertionFrequency=0.1, HardyWeinberg=True)
with open("c1P200C00-00.pickle", "wb") as f:
    pickle.dump((sim), f)
    
sim   = runBatch(numberOfSimulations=100000, NumberOfTransposonInsertions=1, baseSelection=0.00, 
                    baseTransposition=0, numberOfChromosomes=1, baseRecombinationRate=0.0, 
                    NumberOfIndividual=300, insertionFrequency=0.1, HardyWeinberg=True)
with open("c1P300C00-00.pickle", "wb") as f:
    pickle.dump((sim), f)
