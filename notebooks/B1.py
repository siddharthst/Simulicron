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

sim   = runBatch(numberOfSimulations=100000, NumberOfTransposonInsertions=2, baseSelection=0, 
                    baseTransposition=0, numberOfChromosomes=1, baseRecombinationRate=0.0, 
                    consecutiveTransposons=True, changeRecombination=True, baseTrRecombination=0.0,
                    NumberOfIndividual=10)
with open("b1P10C0-0.pickle", "wb") as f:
    pickle.dump((sim), f)

sim   = runBatch(numberOfSimulations=100000, NumberOfTransposonInsertions=2, baseSelection=0, 
                    baseTransposition=0, numberOfChromosomes=1, baseRecombinationRate=0.0, 
                    consecutiveTransposons=True, changeRecombination=True, baseTrRecombination=0.1,
                    NumberOfIndividual=10)
with open("b1P10C0-1.pickle", "wb") as f:
    pickle.dump((sim), f)
    
sim   = runBatch(numberOfSimulations=100000, NumberOfTransposonInsertions=2, baseSelection=0, 
                    baseTransposition=0, numberOfChromosomes=1, baseRecombinationRate=0.0, 
                    consecutiveTransposons=True, changeRecombination=True, baseTrRecombination=0.5,
                    NumberOfIndividual=10)
with open("b1P10C0-5.pickle", "wb") as f:
    pickle.dump((sim), f)

   
# 30
    
sim   = runBatch(numberOfSimulations=100000, NumberOfTransposonInsertions=2, baseSelection=0, 
                    baseTransposition=0, numberOfChromosomes=1, baseRecombinationRate=0.0, 
                    consecutiveTransposons=True, changeRecombination=True, baseTrRecombination=0.0,
                    NumberOfIndividual=30)
with open("b1P30C0-0.pickle", "wb") as f:
    pickle.dump((sim), f)

sim   = runBatch(numberOfSimulations=100000, NumberOfTransposonInsertions=2, baseSelection=0, 
                    baseTransposition=0, numberOfChromosomes=1, baseRecombinationRate=0.0, 
                    consecutiveTransposons=True, changeRecombination=True, baseTrRecombination=0.1,
                    NumberOfIndividual=30)
with open("b1P30C0-1.pickle", "wb") as f:
    pickle.dump((sim), f)
    
sim   = runBatch(numberOfSimulations=100000, NumberOfTransposonInsertions=2, baseSelection=0, 
                    baseTransposition=0, numberOfChromosomes=1, baseRecombinationRate=0.0, 
                    consecutiveTransposons=True, changeRecombination=True, baseTrRecombination=0.5,
                    NumberOfIndividual=30)
with open("b1P30C0-5.pickle", "wb") as f:
    pickle.dump((sim), f)
    
#50

sim   = runBatch(numberOfSimulations=100000, NumberOfTransposonInsertions=2, baseSelection=0, 
                    baseTransposition=0, numberOfChromosomes=1, baseRecombinationRate=0.0, 
                    consecutiveTransposons=True, changeRecombination=True, baseTrRecombination=0.0,
                    NumberOfIndividual=50)
with open("b1P50C0-0.pickle", "wb") as f:
    pickle.dump((sim), f)

sim   = runBatch(numberOfSimulations=100000, NumberOfTransposonInsertions=2, baseSelection=0, 
                    baseTransposition=0, numberOfChromosomes=1, baseRecombinationRate=0.0, 
                    consecutiveTransposons=True, changeRecombination=True, baseTrRecombination=0.1,
                    NumberOfIndividual=50)
with open("b1P50C0-1.pickle", "wb") as f:
    pickle.dump((sim), f)
    
sim   = runBatch(numberOfSimulations=100000, NumberOfTransposonInsertions=2, baseSelection=0, 
                    baseTransposition=0, numberOfChromosomes=1, baseRecombinationRate=0.0, 
                    consecutiveTransposons=True, changeRecombination=True, baseTrRecombination=0.5,
                    NumberOfIndividual=50)
with open("b1P50C0-5.pickle", "wb") as f:
    pickle.dump((sim), f)
    
# 70

sim   = runBatch(numberOfSimulations=100000, NumberOfTransposonInsertions=2, baseSelection=0, 
                    baseTransposition=0, numberOfChromosomes=1, baseRecombinationRate=0.0, 
                    consecutiveTransposons=True, changeRecombination=True, baseTrRecombination=0.0,
                    NumberOfIndividual=70)
with open("b1P70C0-0.pickle", "wb") as f:
    pickle.dump((sim), f)

sim   = runBatch(numberOfSimulations=100000, NumberOfTransposonInsertions=2, baseSelection=0, 
                    baseTransposition=0, numberOfChromosomes=1, baseRecombinationRate=0.0, 
                    consecutiveTransposons=True, changeRecombination=True, baseTrRecombination=0.1,
                    NumberOfIndividual=70)
with open("b1P70C0-1.pickle", "wb") as f:
    pickle.dump((sim), f)
    
sim   = runBatch(numberOfSimulations=100000, NumberOfTransposonInsertions=2, baseSelection=0, 
                    baseTransposition=0, numberOfChromosomes=1, baseRecombinationRate=0.0, 
                    consecutiveTransposons=True, changeRecombination=True, baseTrRecombination=0.5,
                    NumberOfIndividual=70)
with open("b1P70C0-5.pickle", "wb") as f:
    pickle.dump((sim), f)
    
# 100

sim   = runBatch(numberOfSimulations=100000, NumberOfTransposonInsertions=2, baseSelection=0, 
                    baseTransposition=0, numberOfChromosomes=1, baseRecombinationRate=0.0, 
                    consecutiveTransposons=True, changeRecombination=True, baseTrRecombination=0.0,
                    NumberOfIndividual=100)
with open("b1P100C0-0.pickle", "wb") as f:
    pickle.dump((sim), f)

sim   = runBatch(numberOfSimulations=100000, NumberOfTransposonInsertions=2, baseSelection=0, 
                    baseTransposition=0, numberOfChromosomes=1, baseRecombinationRate=0.0, 
                    consecutiveTransposons=True, changeRecombination=True, baseTrRecombination=0.1,
                    NumberOfIndividual=100)
with open("b1P100C0-1.pickle", "wb") as f:
    pickle.dump((sim), f)
    
sim   = runBatch(numberOfSimulations=100000, NumberOfTransposonInsertions=2, baseSelection=0, 
                    baseTransposition=0, numberOfChromosomes=1, baseRecombinationRate=0.0, 
                    consecutiveTransposons=True, changeRecombination=True, baseTrRecombination=0.5,
                    NumberOfIndividual=100)
with open("b1P100C0-5.pickle", "wb") as f:
    pickle.dump((sim), f)

# 200

sim   = runBatch(numberOfSimulations=100000, NumberOfTransposonInsertions=2, baseSelection=0, 
                    baseTransposition=0, numberOfChromosomes=1, baseRecombinationRate=0.0, 
                    consecutiveTransposons=True, changeRecombination=True, baseTrRecombination=0.0,
                    NumberOfIndividual=200)
with open("b1P200C0-0.pickle", "wb") as f:
    pickle.dump((sim), f)

sim   = runBatch(numberOfSimulations=100000, NumberOfTransposonInsertions=2, baseSelection=0, 
                    baseTransposition=0, numberOfChromosomes=1, baseRecombinationRate=0.0, 
                    consecutiveTransposons=True, changeRecombination=True, baseTrRecombination=0.1,
                    NumberOfIndividual=200)
with open("b1P200C0-1.pickle", "wb") as f:
    pickle.dump((sim), f)
    
sim   = runBatch(numberOfSimulations=100000, NumberOfTransposonInsertions=2, baseSelection=0, 
                    baseTransposition=0, numberOfChromosomes=1, baseRecombinationRate=0.0, 
                    consecutiveTransposons=True, changeRecombination=True, baseTrRecombination=0.5,
                    NumberOfIndividual=200)
with open("b1P200C0-5.pickle", "wb") as f:
    pickle.dump((sim), f)

# 300

sim   = runBatch(numberOfSimulations=100000, NumberOfTransposonInsertions=2, baseSelection=0, 
                    baseTransposition=0, numberOfChromosomes=1, baseRecombinationRate=0.0, 
                    consecutiveTransposons=True, changeRecombination=True, baseTrRecombination=0.0,
                    NumberOfIndividual=300)
with open("b1P300C0-0.pickle", "wb") as f:
    pickle.dump((sim), f)

sim   = runBatch(numberOfSimulations=100000, NumberOfTransposonInsertions=2, baseSelection=0, 
                    baseTransposition=0, numberOfChromosomes=1, baseRecombinationRate=0.0, 
                    consecutiveTransposons=True, changeRecombination=True, baseTrRecombination=0.1,
                    NumberOfIndividual=300)
with open("b1P300C0-1.pickle", "wb") as f:
    pickle.dump((sim), f)
    
sim   = runBatch(numberOfSimulations=100000, NumberOfTransposonInsertions=2, baseSelection=0, 
                    baseTransposition=0, numberOfChromosomes=1, baseRecombinationRate=0.0, 
                    consecutiveTransposons=True, changeRecombination=True, baseTrRecombination=0.5,
                    NumberOfIndividual=300)
with open("b1P300C0-5.pickle", "wb") as f:
    pickle.dump((sim), f)

# 400

sim   = runBatch(numberOfSimulations=100000, NumberOfTransposonInsertions=2, baseSelection=0, 
                    baseTransposition=0, numberOfChromosomes=1, baseRecombinationRate=0.0, 
                    consecutiveTransposons=True, changeRecombination=True, baseTrRecombination=0.0,
                    NumberOfIndividual=400)
with open("b1P400C0-0.pickle", "wb") as f:
    pickle.dump((sim), f)

sim   = runBatch(numberOfSimulations=100000, NumberOfTransposonInsertions=2, baseSelection=0, 
                    baseTransposition=0, numberOfChromosomes=1, baseRecombinationRate=0.0, 
                    consecutiveTransposons=True, changeRecombination=True, baseTrRecombination=0.1,
                    NumberOfIndividual=400)
with open("b1P400C0-1.pickle", "wb") as f:
    pickle.dump((sim), f)
    
sim   = runBatch(numberOfSimulations=100000, NumberOfTransposonInsertions=2, baseSelection=0, 
                    baseTransposition=0, numberOfChromosomes=1, baseRecombinationRate=0.0, 
                    consecutiveTransposons=True, changeRecombination=True, baseTrRecombination=0.5,
                    NumberOfIndividual=400)
with open("b1P400C0-5.pickle", "wb") as f:
    pickle.dump((sim), f)
    
# 500

sim   = runBatch(numberOfSimulations=100000, NumberOfTransposonInsertions=2, baseSelection=0, 
                    baseTransposition=0, numberOfChromosomes=1, baseRecombinationRate=0.0, 
                    consecutiveTransposons=True, changeRecombination=True, baseTrRecombination=0.0,
                    NumberOfIndividual=500)
with open("b1P500C0-0.pickle", "wb") as f:
    pickle.dump((sim), f)

sim   = runBatch(numberOfSimulations=100000, NumberOfTransposonInsertions=2, baseSelection=0, 
                    baseTransposition=0, numberOfChromosomes=1, baseRecombinationRate=0.0, 
                    consecutiveTransposons=True, changeRecombination=True, baseTrRecombination=0.1,
                    NumberOfIndividual=500)
with open("b1P500C0-1.pickle", "wb") as f:
    pickle.dump((sim), f)
    
sim   = runBatch(numberOfSimulations=100000, NumberOfTransposonInsertions=2, baseSelection=0, 
                    baseTransposition=0, numberOfChromosomes=1, baseRecombinationRate=0.0, 
                    consecutiveTransposons=True, changeRecombination=True, baseTrRecombination=0.5,
                    NumberOfIndividual=500)
with open("b1P500C0-5.pickle", "wb") as f:
    pickle.dump((sim), f)
