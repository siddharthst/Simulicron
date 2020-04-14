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


# 100
sim   = runBatch(numberOfSimulations=1, NumberOfTransposonInsertions=3, baseSelection=0, 
                    baseExcision=0.025, numberOfChromosomes=1, baseRecombinationRate=0.0, 
                    NumberOfIndividual=1000, InsertIntoAll=True, NumberOfGenerations=100)
with open("d1P1000u025.pickle", "wb") as f:
    pickle.dump((sim), f)

# 100
sim   = runBatch(numberOfSimulations=1, NumberOfTransposonInsertions=3, baseSelection=0, 
                    baseExcision=0.035, numberOfChromosomes=1, baseRecombinationRate=0.0, 
                    NumberOfIndividual=1000, InsertIntoAll=True, NumberOfGenerations=100)
with open("d1P1000u035.pickle", "wb") as f:
    pickle.dump((sim), f)

# 100
sim   = runBatch(numberOfSimulations=1, NumberOfTransposonInsertions=3, baseSelection=0, 
                    baseExcision=0.045, numberOfChromosomes=1, baseRecombinationRate=0.0, 
                    NumberOfIndividual=1000, InsertIntoAll=True, NumberOfGenerations=100)
with open("d1P1000u045.pickle", "wb") as f:
    pickle.dump((sim), f)

# 100
sim   = runBatch(numberOfSimulations=1, NumberOfTransposonInsertions=3, baseSelection=0, 
                    baseExcision=0.055, numberOfChromosomes=1, baseRecombinationRate=0.0, 
                    NumberOfIndividual=1000, InsertIntoAll=True, NumberOfGenerations=100)
with open("d1P1000u055.pickle", "wb") as f:
    pickle.dump((sim), f)




               
