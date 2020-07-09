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
sim   = runBatch(numberOfSimulations=1, NumberOfTransposonInsertions=5, baseSelection=-0.02, 
                    baseExcision=0.02, numberOfChromosomes=5, baseRecombinationRate=0.1, 
                    NumberOfIndividual=1000, InsertIntoAll=True, NumberOfGenerations=200)
with open("d2P1000u02.pickle", "wb") as f:
    pickle.dump((sim), f)

# 100
sim   = runBatch(numberOfSimulations=1, NumberOfTransposonInsertions=5, baseSelection=-0.02, 
                    baseExcision=0.03, numberOfChromosomes=5, baseRecombinationRate=0.1, 
                    NumberOfIndividual=1000, InsertIntoAll=True, NumberOfGenerations=200)
with open("d2P1000u03.pickle", "wb") as f:
    pickle.dump((sim), f)

# 100
sim   = runBatch(numberOfSimulations=1, NumberOfTransposonInsertions=5, baseSelection=-0.02, 
                    baseExcision=0.04, numberOfChromosomes=5, baseRecombinationRate=0.1, 
                    NumberOfIndividual=1000, InsertIntoAll=True, NumberOfGenerations=200)
with open("d2P1000u04.pickle", "wb") as f:
    pickle.dump((sim), f)


sim   = runBatch(numberOfSimulations=1, NumberOfTransposonInsertions=5, baseSelection=-0.02, 
                    baseExcision=0.05, numberOfChromosomes=5, baseRecombinationRate=0.1, 
                    NumberOfIndividual=1000, InsertIntoAll=True, NumberOfGenerations=200)
with open("d2P1000u05.pickle", "wb") as f:
    pickle.dump((sim), f)

sim   = runBatch(numberOfSimulations=1, NumberOfTransposonInsertions=5, baseSelection=-0.02, 
                    baseExcision=0.06, numberOfChromosomes=5, baseRecombinationRate=0.1, 
                    NumberOfIndividual=1000, InsertIntoAll=True, NumberOfGenerations=200)
with open("d2P1000u06.pickle", "wb") as f:
    pickle.dump((sim), f)
