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
sim   = runBatch(numberOfSimulations=1000, NumberOfTransposonInsertions=1, baseSelection=0, 
                    baseExcision=0.02, numberOfChromosomes=1, baseRecombinationRate=0.0, 
                    NumberOfIndividual=100, InsertIntoAll=True)
with open("d1P100u02.pickle", "wb") as f:
    pickle.dump((sim), f)

# 200
sim   = runBatch(numberOfSimulations=1000, NumberOfTransposonInsertions=1, baseSelection=0, 
                    baseExcision=0.02, numberOfChromosomes=1, baseRecombinationRate=0.0, 
                    NumberOfIndividual=200, InsertIntoAll=True)
with open("d1P200u02.pickle", "wb") as f:
    pickle.dump((sim), f)
    

        
