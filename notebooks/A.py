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
    generateFitness,
    recombination,
    transposition,
    runSim,
    runBatch,
)

simP10   = runBatch(numberOfSimulations=100000, NumberOfTransposonInsertions=1, baseSelection=0, baseTransposition=0, numberOfChromosomes=1, baseRecombinationRate=0.0, NumberOfIndividual=10)
with open("a1simP10.pickle", "wb") as f:
    pickle.dump((simP10), f)

simP30   = runBatch(numberOfSimulations=100000, NumberOfTransposonInsertions=1, baseSelection=0, baseTransposition=0, numberOfChromosomes=1, baseRecombinationRate=0.0, NumberOfIndividual=30)
with open("a1simP30.pickle", "wb") as f:
    pickle.dump((simP30), f)
    
simP50   = runBatch(numberOfSimulations=100000, NumberOfTransposonInsertions=1, baseSelection=0, baseTransposition=0, numberOfChromosomes=1, baseRecombinationRate=0.0, NumberOfIndividual=50)
with open("a1simP50.pickle", "wb") as f:
    pickle.dump((simP50), f)
    
simP70   = runBatch(numberOfSimulations=100000, NumberOfTransposonInsertions=1, baseSelection=0, baseTransposition=0, numberOfChromosomes=1, baseRecombinationRate=0.0, NumberOfIndividual=70)
with open("a1simP70.pickle", "wb") as f:
    pickle.dump((simP70), f)
    
simP100  = runBatch(numberOfSimulations=100000, NumberOfTransposonInsertions=1, baseSelection=0, baseTransposition=0, numberOfChromosomes=1, baseRecombinationRate=0.0, NumberOfIndividual=100)
with open("a1simP100.pickle", "wb") as f:
    pickle.dump((simP100), f)
    
simP150  = runBatch(numberOfSimulations=100000, NumberOfTransposonInsertions=1, baseSelection=0, baseTransposition=0, numberOfChromosomes=1, baseRecombinationRate=0.0, NumberOfIndividual=150)
with open("a1simP150.pickle", "wb") as f:
    pickle.dump((simP150), f)

simP200  = runBatch(numberOfSimulations=100000, NumberOfTransposonInsertions=1, baseSelection=0, baseTransposition=0, numberOfChromosomes=1, baseRecombinationRate=0.0, NumberOfIndividual=200)
with open("a1simP200.pickle", "wb") as f:
    pickle.dump((simP200), f)
    
simP300  = runBatch(numberOfSimulations=100000, NumberOfTransposonInsertions=1, baseSelection=0, baseTransposition=0, numberOfChromosomes=1, baseRecombinationRate=0.0, NumberOfIndividual=300)
with open("a1simP300.pickle", "wb") as f:
    pickle.dump((simP300), f)
    
simP400  = runBatch(numberOfSimulations=100000, NumberOfTransposonInsertions=1, baseSelection=0, baseTransposition=0, numberOfChromosomes=1, baseRecombinationRate=0.0, NumberOfIndividual=400)
with open("a1simP400.pickle", "wb") as f:
    pickle.dump((simP400), f)

simP500  = runBatch(numberOfSimulations=100000, NumberOfTransposonInsertions=1, baseSelection=0, baseTransposition=0, numberOfChromosomes=1, baseRecombinationRate=0.0, NumberOfIndividual=500)
with open("a1simP500.pickle", "wb") as f:
    pickle.dump((simP500), f)
    
simP10   = runBatch(numberOfSimulations=100000, NumberOfTransposonInsertions=5, baseSelection=0, baseTransposition=0, numberOfChromosomes=1, baseRecombinationRate=0.0, NumberOfIndividual=10)
with open("a2simP10.pickle", "wb") as f:
    pickle.dump((simP10), f)

simP30   = runBatch(numberOfSimulations=100000, NumberOfTransposonInsertions=15, baseSelection=0, baseTransposition=0, numberOfChromosomes=1, baseRecombinationRate=0.0, NumberOfIndividual=30)
with open("a2simP30.pickle", "wb") as f:
    pickle.dump((simP30), f)
    
simP50   = runBatch(numberOfSimulations=100000, NumberOfTransposonInsertions=25, baseSelection=0, baseTransposition=0, numberOfChromosomes=1, baseRecombinationRate=0.0, NumberOfIndividual=50)
with open("a2simP50.pickle", "wb") as f:
    pickle.dump((simP50), f)
    
simP70   = runBatch(numberOfSimulations=100000, NumberOfTransposonInsertions=35, baseSelection=0, baseTransposition=0, numberOfChromosomes=1, baseRecombinationRate=0.0, NumberOfIndividual=70)
with open("a2simP70.pickle", "wb") as f:
    pickle.dump((simP70), f)
    
simP100  = runBatch(numberOfSimulations=100000, NumberOfTransposonInsertions=50, baseSelection=0, baseTransposition=0, numberOfChromosomes=1, baseRecombinationRate=0.0, NumberOfIndividual=100)
with open("a2simP100.pickle", "wb") as f:
    pickle.dump((simP100), f)
    
simP150  = runBatch(numberOfSimulations=100000, NumberOfTransposonInsertions=75, baseSelection=0, baseTransposition=0, numberOfChromosomes=1, baseRecombinationRate=0.0, NumberOfIndividual=150)
with open("a2simP150.pickle", "wb") as f:
    pickle.dump((simP150), f)

simP200  = runBatch(numberOfSimulations=100000, NumberOfTransposonInsertions=100, baseSelection=0, baseTransposition=0, numberOfChromosomes=1, baseRecombinationRate=0.0, NumberOfIndividual=200)
with open("a2simP200.pickle", "wb") as f:
    pickle.dump((simP200), f)
    
simP300  = runBatch(numberOfSimulations=100000, NumberOfTransposonInsertions=150, baseSelection=0, baseTransposition=0, numberOfChromosomes=1, baseRecombinationRate=0.0, NumberOfIndividual=300)
with open("a2simP300.pickle", "wb") as f:
    pickle.dump((simP300), f)
    
simP400  = runBatch(numberOfSimulations=100000, NumberOfTransposonInsertions=200, baseSelection=0, baseTransposition=0, numberOfChromosomes=1, baseRecombinationRate=0.0, NumberOfIndividual=400)
with open("a2simP400.pickle", "wb") as f:
    pickle.dump((simP400), f)

simP500  = runBatch(numberOfSimulations=100000, NumberOfTransposonInsertions=250, baseSelection=0, baseTransposition=0, numberOfChromosomes=1, baseRecombinationRate=0.0, NumberOfIndividual=500)
with open("a2simP500.pickle", "wb") as f:
    pickle.dump((simP500), f)

simP700  = runBatch(numberOfSimulations=100000, NumberOfTransposonInsertions=1, baseSelection=0, baseTransposition=0, numberOfChromosomes=1, baseRecombinationRate=0.0, NumberOfIndividual=700)
with open("a1simP700.pickle", "wb") as f:
    pickle.dump((simP700), f)

simP700  = runBatch(numberOfSimulations=100000, NumberOfTransposonInsertions=350, baseSelection=0, baseTransposition=0, numberOfChromosomes=1, baseRecombinationRate=0.0, NumberOfIndividual=700)
with open("a2simP700.pickle", "wb") as f:
    pickle.dump((simP700), f)
    
simP1000 = runBatch(numberOfSimulations=100000, NumberOfTransposonInsertions=1, baseSelection=0, baseTransposition=0, numberOfChromosomes=1, baseRecombinationRate=0.0, NumberOfIndividual=1000)
with open("a1simP1000.pickle", "wb") as f:
    pickle.dump((simP1000), f)

simP1000 = runBatch(numberOfSimulations=100000, NumberOfTransposonInsertions=500, baseSelection=0, baseTransposition=0, numberOfChromosomes=1, baseRecombinationRate=0.0, NumberOfIndividual=1000)
with open("a2simP1000.pickle", "wb") as f:
    pickle.dump((simP1000), f)
    
    
