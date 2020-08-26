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
from itertools import repeat

# Simulation imports
from popSim import runSim
from generateSim import generatePopulation, generateGenome, initHGT
from stats import stats
from regulation import regulation
from checkCopyNumber import checkCopyNumber
from fitness import calculateFitness
from transposition import transposition
from recombination import recombination

# Current multiprocessing implementation
from multiprocessing import Process
import concurrent.futures 

# Define the simulation parameters
# max and min
parameters = {
    "generations": 10000,
    "selectionPenaltyMin": 0.0,
    "selectionPenaltyMax": 1.0,
    "etaMin": 0.0,
    "etaMax": 1.0,
    "tauMin": 0.5,
    "tauMax": 1.0,
    "ExcisionRateMainMin": 0.0,
    "ExcisionRateMainMax": 1.0,
    "ExcisionRateHGTMin": 0.0,
    "ExcisionRateHGTMax": 1.0,
    "FrequencyOfInsertionMainMin": 0.1,
    "FrequencyOfInsertionMainMax": 1.0,
    "FrequencyOfInsertionHGTMax": 0.1,
    "FrequencyOfInsertionHGTMin": 1.0,
    "HGTgenerationMin": 5,
    "HGTgenerationMax": 8000,
    "maxProcceses": 30,
}

# Wrapper function for multiprocessing
def worker(parameters):
    # Generate genome and population
    selectionCoef = np.random.uniform(
        parameters["selectionPenaltyMin"],
        parameters["selectionPenaltyMax"],
    )
    FrequencyOfInsertionMain = np.random.uniform(
        parameters["FrequencyOfInsertionMainMin"],
        parameters["FrequencyOfInsertionMainMax"],
    )
    FrequencyOfInsertionHGT = np.random.uniform(
        parameters["ExcisionRateMainMin"],
        parameters["ExcisionRateMainMax"],
    )
    ExcisionRateMain = np.random.uniform(
        parameters["ExcisionRateMainMin"],
        parameters["ExcisionRateMainMax"],
    )
    ExcisionRateHGT = np.random.uniform(
        parameters["ExcisionRateHGTMin"],
        parameters["ExcisionRateHGTMax"],
    )
    HMTgen = np.random.randint(
        parameters["HGTgenerationMin"],
        high=parameters["HGTgenerationMax"],
    )
    tau = np.random.uniform(
        parameters["tauMin"], parameters["tauMax"],
    )
    eta = np.random.uniform(
        parameters["etaMin"], parameters["etaMax"],
    )
    genome, piset, piIndice, rates = generateGenome(
        numberOfInsertionSites=10000,
        numberOfChromosomes=6,
        baseRecombinationRate=0.1,
        baseSelection=selectionCoef,
        baseTau=tau,
    )
    population, transposons, TEset = generatePopulation(
        genome,
        piIndice,
        NumberOfIndividual=1000,
        NumberOfTransposonTypes=2,
        NumberOfInsertionsPerType=[1, 0],
        FrequencyOfInsertions=[
            FrequencyOfInsertionMain,
            FrequencyOfInsertionHGT,
        ],
        ExcisionRates=[ExcisionRateMain, ExcisionRateHGT],
        RepairRates=[1, 1],
        InsertionRates=[1, 1],
    )
    result = runSim(
        genomeMatrix=genome,
        populationMatrix=population,
        transposonMatrix=transposons,
        TEset=TEset,
        NumberOfTransposonInsertions=2,
        generations=parameters["generations"],
        genMap=rates,
        piRNAindices=piIndice,
        simHGT=None,
        HMTgen=HMTgen,
        NumberOfInsertionsPerType=None,
        FrequencyOfInsertions=[
            FrequencyOfInsertionMain,
            FrequencyOfInsertionHGT,
        ],
        ExcisionRates=[ExcisionRateMain, ExcisionRateHGT],
        RepairRates=[1, 1],
        InsertionRates=[1, 1],
        eta=eta,
        tau=tau,
        selPen=selectionCoef,
    )
    with open('%030x' % random.randrange(16**30) + ".pickle", "wb") as f:
        pickle.dump((result), f)

    return 0

countSimulations = 1 
with concurrent.futures.ProcessPoolExecutor(max_workers=parameters["maxProcceses"]) as executor:
    futures = [executor.submit(worker, arg) for arg in repeat(parameters,1000)]
    for future in concurrent.futures.as_completed(futures):
        print (future)

