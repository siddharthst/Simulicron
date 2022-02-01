# Init
import os
import sys

module_path = os.path.abspath(os.path.join("../../src/simulicronalpha/"))
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
    "generations": 4000,
    "individuals": 500,
    "selectionPenaltyMin": -0.001,
    "selectionPenaltyMax": -0.001,
    "etaMin": 0,
    "etaMax": 0,
    "tauMin": 1.0,
    "tauMax": 1.0,
    "ExcisionRateMainMin": 0.05,
    "ExcisionRateMainMax": 0.05,
    "ExcisionRateHGTMin": 0.05,
    "ExcisionRateHGTMax": 0.05,
    "FrequencyOfInsertionMainMin": 0.2,
    "FrequencyOfInsertionMainMax": 0.2,
    "FrequencyOfInsertionHGTMin": 0.2,
    "FrequencyOfInsertionHGTMax": 0.2,
    "HGTgenerationMin": 1000,
    "HGTgenerationMax": 1000,
    "maxProcceses": 60,
    "epistasisCoefficientMax": 0,
    "epistasisCoefficientMin": 0,
    "saveSuffix": "_0",
    "epistasisCoefficientArray":[0]
}

# Wrapper function for multiprocessing
def worker(parameters):
    np.random.seed(int.from_bytes(os.urandom(4), byteorder='little'))
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
        parameters["FrequencyOfInsertionHGTMin"],
        parameters["FrequencyOfInsertionHGTMax"],
    )
    ExcisionRateMain = np.random.uniform(
        parameters["ExcisionRateMainMin"],
        parameters["ExcisionRateMainMax"],
    )
    ExcisionRateHGT = np.random.uniform(
        parameters["ExcisionRateHGTMin"],
        parameters["ExcisionRateHGTMax"],
    )
    HMTgen = random.randint(
        parameters["HGTgenerationMin"],
        parameters["HGTgenerationMax"],
    )
    tau = np.random.uniform(
        parameters["tauMin"], parameters["tauMax"],
    )
    eta = np.random.uniform(
        parameters["etaMin"], parameters["etaMax"],
    )
    epistasisCoefficient = random.choice(parameters["epistasisCoefficientArray"])
    genome, piset, piIndice, rates = generateGenome(
        numberOfInsertionSites=1000,
        numberOfChromosomes=5,
        baseRecombinationRate=0.1,
        baseSelection=selectionCoef,
        baseTau=tau,
    )
    population, transposons, TEset = generatePopulation(
        genome,
        piIndice,
        NumberOfIndividual=parameters["individuals"],
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
        epistasisCoefficient = epistasisCoefficient,
        fitnessFunction = 2,
    )
    with open("./Results_1000/"+'%030x' % random.randrange(16**30) + parameters["saveSuffix"] + ".pickle", "wb") as f:
        pickle.dump((result), f)

    return 0

EpistatsisCoefficients = [-100,-10,-1,-0.01,0,0.01,0.1,1.0]

# Base simulation with just selection 
parameters["selectionPenaltyMin"] = -0.01
parameters["selectionPenaltyMax"] = 0.0
parameters["saveSuffix"] = "_1"
with concurrent.futures.ProcessPoolExecutor(max_workers=parameters["maxProcceses"]) as executor:
    futures = [executor.submit(worker, arg) for arg in repeat(parameters,250)]
    for future in concurrent.futures.as_completed(futures):
        print (future)
parameters["selectionPenaltyMin"] = -0.005
parameters["selectionPenaltyMax"] = -0.005

# Base simulation with just epistasis 
parameters["epistasisCoefficientArray"] = EpistatsisCoefficients
parameters["saveSuffix"] = "_2"
with concurrent.futures.ProcessPoolExecutor(max_workers=parameters["maxProcceses"]) as executor:
    futures = [executor.submit(worker, arg) for arg in repeat(parameters,250)]
    for future in concurrent.futures.as_completed(futures):
        print (future)
parameters["epistasisCoefficientArray"] = [0]


# Base simulation with just coregulation 
parameters["etaMin"] = 0.0
parameters["etaMax"] = 1.0
parameters["saveSuffix"] = "_3"
with concurrent.futures.ProcessPoolExecutor(max_workers=parameters["maxProcceses"]) as executor:
    futures = [executor.submit(worker, arg) for arg in repeat(parameters,250)]
    for future in concurrent.futures.as_completed(futures):
        print (future)
parameters["etaMin"] = 0.0
parameters["etaMax"] = 0.0

# Setting for single TE
parameters["FrequencyOfInsertionHGTMin"] = 0.0,
parameters["FrequencyOfInsertionHGTMax"] = 0.0,
# Single TE simulation with just selection 
parameters["selectionPenaltyMin"] = -0.01
parameters["selectionPenaltyMax"] = 0.0
parameters["saveSuffix"] = "_4"
with concurrent.futures.ProcessPoolExecutor(max_workers=parameters["maxProcceses"]) as executor:
    futures = [executor.submit(worker, arg) for arg in repeat(parameters,250)]
    for future in concurrent.futures.as_completed(futures):
        print (future)
parameters["selectionPenaltyMin"] = -0.005
parameters["selectionPenaltyMax"] = -0.005

# Single TE simulation with just epistasis 
parameters["epistasisCoefficientArray"] = EpistatsisCoefficients
parameters["saveSuffix"] = "_5"
with concurrent.futures.ProcessPoolExecutor(max_workers=parameters["maxProcceses"]) as executor:
    futures = [executor.submit(worker, arg) for arg in repeat(parameters,250)]
    for future in concurrent.futures.as_completed(futures):
        print (future)
parameters["epistasisCoefficientArray"] = [0]


# Single TE simulation with just coregulation 
parameters["etaMin"] = 0.0
parameters["etaMax"] = 1.0
parameters["saveSuffix"] = "_6"
with concurrent.futures.ProcessPoolExecutor(max_workers=parameters["maxProcceses"]) as executor:
    futures = [executor.submit(worker, arg) for arg in repeat(parameters,250)]
    for future in concurrent.futures.as_completed(futures):
        print (future)
parameters["etaMin"] = 0.0
parameters["etaMax"] = 0.0



