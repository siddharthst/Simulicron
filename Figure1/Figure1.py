# Init
import os
import sys
import json

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
    epistasisCoefficient = np.random.uniform(
        parameters["epistasisCoefficientMax"], parameters["epistasisCoefficientMin"],
    )
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
    with open("./Results/"+'%030x' % random.randrange(16**30) + parameters["saveSuffix"] + ".pickle", "wb") as f:
        pickle.dump((result), f)

    return 0

with open('"../../Default.parameters', 'r') as file:
    parameters = json.load(file)

# General simulations
with concurrent.futures.ProcessPoolExecutor(max_workers=parameters["maxProcceses"]) as executor:
    futures = [executor.submit(worker, arg) for arg in repeat(parameters,2000)]
    for future in concurrent.futures.as_completed(futures):
        print (future)

# Set the four corners
epistasis = [parameters["etaMin"],parameters["etaMax"]]
generationOfInfection = [0,1000]
for i in epistasis:
    for k in generationOfInfection:
        parameters["etaMin"] = i
        parameters["etaMax"] = i
        parameters["HGTgenerationMin"] = k
        parameters["HGTgenerationMax"] = k
        with concurrent.futures.ProcessPoolExecutor(max_workers=parameters["maxProcceses"]) as executor:
            futures = [executor.submit(worker, arg) for arg in repeat(parameters,250)]
            for future in concurrent.futures.as_completed(futures):
                print (future)
