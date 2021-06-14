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
    "generations": 5000,
    "individuals": 500,
    "selectionPenalty": 0.001,
    "tau": 1.0,
    "ExcisionRate": 1.0,
    "FrequencyOfInsertion": 1.0,
    "Chromosomes": 6,
    "RecombinationRate": 0.1,
    "NumberOfInsertions": 1,
    "piRNASelection":False,
}

# Wrapper function for multiprocessing
def worker(parameters):
    genome, piset, piIndice, rates = generateGenome(
        numberOfInsertionSites=10000,
        numberOfChromosomes=parameters["Chromosomes"],
        baseRecombinationRate=parameters["RecombinationRate"],
        baseSelection=parameters["selectionPenalty"],
        baseTau=parameters["tau"],
        DisablePiSelection=parameters["piRNASelection"]
    )
    population, transposons, TEset = generatePopulation(
        genome,
        piIndice,
        NumberOfIndividual=parameters["individuals"],
        NumberOfTransposonTypes=1,
        NumberOfInsertionsPerType=[parameters["NumberOfInsertions"]],
        FrequencyOfInsertions=[
            parameters["FrequencyOfInsertion"],
        ],
        ExcisionRates=[parameters["ExcisionRate"]],
        RepairRates=[1],
        InsertionRates=[1],
    )
    result = runSim(
        genomeMatrix=genome,
        populationMatrix=population,
        transposonMatrix=transposons,
        TEset=TEset,
        NumberOfTransposonInsertions=1,
        generations=parameters["generations"],
        genMap=rates,
        piRNAindices=piIndice,
        simHGT=None,
        HMTgen=None,
        NumberOfInsertionsPerType=None,
        FrequencyOfInsertions=None,
        ExcisionRates=None,
        RepairRates=None,
        InsertionRates=None,
        eta=0,
        tau=parameters["tau"],
        selPen=parameters["selectionPenalty"],
    )
    with open("etaVsHG/" + "%030x" % random.randrange(16 ** 30) + ".pickle", "wb") as f:
        pickle.dump((result), f)

    return 0


countSimulations = 1
