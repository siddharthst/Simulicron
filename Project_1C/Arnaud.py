# Init
import os
import sys
import click
from collections import defaultdict


module_path = os.path.abspath(os.path.join("../src/simulicronalpha/"))
if module_path not in sys.path:
    sys.path.append(module_path)

# Imports
import random
import numpy as np
import pandas as pd
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

    copyNumber = result["AvgCopyNum"]
    piOccupancy = result["TEpi"]
    piOccupancy = piOccupancy[1]
    resultToWrite = [copyNumber, piOccupancy]

    with open("%030x" % random.randrange(16 ** 30) + ".txt", "w+") as f:
        for x in zip(*resultToWrite):
            f.write("{0}\t{1}\n".format(*x))

    return 0

commandArgs = defaultdict(list)
for k, v in ((k.lstrip('-'), v) for k,v in (a.split('=') for a in sys.argv[1:])):
    commandArgs[k].append(v)

for key, value in commandArgs.items():
    if key in parameters.keys():
        parameters[key] = type(parameters[key])(value[0])
        value
        print ("Supplied parameter: ", key)
        print ("Supplied value: ", parameters[key])
        print ("Paramter type: ", type(parameters[key]))

run = worker(parameters)
