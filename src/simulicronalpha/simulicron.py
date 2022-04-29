# README
# Arguments can be passed with ArgName=ArgValue
# Arguments (default):
# generations = Number of generations (5000)
# individuals = Population size (500)
# selectionPenalty  = Selection penalty (0.001)
# tau = piRNA regulation strength (1.0)
# ExcisionRate = Transposition rate (1.0)
# FrequencyOfInsertion = Frequency of individuals with insertion(1.0)
# Chromosomes = Number of chromosomes (6)
# RecombinationRate = Recombination rate (0.1)
# NumberOfInsertions = Number of insertions per individual (1)
# piRNASelection = piRNA selection penalty (False) - Anything other than
# default will be a float defining the value - can be 0
# piPercentage = Percentage of genome which will make up for piRNA (3)
# numberOfPiRNA = Number of piRNA clusters (6)
# disablePiRecombination = Disable recombination in piRNA cluster (False)
# FileName = Output file name (DefaultOut.txt)




# Init
import os
import sys
from collections import defaultdict
from distutils.util import strtobool

import pathlib


module_path = os.path.abspath(os.path.join(pathlib.Path(__file__).parent.absolute(), "../src/simulicronalpha/"))
if module_path not in sys.path:
    sys.path.append(module_path)

# Imports
import random

# Simulation imports
from popSim import runSim
from generateSim import generatePopulation, generateGenome

# Define the simulation parameters
# max and min
parameters = {
    "generations": 5000,
    "individuals": 500,
    "loci": 1000,
    "selectionPenalty": 0.001,
    "tau": 1.0,
    "ExcisionRate": 1.0,
    "FrequencyOfInsertion": 1.0,
    "Chromosomes": 6,
    "RecombinationRate": 0.1,
    "NumberOfInsertions": 1,
    "piRNASelection": False,
    "piPercentage": 3,
    "numberOfPiRNA": 6,
    "disablePiRecombination": True,
    "regulationStr": 0.0,
    "FileName": "DefaultOut.txt",
}

# Wrapper function for multiprocessing
def worker(parameters):
    genome, piset, piIndice, rates = generateGenome(
        numberOfInsertionSites=parameters["loci"],
        numberOfChromosomes=parameters["Chromosomes"],
        baseRecombinationRate=parameters["RecombinationRate"],
        baseSelection=parameters["selectionPenalty"],
        baseTau=parameters["tau"],
        enablePiSelection=parameters["piRNASelection"],
        piPercentage = parameters["piPercentage"],
        numberOfPiRNA = parameters["numberOfPiRNA"],
        disablePiRecombination = parameters["disablePiRecombination"]
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
        regulationStr=parameters["regulationStr"],
        SingleFamily=True,
    )

    copyNumber = result["AvgCopyNum"]
    piOccupancy = result["TECoreOverlap"]
    resultToWrite = [copyNumber, piOccupancy]

    with open(parameters["FileName"], "w+") as f:
        for x in zip(*resultToWrite):
            f.write("{0}\t{1}\n".format(*x))

    return 0


commandArgs = defaultdict(list)
for k, v in ((k.lstrip("-"), v) for k, v in (a.split("=") for a in sys.argv[1:])):
    commandArgs[k].append(v)

for key, value in commandArgs.items():
    if key in parameters.keys():
        if type(parameters[key]) == type(True):
            parameters[key] = bool(strtobool(value[0]))
        else:
            parameters[key] = type(parameters[key])(value[0])
        value
        print("Supplied parameter: ", key)
        print("Supplied value: ", parameters[key])
        print("Paramter type: ", type(parameters[key]))

run = worker(parameters)
