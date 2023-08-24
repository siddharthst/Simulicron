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
import copy

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
    genome, piset, piIndice, rates = generateGenome(
        numberOfInsertionSites = parameters["NumberOfInsertionSites"],
        numberOfChromosomes    = parameters["NumberOfChromosomes"],
        baseRecombinationRate  = parameters["RecombinationRate"],
        baseSelection          = parameters["selectionPenalty"],
        baseTau                = parameters["tau"],
        disablePiRecombination = True,
    )
        
    population, transposons, TEset = generatePopulation(
        genome,
        piIndice,
        NumberOfIndividual     = parameters["Individuals"],
        NumberOfTransposonTypes = 2,
        NumberOfInsertionsPerType = [1, 0],
        FrequencyOfInsertions  =[
            parameters["FrequencyOfInsertionMain"],
            parameters["FrequencyOfInsertionHGT"],
        ],
        ExcisionRates          = [
            parameters["ExcisionRateMain"], 
            parameters["ExcisionRateHGT"],
        ],
        RepairRates            = [1, 1],
        InsertionRates         = [1, 1],
    )
    
    # run the simulation
    result = runSim(
        genomeMatrix           = genome,
        populationMatrix       = population,
        transposonMatrix       = transposons,
        TEset                  = TEset,
        NumberOfTransposonInsertions=2,
        generations            = parameters["Generations"],
        genMap                 = rates,
        piRNAindices           = piIndice,
        simHGT                 = None,
        HMTgen                 = parameters["HGTgeneration"],
        NumberOfInsertionsPerType = None,
        FrequencyOfInsertions  =[
            parameters["FrequencyOfInsertionMain"],
            parameters["FrequencyOfInsertionHGT"],
        ],
        ExcisionRates          = [
            parameters["ExcisionRateMain"], 
            parameters["ExcisionRateHGT"],
        ],
        RepairRates            = [1, 1],
        InsertionRates         = [1, 1],
        eta                    = parameters["eta"],
        tau                    = parameters["tau"],
        selPen                 = parameters["selectionPenalty"],
        epistasisCoefficient   = parameters["epistasisCoefficient"],
        fitnessFunction        = 2
    )
    
    ff = "./Results/"+'%030x' % random.randrange(16**30) + "-" + parameters["saveSuffix"] + ".pickle"
    with open(ff, "wb") as f:
        pickle.dump((result), f)
        
    return ff


with open('../../Default.parameters', 'r') as file:
    parameters = json.load(file)


parameters["Generations"] = 100
parameters["Individuals"] = 100
parameters["NumberOfInsertionSites"] = 20

etas          = np.arange(0, 1.0001, 0.05).tolist()
HTgenerations = np.arange(0, 200.01, 10).tolist()

replicates    = 3

def makepar(pp):
    # pp[0] is HTgen, pp[1] is eta
    mypar = copy.deepcopy(parameters) # Otherwise the reference only is copied
    mypar["HGTgeneration"] = pp[0]
    mypar["eta"]           = pp[1]
    mypar["saveSuffix"]    = "HT" + str(pp[0]) + f"-eta{pp[1]:.2f}"	
    return mypar
	
allpar = [makepar([h,e]) for h in HTgenerations for e in etas]
allparrep = sum(list(repeat(allpar,replicates)), [])

with concurrent.futures.ProcessPoolExecutor(max_workers=parameters["maxProcceses"]) as executor:
    futures = [executor.submit(worker, arg) for arg in allparrep]
    for future in concurrent.futures.as_completed(futures):
                print (future.result())



