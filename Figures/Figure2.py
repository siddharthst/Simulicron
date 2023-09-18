# Init
import os
import sys
import json

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
    
    ff = "./Results-fig2/"+'%030x' % random.randrange(16**30) + "-" + parameters["saveSuffix"] + ".pickle"
    with open(ff, "wb") as f:
        pickle.dump((result), f)
        
    return ff





with open('./Default.parameters', 'r') as file:
    parameters = json.load(file)

etas          = np.linspace(0.0, 1.0, 11).tolist()
sels          = np.linspace(0.0, 0.01, 11).tolist()
HTgenerations = [0,300,parameters["Generations"]+1]
replicates    = 40


# ~ parameters["Generations"] = 10
# ~ parameters["Individuals"] = 5
# ~ parameters["NumberOfInsertionSites"] = 20

for HTgen in HTgenerations: 
    # Figure 2A
    for eta in etas:
        parameters["HGTgeneration"] = HTgen
        parameters["eta"]           = eta
        parameters["saveSuffix"]    = "HT" + str(HTgen) + f"-eta{eta:.2f}"
        
        with concurrent.futures.ProcessPoolExecutor(max_workers=parameters["maxProcceses"]) as executor:
            futures = [executor.submit(worker, arg) for arg in repeat(parameters,replicates)]
            for future in concurrent.futures.as_completed(futures):
                print (future.result())
    
    # Figure 2B
    for s in sels:
        parameters["HGTgeneration"] = HTgen
        parameters["eta"]           = 0.9
        parameters["selectionPenalty"] = -s # In the simulation program, selection coefficients are negative
        parameters["saveSuffix"]    = "HT" + str(HTgen) + f"-sel{s:.3f}"
        
        with concurrent.futures.ProcessPoolExecutor(max_workers=parameters["maxProcceses"]) as executor:
            futures = [executor.submit(worker, arg) for arg in repeat(parameters,replicates)]
            for future in concurrent.futures.as_completed(futures):
                print (future.result())
