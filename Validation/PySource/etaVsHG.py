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

# Wrapper function for multiprocessing
def worker(u):
    genome, piset, piIndice, rates = generateGenome(
        numberOfInsertionSites=10000,
        numberOfChromosomes=6,
        baseRecombinationRate=0.1,
        baseSelection=0,
        baseTau=1,
    )
    population, transposons, TEset = generatePopulation(
        genome,
        piIndice,
        NumberOfIndividual=500,
        NumberOfTransposonTypes=1,
        NumberOfInsertionsPerType=[1],
        FrequencyOfInsertions=[1.0],
        ExcisionRates=[u],
        RepairRates=[1],
        InsertionRates=[1],
    )
    result = runSim(
        genomeMatrix=genome,
        populationMatrix=population,
        transposonMatrix=transposons,
        TEset=TEset,
        NumberOfTransposonInsertions=1,
        generations=5000,
        genMap=rates,
        piRNAindices=piIndice,
        simHGT=None,
        HMTgen=None,
        NumberOfInsertionsPerType=None,
        FrequencyOfInsertions=[1.0],
        ExcisionRates=[u],
        RepairRates=[1],
        InsertionRates=[1],
        eta=0,
        tau=1,
        selPen=0,
    )
    with open("R1_"+ '%030x' % random.randrange(16**30) + ".pickle", "wb") as f:
        pickle.dump((result), f)

    return 0

countSimulations = 1
u = 0.01
with concurrent.futures.ProcessPoolExecutor(max_workers=8) as executor:
    futures = [executor.submit(worker, arg) for arg in repeat(u,20)]
    for future in concurrent.futures.as_completed(futures):
        print (future)

u = 0.1
with concurrent.futures.ProcessPoolExecutor(max_workers=8) as executor:
    futures = [executor.submit(worker, arg) for arg in repeat(u,20)]
    for future in concurrent.futures.as_completed(futures):
        print (future)

u = 1.0
with concurrent.futures.ProcessPoolExecutor(max_workers=8) as executor:
    futures = [executor.submit(worker, arg) for arg in repeat(u,20)]
    for future in concurrent.futures.as_completed(futures):
        print (future)
