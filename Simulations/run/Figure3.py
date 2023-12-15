# Init
import os
import sys
import json

module_path = os.path.abspath(os.path.join("../src/"))
if module_path not in sys.path:
    sys.path.append(module_path)

# Imports
import random
import numpy as np
import pandas as pd
import warnings
from numpy import concatenate as c
from itertools import repeat

# Simulation imports
from stats import stats
from regulation import regulation
from checkCopyNumber import checkCopyNumber
from fitness import calculateFitness
from transposition import transposition
from recombination import recombination

# Current multiprocessing implementation
from multiprocessing import Process
import concurrent.futures 


# Common functions for different figures
from common import worker

outputdir = "../results/Results-fig3/"

with open('./Default.parameters', 'r') as file:
    parameters = json.load(file)

etas          = np.linspace(0.0, 1.0, 11).tolist()
sels          = np.linspace(0.0, 0.01, 11).tolist()
HTgenerations = [0,200,parameters["Generations"]+1]
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
            futures = [executor.submit(worker, arg, outputdir) for arg in repeat(parameters,replicates)]
            for future in concurrent.futures.as_completed(futures):
                print (future.result())
    
    # Figure 2B
    for s in sels:
        parameters["HGTgeneration"] = HTgen
        parameters["eta"]           = 0.9
        parameters["selectionPenalty"] = -s # In the simulation program, selection coefficients are negative
        parameters["saveSuffix"]    = "HT" + str(HTgen) + f"-sel{s:.3f}"
        
        with concurrent.futures.ProcessPoolExecutor(max_workers=parameters["maxProcceses"]) as executor:
            futures = [executor.submit(worker, arg, outputdir) for arg in repeat(parameters,replicates)]
            for future in concurrent.futures.as_completed(futures):
                print (future.result())
