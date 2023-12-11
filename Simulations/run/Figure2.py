# Init
import os
import sys
import json
import math
import copy

module_path = os.path.abspath(os.path.join("../src/"))
if module_path not in sys.path:
    sys.path.append(module_path)

# Imports
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


def makepar(pp):
    # pp[0] is HTgen, pp[1] is eta
    mypar = copy.deepcopy(parameters) # Otherwise the reference only is copied
    mypar["HGTgeneration"] = pp[0]
    mypar["eta"]           = pp[1]
    mypar["saveSuffix"]    = "HT" + str(pp[0]) + f"-eta{pp[1]:.3f}"	
    return mypar
    
outputdir = "../results/Results-fig2/"

with open('./Default.parameters', 'r') as file:
    parameters = json.load(file)

method             = "grid"   # Alternative: random / grid
maxHT              = 200
replicatesCond    = 400
replicatesEach    = 7

if method == "grid": 
    etas          = np.linspace(0, 1,     int(math.sqrt(replicatesCond))).tolist()
    HTgenerations = np.linspace(0, maxHT, int(math.sqrt(replicatesCond))).tolist()
    HTgenerations = [round(x) for x in HTgenerations]
    
    allpar = [makepar([h,e]) for h in HTgenerations for e in etas]
else:
    allpar = [makepar([round(random.uniform(0.0, maxHT)), random.uniform(0.0, 1.0)]) for i in range(replicatesCond)]



allparrep = sum(list(repeat(allpar,replicatesEach)), [])

with concurrent.futures.ProcessPoolExecutor(max_workers=parameters["maxProcceses"]) as executor:
    futures = [executor.submit(worker, arg, outputdir) for arg in allparrep]
    for future in concurrent.futures.as_completed(futures):
                print (future.result())
