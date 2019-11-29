# Imports
import inspect
import random
import tempfile
import uuid
from functools import wraps
import collections

import coloredlogs
import dask as dd
import logging
import numpy as np
import pandas as pd
from coolname import generate_slug

# Logger
logger = logging.getLogger(__name__)
coloredlogs.install(level='DEBUG')

Genome = pd.DataFrame(
    columns=['Position', 'RecombinationRate', 'InsertionSite'])


def generateGenome(genomeSize, nInsertionSites, nChrom):
    """Method to generate genome
    
    :param genomeSize: Size of genome in cM
    :type genomeSize: int
    :param nInsertionSites: Number of insertion sites present in genome.
    :type nInsertionSites: int
    :param nChrom: Number of chromosomes present in genome
    :type nChrom: int
    :return: Dataframe containing genome information
    :rtype: DataFrame
    """
    genome = list(range(1, genomeSize + 1))
    InsertionSite = [0] * len(genome)
    RecombinationRate = [0.01] * len(genome)
    insertionLocation = np.random.choice(genomeSize,
                                         nInsertionSites,
                                         replace=False)
    chrLocation = np.random.choice(genomeSize, nChrom, replace=False)

    insertionCounter = 1
    for i in insertionLocation:
        InsertionSite[i] = insertionCounter
        insertionCounter += 1

    insertionCounter = 1
    for i in chrLocation:
        if (InsertionSite[i] != 0):
            RecombinationRate[i + 1] = 0.5
            insertionCounter += 1
        else:
            RecombinationRate[i] = 0.5
            insertionCounter += 1

    genomeDict = {
        'Position': genome,
        'RecombinationRate': RecombinationRate,
        'InsertionSite': InsertionSite
    }
    genome = pd.DataFrame(genomeDict)
    return (genome)

def recombination(genomeFrame, genomeOrg):
    """[summary]
    
    :param genomeFrame: DataFrame containing genome information
    :type genomeFrame: DataFrame
    :param genomeOrg: Pandas series containing one record of population DataFrame
    :type genomeOrg: Pandas series
    :return: List of transposons after recombination
    :rtype: int list
    """
    # Create insertion list for the progeny
    TEprogeny = []
    # Create a copy of genomeFrame
    genomeCopy = genomeFrame.copy(deep=True)
    Insertion_Father = genomeOrg['Insertion_Father']
    Insertion_Mother = genomeOrg['Insertion_Mother']
    if (Insertion_Father[0] == 0 and Insertion_Mother[0] == 0):
        return (random.choice(
            [genomeOrg['Insertion_Father'], genomeOrg['Insertion_Mother']]))
    else:
        initParent = random.choice(["M", "F"])
        surrogateParent = "M" if ("M" != initParent) else "F"
        RandArray = np.random.uniform(0, 1.0, genomeFrame.shape[0])
        switch = genomeCopy['RecombinationRate'] > RandArray
        #counter = collections.Counter(switch)
        #print(counter)
        genomeCopy['Progeny'] = np.where(switch.cumsum() % 2 == 0, initParent,
                                         surrogateParent)
        for i in Insertion_Father:
            se = genomeCopy[genomeCopy['InsertionSite'] ==
                            i]['Progeny'].values[0]
            if (se == "M"):
                TEprogeny.append(i)
        for i in Insertion_Mother:
            se = genomeCopy[genomeCopy['InsertionSite'] ==
                            i]['Progeny'].values[0]
            if (se == "F"):
                TEprogeny.append(i)
        #print(genomeCopy['Progeny'].value_counts())
    return (TEprogeny)